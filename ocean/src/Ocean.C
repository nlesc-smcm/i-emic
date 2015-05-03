//=====================================================================
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Time.h>

#include <EpetraExt_BlockMapIn.h>
#include <EpetraExt_BlockMapOut.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_VectorOut.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>

#include <Ifpack_Preconditioner.h>

//=====================================================================
#include "Ocean.H"
#include "THCM.H"
#include "THCMdefs.H"
#include "TRIOS_Domain.H"
#include "TRIOS_BlockPreconditioner.H"
#include "GlobalDefinitions.H"

//=====================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//=====================================================================
// Fortran stuff:
extern "C"
{ 
	_SUBROUTINE_(write_data)(double*, int*, int*); 
} 

//=====================================================================
// Constructor:
Ocean::Ocean(RCP<Epetra_Comm> Comm)
	:
	comm_(Comm),
	solverInitialized_(false),
	recomputePreconditioner_(true),
	recomputeBound_(10),
	useScaling_(false)
{
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");


	// Setup THCM parameters:
	RCP<Teuchos::ParameterList> globalParamList =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("../ocean/parameters/thcm_params.xml",
								globalParamList.ptr());
	Teuchos::ParameterList &thcmList =
		globalParamList->sublist("THCM");
	DEBUG(*globalParamList);
	DEBUG(thcmList);

	// Create THCM object
    thcm_ = rcp(new THCM(thcmList, comm_));

	// Obtain solution vector from THCM
	//  THCM is implemented as a Singleton, which allows only a single
	//  instance at a time. The Ocean class can access THCM with a call
	//  to THCM::Instance()
	state_ = THCM::Instance().getSolution();
	INFO("Ocean: Obtained solution from THCM");

	state_->PutScalar(0.0);
	INFO("Ocean: Initialized solution -> Ocean::state = zeros...");
	
	// Obtain Jacobian from THCM    
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	jac_ = THCM::Instance().getJacobian();
	INFO("Ocean: Obtained jacobian from THCM");

	// Initialize a few datamembers
	sol_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	rhs_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));

	// Create timing object
	timer_ = rcp(new Epetra_Time(*Comm));
}

//=====================================================================
void Ocean::RandomizeState(double scaling)
{
	state_->Random();
	state_->Scale(scaling);
	INFO("Ocean: Initialized solution vector");
}

//=====================================================================
void Ocean::DumpState()
{
	// This function will probably break (?) on a distributed memory system.
	// For now it is convenient.
	// Use some of Jonas' utilities to gather the solution in the right way
	Teuchos::RCP<Epetra_MultiVector> fullSol = Utils::Gather(*state_, 0);
	int filename = 3;
	int label    = 2;	
	int length   = fullSol->GlobalLength();
	double *solutionArray = new double[length]; 
	if (comm_->MyPID() == 0)
	{
		std::cout << "Writing to fort." << filename
				  << " at label " << label << "." << std::endl;

		// Using operator() to access first vector in multivector
		(*fullSol)(0)->ExtractCopy(solutionArray); //  
		FNAME(write_data)(solutionArray, &filename, &label);
	}
	delete [] solutionArray;
}

//=====================================================================
void Ocean::InitializeSolver()
{

	// Belos::LinearProblem setup

	problem_ = rcp(new Belos::LinearProblem
				   <double, Epetra_MultiVector, Epetra_Operator>
				   (jac_, sol_, rhs_) );

	// Block preconditioner 
	Teuchos::RCP<Teuchos::ParameterList> solverParams =
		Teuchos::rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("../ocean/parameters/solver_params.xml",
								solverParams.ptr());	

	RCP<TRIOS::Domain> domain = THCM::Instance().GetDomain();

	precPtr_ =
		Teuchos::rcp(new TRIOS::BlockPreconditioner(jac_, domain,
													*solverParams));

	//
	INFO("Ocean: initializing preconditioner...");
	precPtr_->Initialize();
	INFO("Ocean: initializing preconditioner... done");
	
	RCP<Belos::EpetraPrecOp> belosPrec =
		rcp(new Belos::EpetraPrecOp(precPtr_));
	problem_->setRightPrec(belosPrec);
	
	// Belos parameter setup
   belosParamList_ = rcp(new Teuchos::ParameterList());
   belosParamList_->set("Block Size", 1);
   belosParamList_->set("Flexible Gmres", true);
   belosParamList_->set("Adaptive Block Size", false);
   belosParamList_->set("Num Blocks",1000);
   belosParamList_->set("Maximum Restarts",0);
   belosParamList_->set("Orthogonalization","DGKS");
   belosParamList_->set("Output Frequency",1);
   belosParamList_->set("Maximum Iterations", 1000);
   belosParamList_->set("Convergence Tolerance", 1.0e-3); 
   belosParamList_->set("Explicit Residual Test", false); 
   belosParamList_->set("Verbosity", Belos::FinalSummary);
   belosParamList_->set("Implicit Residual Scaling", "Norm of RHS");

	// Belos block GMRES setup
	belosSolver_ =
		rcp(new Belos::BlockGmresSolMgr
			<double, Epetra_MultiVector, Epetra_Operator>
			(problem_, belosParamList_));
	solverInitialized_ = true;

}

//=====================================================================
void Ocean::Solve()
{
	if (!solverInitialized_)
		InitializeSolver();
	
	// 
	sol_->PutScalar(0.0);

	if (useScaling_)
		ScaleProblem();

	double time;
	if (recomputePreconditioner_)
	{
		INFO("Ocean: Computing preconditioner...");
		timer_->ResetStartTime();
		precPtr_->Compute();
		time = timer_->ElapsedTime(); 
		INFO("Ocean: Computing preconditioner... done("
			 << time << " seconds)" );
		recomputePreconditioner_ = false;
	}
	
	bool set = problem_->setProblem(sol_, rhs_);
	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
							   "*** Belos::LinearProblem failed to setup");
	INFO("Ocean: Perform solve...");
	timer_->ResetStartTime();
	Belos::ReturnType ret = belosSolver_->solve();
	time = timer_->ElapsedTime(); 
	belosIters_ = belosSolver_->getNumIters();	
	INFO("Ocean: Perform solve... done ("
		 << belosIters_ << " iterations, "
		 << time        << " seconds)");

	if (belosIters_ > recomputeBound_)
	{
		INFO("Ocean: Number of iterations exceeds " << recomputeBound_);
		INFO("Ocean:   Enabling computation of preconditioner.");
		recomputePreconditioner_ = true;
	}
	
	if (useScaling_)
		UnscaleProblem();
	
	double nrm;
	sol_->Norm2(&nrm);
	DEBUG(" Ocean::solve()   norm solution: " << nrm);
}

//=====================================================================
void Ocean::ScaleProblem()
{
	RCP<Epetra_Vector> rowScaling = THCM::Instance().getRowScaling();
	RCP<Epetra_Vector> colScaling = THCM::Instance().getColScaling();

	//------------------------------------------------------
	if (rowScalingRecipr_ == Teuchos::null or
		!rowScaling->Map().SameAs(rowScalingRecipr_->Map()))
	{
		rowScalingRecipr_ =
			rcp(new Epetra_Vector(rowScaling->Map()));
	}
	*rowScalingRecipr_ = *rowScaling;
	rowScalingRecipr_->Reciprocal(*rowScaling);
	
	//------------------------------------------------------
	if (colScalingRecipr_ == Teuchos::null or
		!colScaling->Map().SameAs(colScalingRecipr_->Map()))
	{
		colScalingRecipr_ =
			rcp(new Epetra_Vector(colScaling->Map()));
	}
	*colScalingRecipr_ = *colScaling;
	colScalingRecipr_->Reciprocal(*colScaling);

	//------------------------------------------------------
	double nrm;
	sol_->Norm2(&nrm);
	DEBUG("Ocean::scaleProblem() ----->  sol (before scaling): " << nrm);
	
	//------------------------------------------------------
	jac_->LeftScale(*rowScalingRecipr_);
	rhs_->Multiply(1.0, *rowScalingRecipr_, *rhs_, 0.0);
	jac_->RightScale(*colScalingRecipr_);
	sol_->ReciprocalMultiply(1.0, *colScalingRecipr_, *sol_, 0.0);

	//------------------------------------------------------
	sol_->Norm2(&nrm);
	DEBUG("Ocean::scaleProblem() ----->  sol (after scaling): " << nrm);

}

//=====================================================================
void Ocean::UnscaleProblem()
{
	RCP<Epetra_Vector> rowScaling = THCM::Instance().getRowScaling();
	RCP<Epetra_Vector> colScaling = THCM::Instance().getColScaling();

	//------------------------------------------------------
	double nrm;
	sol_->Norm2(&nrm);
	DEBUG("Ocean::unscaleProblem() ----->  sol (before unscaling): " << nrm);
	
	//------------------------------------------------------
	jac_->LeftScale(*rowScaling);
	rhs_->Multiply(1.0, *rowScaling, *rhs_, 0.0);
	jac_->RightScale(*colScaling);
	sol_->ReciprocalMultiply(1.0, *colScaling, *sol_, 0.0);

	//------------------------------------------------------
	sol_->Norm2(&nrm);
	DEBUG("Ocean::unscaleProblem() ----->  sol (after unscaling): " << nrm);

}

//=====================================================================
void Ocean::ComputeRHS()
{
	// evaluate rhs in THCM with the current state
	THCM::Instance().evaluate(*state_, rhs_, false);
	rhs_->Scale(-1.0);
}

//=====================================================================
void Ocean::ComputeJacobian()
{
	// Compute the Jacobian in THCM using the current state
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	// Get the Jacobian from THCM
	jac_ = THCM::Instance().getJacobian();
	DUMP("jac_.txt", *jac_);
}

//=====================================================================
double Ocean::GetNormRHS()
{
	double nrm;
	rhs_->Norm2(&nrm);
	return nrm;
}

//=====================================================================
double Ocean::GetNormState()
{
	double nrm;
	state_->Norm2(&nrm);
	return nrm;
}

//=====================================================================
void Ocean::SaveStateToFile(std::string const &name)
{
	std::string vectorFile    = name + "vec";
	std::string mapFile       = name + "map";
	std::string parameterFile = name + "pars";

	INFO("Writing to " << vectorFile);
	INFO("           " << mapFile);
	INFO("           " << parameterFile);
	
	EpetraExt::VectorToMatrixMarketFile
		(vectorFile.c_str(), *state_);
	EpetraExt::BlockMapToMatrixMarketFile
		(mapFile.c_str(), state_->Map());
}

//=====================================================================
void LoadStateFromFile(std::string const &name)
{
	
}

//=====================================================================
// OceanTheta
//=====================================================================
OceanTheta::OceanTheta(Teuchos::RCP<Epetra_Comm> Comm)
	:
	Ocean(Comm),
	theta_(1.0),
	timestep_(1.0e-03)
{

	// Initialize a few datamembers
	oldState_ = rcp(new Epetra_Vector(*state_));
	stateDot_ = rcp(new Epetra_Vector(*state_));
	oldRhs_   = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
}

//=====================================================================
void OceanTheta::Store()
{
	DEBUG("Storing the model");
	
	if (oldState_ == Teuchos::null or
		!oldState_->Map().SameAs(state_->Map()))
	{
		oldState_ =
			rcp(new Epetra_Vector(state_->Map()));
	}
	
	*oldState_ = *state_;

	if (oldRhs_ == Teuchos::null or
		!oldRhs_->Map().SameAs(rhs_->Map()))
	{
		oldRhs_ =
			rcp(new Epetra_Vector(rhs_->Map()));
	}
	
	*oldRhs_   = *rhs_;
}

//=====================================================================
void OceanTheta::Restore()
{
	DEBUG("Restoring the model");

	if (state_ == Teuchos::null or
		!state_->Map().SameAs(oldState_->Map()))
	{
		state_ =
			rcp(new Epetra_Vector(oldState_->Map()));
	}
	
	*state_ = *oldState_;

	if (rhs_ == Teuchos::null or
		!rhs_->Map().SameAs(oldRhs_->Map()))
	{
		rhs_ =
			rcp(new Epetra_Vector(oldRhs_->Map()));
	}
	
	*rhs_   = *oldRhs_;
}

//=====================================================================
void OceanTheta::ComputeRHS()
{
	THCM::Instance().evaluate(*state_, rhs_, false);

    // Calculate mass matrix
	THCM::Instance().evaluateB();

    // Get the mass matrix from THCM 
	massMatrix_ = rcp(&THCM::Instance().DiagB(), false);

    // Calculate d/dt x = (xnew - xold)/dt
	stateDot_->Update(1.0 / timestep_, *state_, -1.0 / timestep_,
					  *oldState_, 0.0);
	
	// Obtain number of local elements
	int numMyElements     = stateDot_->Map().NumMyElements();

    // Get a list of the global element IDs owned by the calling proc
	int *myGlobalElements = stateDot_->Map().MyGlobalElements();

	// Scale xdot with the values in the mass matrix
	double value;
	for (int i = 0; i != numMyElements; ++i)
	{
		value = (*massMatrix_)[i] * (*stateDot_)[i];
		stateDot_->ReplaceGlobalValues(1, &value, myGlobalElements + i);
	}

    // The final theta timestepping rhs is given by
	// -1 * (B d/dt x + theta*F(x) + (theta-1) * F(x_old))
	rhs_->Update(theta_ - 1.0, *oldRhs_, theta_);
	rhs_->Update(1.0, *stateDot_, 1.0);
	rhs_->Scale(-1.0);
}

//=====================================================================
void OceanTheta::ComputeJacobian()
{

    // Check theta
	if (theta_ < 0 || theta_ > 1)
	{
		WARNING("Ocean: Incorrect theta: " << theta_,
				__FILE__, __LINE__);
	}	

    // First compute the Jacobian in THCM using the current state
	THCM::Instance().evaluate(*state_, Teuchos::null, true);

    // Get the plain Jacobian from THCM
	jac_ = THCM::Instance().getJacobian();

    // Scale it with theta
	jac_->Scale(theta_);

    // Get the mass matrix from THCM (which is actually just a
	//    vector with diagonal elements)
	// Wrap it in a non-owning RCP
	massMatrix_ = rcp(&THCM::Instance().DiagB(), false);

    // Get the number of local elements
	int numMyElements     =	jac_->Map().NumMyElements();

    // Get a list of the global element IDs owned by the calling proc
	int *myGlobalElements = jac_->Map().MyGlobalElements();

    // Add to the Jacobian the values B[i]/dt
	double value;
	for (int i = 0; i != numMyElements; ++i)
	{
		value = (*massMatrix_)[i] / timestep_;
		jac_->SumIntoGlobalValues(myGlobalElements[i], 1,
								  &value, myGlobalElements + i);
	}
	jac_->FillComplete();
}


