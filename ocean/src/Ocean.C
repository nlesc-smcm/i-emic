//=====================================================================
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
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
	_SUBROUTINE_(write_data)(double*, int*, int*);  // file inout.f
} //extern
//=====================================================================

Ocean::Ocean(RCP<Epetra_Comm> Comm)
	:
	comm_(Comm),
	solverInitialized_(false)
{
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");

	DEBUG("Entering Ocean constructor...");   

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
	INFO("  Obtained solution from THCM");

	// Randomize state vector 
	double randScale = 1.0e-6;
	randomizeState(0.0);
	INFO("  Randomized solution and scaled with a factor " << randScale);
	
	// Obtain Jacobian from THCM    
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	jac_ = THCM::Instance().getJacobian();
	INFO("  Obtained jacobian from THCM");

	// Initialize a few datamembers
	sol_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	rhs_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));

	DEBUG("Leaving Ocean constructor...");
	
}
//=====================================================================
void Ocean::randomizeState(double scaling)
{
	DEBUG("Entering Ocean::randomizeState()...");
	state_->Random();
	state_->Scale(scaling);
	INFO("Initialized solution vector");
	DEBUG("Leaving  Ocean::randomizeState()...");
}
//=====================================================================
void Ocean::dumpState()
{
	// This function will probably break (?) on a distributed memory system.
	// For now it is convenient.
	// Use some HYMLS functionality to gather the solution in the right way
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
void Ocean::initializeSolver()
{
	DEBUG("Entering Ocean::initializeSolver()...");

	// Belos::LinearProblem setup

	problem_ =
		rcp(new Belos::LinearProblem
			<double, Epetra_MultiVector,Epetra_Operator>
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
   //belosParamList_->set("Explicit Residual Scaling", "Norm of Preconditioned Initial Residual");

	// Belos block GMRES setup
	belosSolver_ =
		rcp(new Belos::BlockGmresSolMgr
			<double, Epetra_MultiVector, Epetra_Operator>
			(problem_, belosParamList_));
	solverInitialized_ = true;

	DEBUG("Leaving Ocean::initializeSolver()...");
}
//=====================================================================
void Ocean::solve()
{
	DEBUG("Entering Ocean::solve()");
	if (!solverInitialized_)
		initializeSolver();

	sol_->PutScalar(0.0);
	// scaleProblem();
	precPtr_->Compute();
	bool set = problem_->setProblem(sol_, rhs_);
	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
							   "*** Belos::LinearProblem failed to setup");
	Belos::ReturnType ret = belosSolver_->solve();
	// unscaleProblem();
		
	double nrm;
	sol_->Norm2(&nrm);
	DEBUG(" Ocean::solve()   norm solution: " << nrm);
	DEBUG("Leaving Ocean::solve()");
}
 //=====================================================================
void Ocean::scaleProblem()
{
	DEBUG("Entering Ocean::scaleProblem()");
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

	DEBUG("Leaving Ocean::scaleProblem()");
}
 //=====================================================================
void Ocean::unscaleProblem()
{
	DEBUG("Entering Ocean::unscaleProblem()");
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

	DEBUG("Leaving Ocean::unscaleProblem()");
}
//=====================================================================
void Ocean::computeRHS()
{
	DEBUG("Entering Ocean::computeRHS()");
	// evaluate rhs in THCM with the current state
	THCM::Instance().evaluate(*state_, rhs_, false);
	rhs_->Scale(-1.0);
	DEBUG("Leaving Ocean::computeRHS()");
}
//=====================================================================
void Ocean::computeJacobian()
{
	DEBUG("Entering Ocean::computeJacobian()...");
	// Compute the Jacobian in THCM using the current state
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	// Get the Jacobian from THCM
	jac_ = THCM::Instance().getJacobian();
	DUMP("jac_.txt", *jac_);
	DEBUG("Leaving  Ocean::computeJacobian()...");
}
//=====================================================================
double Ocean::getNormRHS()
{
	double nrm;
	rhs_->Norm2(&nrm);
	return nrm;
}
//=====================================================================
double Ocean::getNormState()
{
	double nrm;
	state_->Norm2(&nrm);
	return nrm;
}
//=====================================================================
// OceanTheta
//=====================================================================
OceanTheta::OceanTheta(Teuchos::RCP<Epetra_Comm> Comm)
	:
	Ocean(Comm),
	theta_(1.0),
	timestep_(1.0e-01)
{
	DEBUG("Entering OceanTheta constructor");

	// Initialize a few datamembers
	oldState_ = rcp(new Epetra_Vector(*state_));
	stateDot_ = rcp(new Epetra_Vector(*state_));
	oldRhs_   = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	DEBUG("Leaving OceanTheta constructor");
}
//=====================================================================
void OceanTheta::parkModel()
{
	DEBUG("Entering Ocean::parkModel()");
	INFO("Parking the model"); 
	*oldState_ = *state_;
	*oldRhs_   = *rhs_;
	DEBUG("Leaving Ocean::parkModel()");
}
//=====================================================================
void OceanTheta::computeRHS()
{
	DEBUG("Entering OceanTheta::computeRHS()");
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

	DEBUG("Leaving OceanTheta::computeRHS()");
}
//=====================================================================
void OceanTheta::computeJacobian()
{
	DEBUG("Entering OceanTheta::computeJacobian()...");

    // Check theta
	if (theta_ < 0 || theta_ > 1)
	{
		INFO("Incorrect theta: " << theta_);
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
	// DUMP("jac_.txt", *jac_);
	DEBUG("Leaving  OceanTheta::computeJacobian()...");
}
