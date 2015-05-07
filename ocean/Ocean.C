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
// Fortran: get access to write_data function
extern "C"
{ 
	_SUBROUTINE_(write_data)(double*, int*, int*); 
} 

//=====================================================================
// Constructor:
Ocean::Ocean(RCP<Epetra_Comm> Comm)
	:
	comm_(Comm),                     // Setting the communication object
	solverInitialized_(false),       // Solver needs initialization
	recomputePreconditioner_(true),  // --> xml
	recomputeBound_(10),             // --> xml
	useScaling_(false)               // --> xml
{
	// Check if outFile is specified
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");

	// Setup THCM parameters:
	RCP<Teuchos::ParameterList> globalParamList =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("thcm_params.xml",
								globalParamList.ptr());
	Teuchos::ParameterList &thcmList =
		globalParamList->sublist("THCM");

	// Create THCM object
    thcm_ = rcp(new THCM(thcmList, comm_));

	// Obtain solution vector from THCM
	//  THCM is implemented as a Singleton, which allows only a single
	//  instance at a time. The Ocean class can access THCM with a call
	//  to THCM::Instance()
	state_ = THCM::Instance().getSolution();
	INFO("Ocean: Solution obtained from THCM");

	// Initialize solution
	state_->PutScalar(0.0);
	INFO("Ocean: Initialized solution -> Ocean::state = zeros...");
	
	// Obtain Jacobian from THCM
	// --> Not sure if this is necessary here
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

	// Setup block preconditioner parameters
	// --> xml files should have a better home
	Teuchos::RCP<Teuchos::ParameterList> solverParams =
		Teuchos::rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml",
								solverParams.ptr());	

	// Get the domain decomposition from THCM, needed for the preconditioner.
	RCP<TRIOS::Domain> domain = THCM::Instance().GetDomain();

	// Create and initialize block preconditioner
	precPtr_ = 	Teuchos::rcp(new TRIOS::BlockPreconditioner
							 (jac_, domain, *solverParams));
	INFO("Ocean: initializing preconditioner...");
	precPtr_->Initialize();
	INFO("Ocean: initializing preconditioner... done");

	// Set as right preconditioner for Belos solver
	RCP<Belos::EpetraPrecOp> belosPrec =
		rcp(new Belos::EpetraPrecOp(precPtr_));
	problem_->setRightPrec(belosPrec);
	
	// Belos parameter setup
	// --> xml
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
	// --> xml
	
	// Belos block GMRES setup
	belosSolver_ =
		rcp(new Belos::BlockGmresSolMgr
			<double, Epetra_MultiVector, Epetra_Operator>
			(problem_, belosParamList_));

	// Now the solver and preconditioner are initialized we are allowed to
	// perform a solve.
	solverInitialized_ = true;
}

//=====================================================================
void Ocean::Solve(RCP<Epetra_Vector> rhs)
{
	// Check whether solver is initialized, if not perform the
	// initialization here
	if (!solverInitialized_)
		InitializeSolver();
	
	// Set the initial solution in the solver to zeros.
	// --> this should be improved/changed.
	sol_->PutScalar(0.0);

	// If required we scale the problem here
	// --> not sure if this is the correct approach
	if (useScaling_)
		ScaleProblem();

	// Used to store timings
	double time;

	// Depending on the number of iterations we might need to
	// recompute the preconditioner
	if (recomputePreconditioner_)
	{
		INFO("Ocean: Computing preconditioner...");

		// Start timing
		// --> all this timing stuff should be put inside a macro
		timer_->ResetStartTime();

		// Compute preconditioner
		precPtr_->Compute();

		// End timing
		time = timer_->ElapsedTime();
		
		INFO("Ocean: Computing preconditioner... done("
			 << time << " seconds)" );

		// Disable subsequent recomputes
		recomputePreconditioner_ = false;
	}

	// Set the problem, rhs depends may be given as an argument Solve().
	bool set;
	if (rhs == Teuchos::null)
 		set = problem_->setProblem(sol_, rhs_);
	else
		set = problem_->setProblem(sol_, rhs);
	
	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
							   "*** Belos::LinearProblem failed to setup");

	// Start solving J*x = F, where J = jac_, x = sol_ and F = rhs_
	INFO("Ocean: Perform solve...");

	// Start Timing
	timer_->ResetStartTime();

	// Solve
	Belos::ReturnType ret = belosSolver_->solve();

	// Stop Timing
	time = timer_->ElapsedTime();

	//
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
	DEBUG("Ocean::scaleProblem() ----->  sol (before scaling): "
		  << nrm);
	
	//------------------------------------------------------
	jac_->LeftScale(*rowScalingRecipr_);
	rhs_->Multiply(1.0, *rowScalingRecipr_, *rhs_, 0.0);
	jac_->RightScale(*colScalingRecipr_);
	sol_->ReciprocalMultiply(1.0, *colScalingRecipr_, *sol_, 0.0);

	//------------------------------------------------------
	sol_->Norm2(&nrm);
	DEBUG("Ocean::scaleProblem() ----->  sol (after scaling): "
		  << nrm);

}

//=====================================================================
void Ocean::UnscaleProblem()
{
	RCP<Epetra_Vector> rowScaling = THCM::Instance().getRowScaling();
	RCP<Epetra_Vector> colScaling = THCM::Instance().getColScaling();

	//------------------------------------------------------
	double nrm;
	sol_->Norm2(&nrm);
	DEBUG("Ocean::unscaleProblem() ----->  sol (before unscaling): "
		  << nrm);
	
	//------------------------------------------------------
	jac_->LeftScale(*rowScaling);
	rhs_->Multiply(1.0, *rowScaling, *rhs_, 0.0);
	jac_->RightScale(*colScaling);
	sol_->ReciprocalMultiply(1.0, *colScaling, *sol_, 0.0);

	//------------------------------------------------------
	sol_->Norm2(&nrm);
	DEBUG("Ocean::unscaleProblem() ----->  sol (after unscaling): "
		  << nrm);
}

//=====================================================================
void Ocean::ComputeRHS()
{
	// evaluate rhs in THCM with the current state
	THCM::Instance().evaluate(*state_, rhs_, false);
}

//=====================================================================
void Ocean::ComputeJacobian()
{
	// Compute the Jacobian in THCM using the current state
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	// Get the Jacobian from THCM
	jac_ = THCM::Instance().getJacobian();
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

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::GetSolution(char mode)
{
	if (mode == 'C')
	{
		Teuchos::RCP<Epetra_Vector> copySol =
			Teuchos::rcp(new Epetra_Vector(*sol_));
		return copySol;
	}
	else if (mode == 'V')
	{
		return sol_;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::GetState(char mode)
{
	if (mode == 'C')
	{
		Teuchos::RCP<Epetra_Vector> copyState =
			Teuchos::rcp(new Epetra_Vector(*state_));
		return copyState;
	}
	else if (mode == 'V')
	{
		return state_;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::GetRHS(char mode)
{
	if (mode == 'C')
	{
		Teuchos::RCP<Epetra_Vector> copyRhs =
			Teuchos::rcp(new Epetra_Vector(*rhs_));
		return copyRhs;
	}
	else if (mode == 'V')
	{
		return rhs_;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//=====================================================================
// NOT DONE YET
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
// NOT DONE YET
void LoadStateFromFile(std::string const &name)
{
	
}
