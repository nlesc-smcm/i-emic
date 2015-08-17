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
#include "SuperVector.H"
#include "THCM.H"
#include "THCMdefs.H"
#include "TRIOS_Domain.H"
#include "TRIOS_BlockPreconditioner.H"
#include "GlobalDefinitions.H"

//=====================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//=====================================================================
// Get access to a few THCM functions
extern "C" _SUBROUTINE_(write_data)(double*, int*, int*);
extern "C" _SUBROUTINE_(getparcs)(int*, double*);
extern "C" _SUBROUTINE_(setparcs)(int*,double*);
extern "C" _SUBROUTINE_(getooa)(double*);

//=====================================================================
// Constructor:
Ocean::Ocean(RCP<Epetra_Comm> Comm)
	:
	comm_(Comm),                     // Setting the communication object
	solverInitialized_(false)        // Solver needs initialization
{
	// Check if outFile is specified
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");
	
	// Setup Ocean and THCM parameters:
	RCP<Teuchos::ParameterList> oceanParamList =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("ocean_params.xml",
								oceanParamList.ptr());

	recomputePreconditioner_ =
		oceanParamList->get("recomputePreconditioner", true);
	recomputeBound_          =
		oceanParamList->get("recomputeBound", 50);
	useScaling_              =
		oceanParamList->get("useScaling", false);

	
	Teuchos::ParameterList &thcmList =
		oceanParamList->sublist("THCM");
	
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

	// Get domain object and set the problem dimensions
	domain_ = THCM::Instance().GetDomain();
	N_ = domain_->GlobalN();
 	M_ = domain_->GlobalM();
	L_ = domain_->GlobalL();
	
	// Initialize sst
	sst_ = std::make_shared<std::vector<double> >(N_*M_, 0.0);

	// Create fullSol_ array
	fullSol_ = new double[state_->GlobalLength()];

	// Set THCM continuation parameter and its bounds
	parIdent_ = 19;
	parValue_ = 0.0;
	parStart_ = 0.0;
	parEnd_   = 1.0;

	// Put the initial value into the model.
	setPar(parStart_);	
}

//=====================================================================
// destructor
Ocean::~Ocean()
{
	INFO("Ocean destructor called..."); 
	delete [] fullSol_;
}

//=====================================================================
void Ocean::randomizeState(double scaling)
{
	state_->Random();
	state_->Scale(scaling);
	INFO("Ocean: Initialized solution vector");
}

//=====================================================================
void Ocean::dumpState()
{
	// This function may break on a distributed memory system.
	Teuchos::RCP<Epetra_MultiVector> solution = Utils::Gather(*state_, 0);
	int filename = 3;
	int label    = 2;
	int length   = solution->GlobalLength();
	double *solutionArray = new double[length]; 
	if (comm_->MyPID() == 0)
	{
		std::cout << "Writing to fort." << filename
				  << " at label " << label << "." << std::endl;

		// Using operator() to access first EpetraVector in multivector
		(*solution)(0)->ExtractCopy(solutionArray); //  
		FNAME(write_data)(solutionArray, &filename, &label);
	}
	delete [] solutionArray;
}

//=====================================================================
void Ocean::initializeSolver()
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
	belosParamList_->set("Convergence Tolerance", 5.0e-3); 
	belosParamList_->set("Explicit Residual Test", false); 
	//belosParamList_->set("Verbosity", Belos::FinalSummary);
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
void Ocean::solve(RCP<SuperVector> rhs)
{
	// Check whether solver is initialized, if not perform the
	// initialization here
	if (!solverInitialized_)
		initializeSolver();
	
	// Set the initial solution in the solver to zeros.
	// --> this should be improved/changed but having zero as
	//     initial solution seems to be slightly faster
	sol_->PutScalar(0.0);

	// If required we scale the problem here
	// --> not sure if this is the correct approach
	if (useScaling_)
		scaleProblem();

	// Depending on the number of iterations we might need to
	// recompute the preconditioner
	if (recomputePreconditioner_)
	{
		// Compute preconditioner
		TIMER_START("Ocean: build preconditioner...");
		precPtr_->Compute();		
		TIMER_STOP("Ocean: build preconditioner...");
		recomputePreconditioner_ = false;  // Disable subsequent recomputes
	}

	// Set the problem, rhs may be given as an argument to Solve().
	bool set;
	if (rhs == Teuchos::null)
 		set = problem_->setProblem(sol_, rhs_);
	else
		set = problem_->setProblem(sol_, rhs->getOceanVector());
	
	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
							   "*** Belos::LinearProblem failed to setup");

	// Start solving J*x = F, where J = jac_, x = sol_ and F = rhs_
	TIMER_START("Ocean: solve...");
	Belos::ReturnType ret = belosSolver_->solve();	// Solve
	belosIters_ = belosSolver_->getNumIters();		
	INFO("Ocean: finished solve... " << belosIters_ << " iterations");
	TIMER_STOP("Ocean: solve...");

	if (belosIters_ > recomputeBound_)
	{
		INFO("Ocean: Number of iterations exceeds " << recomputeBound_);
		INFO("Ocean:   Enabling computation of preconditioner.");
		recomputePreconditioner_ = true;
	}
	
	if (useScaling_)
		unscaleProblem();
	
	double nrm;
	sol_->Norm2(&nrm);
	DEBUG(" Ocean::solve()   norm solution: " << nrm);
}

//=====================================================================
void Ocean::scaleProblem()
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
	DEBUG("Ocean::scaleProblem() sol (before scaling): "  << nrm);
	
	//------------------------------------------------------
	jac_->LeftScale(*rowScalingRecipr_);
	rhs_->Multiply(1.0, *rowScalingRecipr_, *rhs_, 0.0);
	jac_->RightScale(*colScalingRecipr_);
	sol_->ReciprocalMultiply(1.0, *colScalingRecipr_, *sol_, 0.0);

	//------------------------------------------------------
	sol_->Norm2(&nrm);
	DEBUG("Ocean::scaleProblem() sol (after scaling): "  << nrm);
}

//=====================================================================
void Ocean::unscaleProblem()
{
	RCP<Epetra_Vector> rowScaling = THCM::Instance().getRowScaling();
	RCP<Epetra_Vector> colScaling = THCM::Instance().getColScaling();

	//------------------------------------------------------
	double nrm;
	sol_->Norm2(&nrm);
	DEBUG("Ocean::unscaleProblem() sol (before unscaling): " << nrm);
	
	//------------------------------------------------------
	jac_->LeftScale(*rowScaling);
	rhs_->Multiply(1.0, *rowScaling, *rhs_, 0.0);
	jac_->RightScale(*colScaling);
	sol_->ReciprocalMultiply(1.0, *colScaling, *sol_, 0.0);

	//------------------------------------------------------
	sol_->Norm2(&nrm);
	DEBUG("Ocean::unscaleProblem() sol (after unscaling): " << nrm);
}

//=====================================================================
void Ocean::computeRHS()
{
	// evaluate rhs in THCM with the current state
 	TIMER_START("Ocean: compute RHS...");
	THCM::Instance().evaluate(*state_, rhs_, false);
	TIMER_STOP("Ocean: compute RHS...");
}

//=====================================================================
void Ocean::computeJacobian()
{
	TIMER_START("Ocean: compute Jacobian...");
	// Compute the Jacobian in THCM using the current state
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	// Get the Jacobian from THCM
	jac_ = THCM::Instance().getJacobian();
	TIMER_STOP("Ocean: compute Jacobian...");
}

//====================================================================
Teuchos::RCP<SuperVector> Ocean::getVector(char mode, RCP<Epetra_Vector> vec)
{
	if (mode == 'C') // copy
	{
		RCP<Epetra_Vector> copy = rcp(new Epetra_Vector(*vec));
		RCP<SuperVector> ptr    = rcp(new SuperVector(copy));
		return ptr;
	}
	else if (mode == 'V') // view
	{
		RCP<SuperVector> ptr = rcp(new SuperVector(vec));
		return ptr;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//====================================================================
Teuchos::RCP<SuperVector> Ocean::getSolution(char mode)
{
	return getVector(mode, sol_);
}

//====================================================================
Teuchos::RCP<SuperVector> Ocean::getState(char mode)
{
	return getVector(mode, state_);
}

//====================================================================
Teuchos::RCP<SuperVector> Ocean::getRHS(char mode)
{
	return getVector(mode, rhs_);
}

//====================================================================
void Ocean::setAtmosphere(std::vector<double> const &atmos)
{
	// This is a job for THCM
	THCM::Instance().setAtmosphere(atmos);
}

//====================================================================
std::shared_ptr<std::vector<double> > Ocean::getAtmosBlock()
{
	double Ooa;
	FNAME(getooa)(&Ooa);
	std::shared_ptr<std::vector<double> > values =
		std::make_shared<std::vector<double> >(N_ * M_, -Ooa);
	return values;
}

//====================================================================
std::shared_ptr<std::vector<int> > Ocean::getSSTRows()
{
	std::shared_ptr<std::vector<int> > rows =
		std::make_shared<std::vector<int> >();
	for (int j = 0; j != M_; ++j)
		for (int i = 0; i != N_; ++i)
			rows->push_back(FIND_ROW2(_NUN_, N_, M_, L_,i,j,L_-1,TT));

	return rows;
}

//====================================================================
// Fill and return a copy of the SST
std::shared_ptr<std::vector<double> > Ocean::getSST()
{
	// extract solution from gathered vector --> this is slow
	// everything should be better distributed
	Teuchos::RCP<Epetra_MultiVector> gathered =
		Utils::AllGather(*state_);
	gathered->ExtractCopy(fullSol_, state_->GlobalLength());
	int row;
	int idx = 0;
	
	// sst_ should be initialized
 	for (int j = 0; j != M_; ++j)
		for (int i = 0; i != N_; ++i)
		{
			// Using macros from THCMdefs.H
			row = FIND_ROW2(_NUN_, N_, M_, L_, i, j, L_-1, TT);
			(*sst_)[idx] = fullSol_[row];
			++idx;
		}
	return sst_;
}

//=====================================================================
// NOT IMPLEMENTED YET
void Ocean::saveStateToFile(std::string const &name)
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
// NOT IMPLEMENTED YET
void Ocean::loadStateFromFile(std::string const &name)
{}

//====================================================================
double Ocean::getPar()
{
	double thcmPar;
	FNAME(getparcs)(&parIdent_, &thcmPar);
	if (thcmPar != parValue_)
	{
		INFO("Ocean::getPar(): Faulty parameter synchronization");
		INFO("              thcm: " << thcmPar);
		INFO("             Ocean: " << parValue_);
		INFO("     fixing this...");
		FNAME(setparcs)(&parIdent_, &parValue_);
		FNAME(getparcs)(&parIdent_, &thcmPar);
		INFO("              thcm: " << thcmPar);
	}					  
	return parValue_;
}

//====================================================================
void Ocean::setPar(double value)
{
	parValue_ = value;
	FNAME(setparcs)(&parIdent_, &value);
}
