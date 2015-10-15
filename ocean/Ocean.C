//=====================================================================
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Time.h>

#include <EpetraExt_HDF5.h>

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
Ocean::Ocean(RCP<Epetra_Comm> Comm, RCP<Teuchos::ParameterList> oceanParamList)
	:
	comm_(Comm),                   // Setting the communication object
	solverInitialized_(false),     // Solver needs initialization
	precInitialized_(false),       // Preconditioner needs initialization
	adaptivePrecCompute_ (oceanParamList->get("Use adaptive preconditioner computation", true)),
	recomputePreconditioner_(true),  // We need a preconditioner to start with
	useScaling_          (oceanParamList->get("Use scaling", false)),
	gmresIters_          (oceanParamList->get("Iterations in FGMRES solver", 100)),
	gmresTol_            (oceanParamList->get("Tolerance in FGMRES solver", 1e-3)),
	recomputeBound_      (oceanParamList->get("Preconditioner recompute bound", 400)),
	inputFile_           (oceanParamList->get("Input file", "ocean.h5")),
	outputFile_          (oceanParamList->get("Output file", "ocean.h5")),
	useExistingState_    (oceanParamList->get("Use existing state", false)),
	parIdent_(19),    // Initialize continuation parameters
	parValue_(0),
	parStart_(0),
	parEnd_(1)
{
	INFO("Ocean: constructor...");	
	
	Teuchos::ParameterList &thcmList =
		oceanParamList->sublist("THCM");

	// Create THCM object
	//  THCM is implemented as a Singleton, which allows only a single
	//  instance at a time. The Ocean class can access THCM with a call
	//  to THCM::Instance()
	thcm_ = rcp(new THCM(thcmList, comm_));
	
	// Obtain solution vector from THCM
	state_ = THCM::Instance().getSolution();
	INFO("Ocean: Solution obtained from THCM");
	
	// Initialize solution
	state_->PutScalar(0.0);
	INFO("Ocean: Initialized solution -> Ocean::state = zeros...");

	// Get domain object and get the problem dimensions
	domain_ = THCM::Instance().GetDomain();
	N_ = domain_->GlobalN();
 	M_ = domain_->GlobalM();
	L_ = domain_->GlobalL();
	
	// If specified we load a pre-existing state and parameter (x,l)
	if (useExistingState_)
		loadStateFromFile(inputFile_);
		
	// Obtain Jacobian from THCM
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	jac_ = THCM::Instance().getJacobian();
	INFO("Ocean: Obtained jacobian from THCM");
	
	// Initialize a few datamembers
	sol_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	rhs_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	
	// Put the correct parameter value in THCM
	setPar(parValue_);

	// Initialize solver
	initializeSolver();
	
	INFO("Ocean: constructor... done");
}

//=====================================================================
// destructor
Ocean::~Ocean()
{
	INFO("Ocean destructor called..."); 
}

//=====================================================================
void Ocean::randomizeState(double scaling)
{
	state_->Random();
	state_->Scale(scaling);
	INFO("Ocean: Initialized solution vector");
}

//====================================================================
void Ocean::preProcess()
{
	// Enable computation of preconditioner
	recomputePreconditioner_ = true;
	INFO("Ocean pre-processing: enabling computation of preconditioner.");
}

//====================================================================
void Ocean::postProcess()
{
	saveStateToFile(outputFile_);
	writeFortFiles();
}

//=====================================================================
void Ocean::initializePreconditioner()
{
	INFO("Ocean: initialize preconditioner...");
	// Setup block preconditioner parameters
	// --> xml files should have a better home
	Teuchos::RCP<Teuchos::ParameterList> precParams =
		Teuchos::rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("ocean_preconditioner_params.xml",
								precParams.ptr());	

	// Create and initialize block preconditioner
	precPtr_ = 	Teuchos::rcp(new TRIOS::BlockPreconditioner
							 (jac_, domain_, *precParams));

	precPtr_->Initialize();
	precPtr_->Compute();

	precInitialized_ = true;
	INFO("Ocean: initialize preconditioner done...");
}	

//=====================================================================
void Ocean::initializeSolver()
{
	INFO("Ocean: initialize solver...");
	if (!precInitialized_)
		initializePreconditioner();

	// Belos::LinearProblem setup
	problem_ = rcp(new Belos::LinearProblem
				   <double, Epetra_MultiVector, Epetra_Operator>
				   (jac_, sol_, rhs_) );

	// Set right preconditioner for Belos solver
	RCP<Belos::EpetraPrecOp> belosPrec =
		rcp(new Belos::EpetraPrecOp(precPtr_));
	problem_->setRightPrec(belosPrec);
	
	// Belos parameter setup
	// --> xml
	int NumGlobalElements = state_->GlobalLength();
	int maxrestarts = 0;
	int blocksize   = 1; // number of vectors in rhs
	int maxiters    = NumGlobalElements/blocksize - 1;
	belosParamList_ = rcp(new Teuchos::ParameterList());
	belosParamList_->set("Block Size", blocksize);
	belosParamList_->set("Flexible Gmres", true);
	belosParamList_->set("Adaptive Block Size", true);
	belosParamList_->set("Num Blocks", gmresIters_);
	belosParamList_->set("Maximum Restarts", maxrestarts);
	belosParamList_->set("Orthogonalization","DGKS");
	belosParamList_->set("Output Frequency", 100);
	belosParamList_->set("Verbosity", Belos::TimingDetails +
			     Belos::Errors +
			     Belos::Warnings +
			     Belos::StatusTestDetails );
	belosParamList_->set("Maximum Iterations", maxiters); 
	belosParamList_->set("Convergence Tolerance", gmresTol_); 
	belosParamList_->set("Explicit Residual Test", false); 
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
	INFO("Ocean: initialize solver... done");
}

//=====================================================================
void Ocean::solve(VectorPtr rhs)
{
	// Check whether solver is initialized, if not perform the
	// initialization here
	if (!solverInitialized_)
		initializeSolver();
	
	// If required we scale the problem here
	if (useScaling_)
		scaleProblem();

	// Depending on the number of iterations we might need to
	// recompute the preconditioner
	if (recomputePreconditioner_)
	{
		// Compute preconditioner
		TIMER_START("Ocean: build preconditioner...");
		INFO("Ocean: build preconditioner...");
		precPtr_->Compute();		
		INFO("Ocean: build preconditioner... done");
		TIMER_STOP("Ocean: build preconditioner...");
		recomputePreconditioner_ = false;  // Disable subsequent recomputes
	}

	// Set the problem, rhs may be given as an argument to solve().
	bool set;
	if (rhs == Teuchos::null)
 		set = problem_->setProblem(sol_, rhs_);
	else
		set = problem_->setProblem(sol_, rhs->getOceanVector());
	
	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
				   "*** Belos::LinearProblem failed to setup");

	// ---------------------------------------------------------------------
	// Start solving J*x = F, where J = jac_, x = sol_ and F = rhs_
	TIMER_START("Ocean: solve...");
	INFO("Ocean: solve...");
	try
	{
		belosSolver_->solve(); 	// Solve
	}
	catch (std::exception const &e)
	{
		INFO("Ocean: exception caught: " << e.what());
	}
	INFO("Ocean: solve... done");
	TIMER_STOP("Ocean: solve...");
	// ---------------------------------------------------------------------

	// Do some post-processing
	int    belosIters = belosSolver_->getNumIters();
	double belosTol   = belosSolver_->achievedTol();
	INFO("Ocean: FGMRES, i = " << belosIters << ", ||r|| = " << belosTol);
	TRACK_ITERATIONS("Ocean: FGMRES iterations...", belosIters);

	// If the number of linear solver iterations exceeds a preset bound
	// we recompute the preconditioner
	if ((belosIters > recomputeBound_) && adaptivePrecCompute_)
	{
		INFO("Ocean: Number of iterations exceeds " << recomputeBound_);
		INFO("Ocean:   Enabling computation of preconditioner.");
		recomputePreconditioner_ = true;
	}

	// If specified, unscale the problem
	if (useScaling_)
		unscaleProblem();
}

//=====================================================================
void Ocean::scaleProblem()
{
	// Not sure if this is the right approach or implemented correctly.
	// Scaling is obtained from THCM and then applied to the problem.
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
Teuchos::RCP<SuperVector> Ocean::applyMatrix(SuperVector const &v)
{
	RCP<Epetra_Vector> result =
		rcp(new Epetra_Vector(*(domain_->GetSolveMap())));
	jac_->Apply(*(v.getOceanVector()), *result);
	return getVector('V', result);
}

//====================================================================
Teuchos::RCP<SuperVector> Ocean::applyPrecon(SuperVector const &v)
{
	if (!precInitialized_)
		initializePreconditioner();
	
	RCP<Epetra_Vector> result =
		rcp(new Epetra_Vector(*(domain_->GetSolveMap())));
	precPtr_->ApplyInverse(*(v.getOceanVector()), *result);
	return getVector('V', result);
}

//====================================================================
void Ocean::setAtmosphere(std::vector<double> const &atmos)
{
	TIMER_START("Ocean: set atmosphere...");
	// This is a job for THCM
	THCM::Instance().setAtmosphere(atmos);
	TIMER_STOP("Ocean: set atmosphere...");
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
Teuchos::RCP<Epetra_CrsMatrix> Ocean::getJacobian()
{
	// This is a job for THCM
	return THCM::Instance().getJacobian();
}

//====================================================================
std::shared_ptr<std::vector<int> > Ocean::getSurfaceTRows()
{
	std::shared_ptr<std::vector<int> > rows =
		std::make_shared<std::vector<int> >();
	for (int j = 0; j != M_; ++j)
		for (int i = 0; i != N_; ++i)
			rows->push_back(FIND_ROW2(_NUN_, N_, M_, L_,i,j,L_-1,TT));

	return rows;
}

//====================================================================
// Fill and return a copy of the surface temperature
std::shared_ptr<std::vector<double> > Ocean::getSurfaceT()
{
	TIMER_START("Ocean: get surface temperature...");
	
	// Get list of global rows corresponding to surface temperature
	std::shared_ptr<std::vector<int> > rows = getSurfaceTRows();
	
	// Get a restricted Epetra_Vector containing only the surface temp values
	Teuchos::RCP<Epetra_Vector> restricted = Utils::RestrictVector(*state_, *rows);
	
	// Gather restricted vector to all procs
	Teuchos::RCP<Epetra_MultiVector> gathered = Utils::AllGather(*restricted);

	// The local array should be allocated
	int numel = rows->size();
	std::shared_ptr<std::vector<double> > surfaceT =
		std::make_shared<std::vector<double> >(numel, 0.0);

	// Get the values
	gathered->ExtractCopy(&(*surfaceT)[0], numel);

	TIMER_STOP("Ocean: get surface temperature...");	
	return surfaceT;
}

//====================================================================
std::shared_ptr<std::vector<int> > Ocean::getLandMask()
{
	return THCM::Instance().getLandMask();
}

//=====================================================================
void Ocean::writeFortFiles()
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
void Ocean::saveStateToFile(std::string const &filename)
{
	TIMER_START("Ocean::saveStateToFile...");
	
	INFO("Writing to " << filename);

 	// Write state, map and continuation parameter
	EpetraExt::HDF5 HDF5(*comm_);
	HDF5.Create(filename);
	HDF5.Write("State", *state_);
	HDF5.Write("Continuation parameter", "Value", parValue_);
	
	TIMER_STOP("Ocean::saveStateToFile...");
}

//=====================================================================
void Ocean::loadStateFromFile(std::string const &filename)
{
	TIMER_START("Ocean::loadStateFromFile...");

	// Check whether file exists
	INFO("Loading from " << filename);

	std::ifstream file(filename);
	if (!file)
	{
		WARNING("Can't open " << filename
				<< " continue with trivial state", __FILE__, __LINE__);
		return;
	}
	else file.close();
	
	// Create HDF5 object
	EpetraExt::HDF5 HDF5(*comm_);
	Epetra_MultiVector *state;

	// Read state
	HDF5.Open(filename);
	HDF5.Read("State", state);

	// Create importer
	// target map: thcm domain SolveMap
	// source map: state with linear map  as read by HDF5.Read
	Teuchos::RCP<Epetra_Import> lin2solve =
		Teuchos::rcp(new Epetra_Import(*(domain_->GetSolveMap()),
									   state->Map() ));
	
	// Import state from HDF5 into state_ datamember
	state_->Import(*((*state)(0)), *lin2solve, Insert);
	
	// Read continuation parameter and put it in THCM
	HDF5.Read("Continuation parameter", "Value", parValue_);
	setPar(parValue_);
	
	TIMER_STOP("Ocean::loadStateFromFile...");
}

//====================================================================
double Ocean::getPar(char mode)
{
	if (mode == 'V')
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
	if (mode == 'D') return parEnd_;
	if (mode == 'S') return parStart_;
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return -1;
	}
}

//====================================================================
void Ocean::setPar(double value)
{
	parValue_ = value;
	FNAME(setparcs)(&parIdent_, &value);
}
