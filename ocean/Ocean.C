//=====================================================================
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Time.h>

#include <EpetraExt_HDF5.h>
#include <EpetraExt_Exception.h>

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
	comm_                (Comm),   // Setting the communication object
	solverInitialized_   (false),  // Solver needs initialization
	precInitialized_     (false),  // Preconditioner needs initialization
	recompPreconditioner_(true),   // We need a preconditioner to start with
	idrSolver_           (*this),  // Initialize IDR solver with current object (ocean);
	inputFile_           (oceanParamList->get("Input file", "ocean.h5")),
	outputFile_          (oceanParamList->get("Output file", "ocean.h5")),
	loadState_           (oceanParamList->get("Load state", false)),
	saveState_           (oceanParamList->get("Save state", false)),

// continuation
	parName_             (oceanParamList->get("Continuation parameter",
											  "Combined Forcing"))
{
	INFO("Ocean: constructor...");	
	
	Teuchos::ParameterList &thcmList =
		oceanParamList->sublist("THCM");
	
	// Create THCM object
	//  THCM is implemented as a Singleton, which allows only a single
	//  instance at a time. The Ocean class can access THCM with a call
	//  to THCM::Instance()
	thcm_ = rcp(new THCM(thcmList, comm_));
	INFO("Ocean: continuation parameter value: " << getPar());
	
	// Obtain solution vector from THCM
	state_ = THCM::Instance().getSolution();
	INFO("Ocean: Solution obtained from THCM");
	
	// Get domain object and get the problem dimensions
	domain_ = THCM::Instance().GetDomain();
	N_ = domain_->GlobalN();
 	M_ = domain_->GlobalM();
	L_ = domain_->GlobalL();

	// If specified we load a pre-existing state and parameter (x,l)
	if (loadState_)
		loadStateFromFile(inputFile_);
	else 			// Initialize with trivial solution
		state_->PutScalar(0.0);

	// Now that we have a state and a parameter we can initialize more datamembers
	initializeOcean();
	
	INFO("Ocean: constructor... done");
}

//=====================================================================
// destructor
Ocean::~Ocean()
{
	INFO("Ocean destructor called..."); 
}

//=====================================================================
// initialize Ocean with trivial state
void Ocean::initializeOcean()
{
	// Obtain Jacobian from THCM
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	jac_ = THCM::Instance().getJacobian();
	INFO("Ocean: Obtained jacobian from THCM");
	
	// Initialize a few datamembers
	sol_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	rhs_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
}

//====================================================================
void Ocean::preProcess()
{
	// Enable computation of preconditioner
	recompPreconditioner_ = true;
	INFO("Ocean pre-processing: enabling computation of preconditioner.");
}

//====================================================================
void Ocean::postProcess()
{
	if (saveState_)
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
	
	precPtr_->Initialize(); // Initialize
	precPtr_->Compute();    // Compute

	precInitialized_ = true;
	INFO("Ocean: initialize preconditioner done...");
}	

//=====================================================================
void Ocean::initializeSolver()
{
	INFO("Ocean: initialize solver...");
	
	solverParams_ = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams_.ptr());
	
	// Get the requested solver type
	solverType_          = solverParams_->get("Ocean solver type", 'I');
	
	// Initialize the preconditioner
	if (!precInitialized_)
		initializePreconditioner();

	// Initialize the requested solver
	if (solverType_ == 'F')
		initializeBelos();
	else if (solverType_ == 'I')
		initializeIDR();
	
	solverInitialized_ = true;
	
	// Now that the solver and preconditioner are initialized we are allowed to
	// perform a solve.
	INFO("Ocean: initialize solver... done");
}

//====================================================================
void Ocean::initializeIDR()
{
	idrSolver_.setParameters(solverParams_);
	idrSolver_.setSolution(getSolution('V'));
	idrSolver_.setRHS(getRHS('V'));
}

//====================================================================
void Ocean::initializeBelos()
{
	// If preconditioner not initialized do it now
	if (!precInitialized_) initializePreconditioner();

	// Belos LinearProblem setup
	problem_ = rcp(new Belos::LinearProblem
				   <double, Epetra_MultiVector, Epetra_Operator>
				   (jac_, sol_, rhs_) );

	// Set right preconditioner for Belos solver
	RCP<Belos::EpetraPrecOp> belosPrec =
		rcp(new Belos::EpetraPrecOp(precPtr_));
	problem_->setRightPrec(belosPrec);
	
	// A few FGMRES parameters are made available in solver_params.xml:
	int gmresIters  = solverParams_->get("FGMRES iterations", 500);
	double gmresTol = solverParams_->get("FGMRES tolerance", 1e-8);
	
	int NumGlobalElements = state_->GlobalLength();
	int maxrestarts       = 0;
	int blocksize         = 1; // number of vectors in rhs
	int maxiters          = NumGlobalElements/blocksize - 1;

	// Create Belos parameterlist
	belosParamList_ = rcp(new Teuchos::ParameterList());
	belosParamList_->set("Block Size", blocksize);
	belosParamList_->set("Flexible Gmres", true);
	belosParamList_->set("Adaptive Block Size", true);
	belosParamList_->set("Num Blocks", gmresIters);
	belosParamList_->set("Maximum Restarts", maxrestarts);
	belosParamList_->set("Orthogonalization","DGKS");
	belosParamList_->set("Output Frequency", 100);
	belosParamList_->set("Verbosity", Belos::TimingDetails +
			     Belos::Errors +
			     Belos::Warnings +
			     Belos::StatusTestDetails );
	belosParamList_->set("Maximum Iterations", maxiters); 
	belosParamList_->set("Convergence Tolerance", gmresTol); 
	belosParamList_->set("Explicit Residual Test", false); 
	belosParamList_->set("Implicit Residual Scaling", "Norm of RHS");
	
	// Belos block FGMRES setup
	belosSolver_ =
		rcp(new Belos::BlockGmresSolMgr
			<double, Epetra_MultiVector, Epetra_Operator>
			(problem_, belosParamList_));

}
//=====================================================================
void Ocean::solve(VectorPtr rhs)
{
	// Check whether solver is initialized, if not perform the
	// initialization here
	if (!solverInitialized_)
		initializeSolver();
	
	// Depending on the number of iterations we might need to
	// recompute the preconditioner
	if (recompPreconditioner_)
		buildPreconditioner();
	
	// Use trivial initial solution
	sol_->Scale(0.0);
	
	// Set the problem, rhs may be given as an argument to solve().
	if (solverType_ == 'F')
	{
		bool set;
		if (rhs == Teuchos::null)
			set = problem_->setProblem(sol_, rhs_);
		else
			set = problem_->setProblem(sol_, rhs->getOceanVector());
		TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
								   "*** Belos::LinearProblem failed to setup");
	}
	else if (solverType_ == 'I')
	{
		// Setup IDR: for our IDRSolver member we supply a view of the
		//            solution and the RHS:
		idrSolver_.setSolution(getSolution('V'));
		if (rhs == Teuchos::null)
			idrSolver_.setRHS(getRHS('V'));
		else
			idrSolver_.setRHS(getVector('V',rhs->getOceanVector()));
	}	

	// ---------------------------------------------------------------------
	// Start solving J*x = F, where J = jac_, x = sol_ and F = rhs
	TIMER_START("Ocean: solve...");
	INFO("Ocean: solve...");
	int    iters;
	double tol;
	if (solverType_ == 'F')
	{
		try
		{
			belosSolver_->solve(); 	// Solve
		}
		catch (std::exception const &e)
		{
			INFO("Ocean: exception caught: " << e.what());
		}
	}
	else if (solverType_ == 'I')
	{
		idrSolver_.solve();
	}
	
	INFO("Ocean: solve... done");
	TIMER_STOP("Ocean: solve...");
	// ---------------------------------------------------------------------

	// Do some post-processing
	if (solverType_ == 'F')
	{
		iters = belosSolver_->getNumIters();
		tol   = belosSolver_->achievedTol();
		INFO("Ocean: FGMRES, i = " << iters << ", ||r|| = " << tol);
		TRACK_ITERATIONS("Ocean: FGMRES iterations...", iters);
	}
	else if (solverType_ == 'I')
	{
		iters = idrSolver_.getNumIters();
		INFO("Ocean: IDR, i = " << iters
			 << " residual = " << idrSolver_.explicitResNorm());
		TRACK_ITERATIONS("Ocean: IDR iterations...", iters);
	}
}

//=====================================================================
double Ocean::explicitResNorm(VectorPtr rhs)
{
	RCP<Epetra_Vector> Ax =
		rcp(new Epetra_Vector(*(domain_->GetSolveMap())));
	jac_->Apply(*sol_, *Ax);
	Ax->Update(1.0, *(rhs->getOceanVector()), -1.0);
	double nrm;
	Ax->Norm2(&nrm);
	return nrm;
}

//=====================================================================
void Ocean::scaleProblem(VectorPtr rhs)
{
	// Not sure if this is the right approach and/or implemented correctly.
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

	if (rhs == Teuchos::null)
		rhs_->Multiply(1.0, *rowScalingRecipr_, *rhs_, 0.0);
	else
		(rhs->getOceanVector())->Multiply(1.0, *rowScalingRecipr_, *rhs_, 0.0);
	
	jac_->RightScale(*colScalingRecipr_);
	sol_->ReciprocalMultiply(1.0, *colScalingRecipr_, *sol_, 0.0);

	//------------------------------------------------------
	sol_->Norm2(&nrm);
	DEBUG("Ocean::scaleProblem() sol (after scaling): "  << nrm);
}

//=====================================================================
void Ocean::unscaleProblem(VectorPtr rhs)
{
	RCP<Epetra_Vector> rowScaling = THCM::Instance().getRowScaling();
	RCP<Epetra_Vector> colScaling = THCM::Instance().getColScaling();

	//------------------------------------------------------
	double nrm;
	sol_->Norm2(&nrm);
	DEBUG("Ocean::unscaleProblem() sol (before unscaling): " << nrm);
	
	//------------------------------------------------------
	jac_->LeftScale(*rowScaling);
	
	if (rhs == Teuchos::null)
		rhs_->Multiply(1.0, *rowScaling, *rhs_, 0.0);
	else
		(rhs->getOceanVector())->Multiply(1.0, *rowScaling, *rhs_, 0.0);
	
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
	TIMER_START("Ocean: apply matrix...");
	RCP<Epetra_Vector> result =
		rcp(new Epetra_Vector(*(domain_->GetSolveMap())));
	jac_->Apply(*(v.getOceanVector()), *result);

	TIMER_STOP("Ocean: apply matrix...");
	return getVector('V', result);
}

//====================================================================
void Ocean::applyMatrix(SuperVector const &v, SuperVector &out)
{
	TIMER_START("Ocean: apply matrix...");
	jac_->Apply(*(v.getOceanVector()), *(out.getOceanVector()));
	TIMER_STOP("Ocean: apply matrix...");
}

//====================================================================
void Ocean::buildPreconditioner()
{
	TIMER_START("Ocean: build preconditioner...");
	INFO("Ocean: build preconditioner...");
	precPtr_->Compute();		
	INFO("Ocean: build preconditioner... done");
	TIMER_STOP("Ocean: build preconditioner...");
	recompPreconditioner_ = false;  // Disable subsequent recomputes
}

//====================================================================
void Ocean::applyPrecon(SuperVector const &v, SuperVector &out)
{
	if (!precInitialized_) // Initialize preconditioner
		initializePreconditioner();
	
	if (recompPreconditioner_) 	// Compute preconditioner
		buildPreconditioner();	

	TIMER_START("Ocean: apply preconditioning...");
	precPtr_->ApplyInverse(*(v.getOceanVector()), *(out.getOceanVector()));
	TIMER_STOP("Ocean: apply preconditioning...");
}

//====================================================================
Teuchos::RCP<SuperVector> Ocean::applyPrecon(SuperVector const &v)
{
	RCP<Epetra_Vector> result =
		rcp(new Epetra_Vector(*(domain_->GetSolveMap())));
	RCP<SuperVector> out = getVector('V', result);
	applyPrecon(v, *out);
	return out;
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
		
		// Copy state
		std::stringstream ss;
		ss << "ocean_state_par" << std::setprecision(4) << std::setfill('_')
		   << std::setw(2) << THCM::Instance().par2int(parName_) << "_"
		   << std::setw(6) << getPar(parName_);
		INFO("copying fort.3 to " << ss.str());
		std::ifstream src("fort.3", std::ios::binary);
		std::ofstream dst(ss.str(), std::ios::binary);
		dst << src.rdbuf();
	}
	delete [] solutionArray;
}

//=====================================================================
int Ocean::saveStateToFile(std::string const &filename)
{	
	INFO("Writing to " << filename);
	
	double nrm;
	state_->Norm2(&nrm);
	INFO("   state: ||x|| = " << nrm);

 	// Write state, map and continuation parameter
	EpetraExt::HDF5 HDF5(*comm_);
	HDF5.Create(filename);
	HDF5.Write("State", *state_);

	// Interface between HDF5 and the THCM parameters,
	// store all the (_NPAR_ = 30) THCM parameters in an HDF5 file.
	std::string parName;
	double parValue;
	for (int par = 1; par <= _NPAR_; ++par)
	{
		parName  = THCM::Instance().int2par(par);
		parValue = getPar(parName);
		INFO("   " << parName << " = " << parValue);
		HDF5.Write("Parameters", parName.c_str(), parValue);
	}
	INFO("Writing to " << filename << " done");
	return 0;
}

// =====================================================================
int Ocean::loadStateFromFile(std::string const &filename)
{

	// Check whether file exists
	INFO("Loading from " << filename);

	std::ifstream file(filename);
	if (!file)
	{
		WARNING("Can't open " << filename
				<< " continue with trivial state", __FILE__, __LINE__);
		state_->PutScalar(0.0);
		initializeOcean();
		return 1;
	}
	else file.close();

	// Obtain state vector from THCM and put in datamember
	state_ = THCM::Instance().getSolution();
	
	// Create HDF5 object 
	EpetraExt::HDF5 HDF5(*comm_);
	Epetra_MultiVector *readState;

	// Read state
	HDF5.Open(filename);
	HDF5.Read("State", readState);

	// Create importer
	// target map: thcm domain SolveMap
	// source map: state with linear map  as read by HDF5.Read
	Teuchos::RCP<Epetra_Import> lin2solve =
		Teuchos::rcp(new Epetra_Import(*(domain_->GetSolveMap()),
									   readState->Map() ));
	
	// Import state from HDF5 into state_ datamember
	state_->Import(*((*readState)(0)), *lin2solve, Insert);

	double nrm;
	state_->Norm2(&nrm);
	INFO("   state: ||x|| = " << nrm);
	
	// Interface between HDF5 and the THCM parameters,
	// put all the (_NPAR_ = 30) THCM parameters back in THCM.
	std::string parName;
	double parValue;
	for (int par = 1; par <= _NPAR_; ++par)
	{
		parName  = THCM::Instance().int2par(par);
		
		// Read continuation parameter and put it in THCM
		try
		{
			HDF5.Read("Parameters", parName.c_str(), parValue);
		}
		catch (EpetraExt::Exception &e)
		{
			e.Print();
			continue;
		}
		
		setPar(parName, parValue);
		INFO("   " << parName << " = " << parValue);
	}
	INFO("Loading from " << filename << " done");
	return 0;
}

//====================================================================
// Adjust locally defined continuation parameter
// This happens when Ocean is managed directly by Continuation
double Ocean::getPar()
{
	return getPar(parName_);
}

//====================================================================
double Ocean::getPar(std::string const &parName)
{
	// We only allow parameters that are available in THCM
	int parIdent = THCM::Instance().par2int(parName);
	if (parIdent > 0 && parIdent <= _NPAR_)
	{
		double thcmPar;
		FNAME(getparcs)(&parIdent, &thcmPar);
		return thcmPar;
	}
	else  // If parameter not available we return 0
		return 0;
}

//====================================================================
void Ocean::setPar(double value)
{
	setPar(parName_, value);
}

//====================================================================
void Ocean::setPar(std::string const &parName, double value)
{
	// We only allow parameters that are available in THCM
	int parIdent = THCM::Instance().par2int(parName);
	if (parIdent > 0 && parIdent <= _NPAR_)
		FNAME(setparcs)(&parIdent, &value);
}
