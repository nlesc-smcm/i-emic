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
#include "THCM.H"
#include "THCMdefs.H"
#include "TRIOS_Domain.H"
#include "TRIOS_BlockPreconditioner.H"
#include "GlobalDefinitions.H"

//==================================================================
#include <math.h>

//=====================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//=====================================================================
// Get access to a few THCM functions
extern "C" _SUBROUTINE_(write_data)(double*, int*, int*);
extern "C" _SUBROUTINE_(getparcs)(int*, double*);
extern "C" _SUBROUTINE_(setparcs)(int*,double*);
extern "C" _SUBROUTINE_(getooa)(double*, double*);

//=====================================================================
// Constructor:
Ocean::Ocean(RCP<Epetra_Comm> Comm, RCP<Teuchos::ParameterList> oceanParamList)
	:
	comm_                (Comm),   // Setting the communication object
	solverInitialized_   (false),  // Solver needs initialization
	precInitialized_     (false),  // Preconditioner needs initialization
	recompPreconditioner_(true),   // We need a preconditioner to start with
	idrSolver_           (*this),  // Initialize IDR solver with current object (ocean);
	inputFile_           (oceanParamList->get("Input file",  "ocean_input.h5")),
	outputFile_          (oceanParamList->get("Output file", "ocean_output.h5")),
	loadState_           (oceanParamList->get("Load state", false)),
	saveState_           (oceanParamList->get("Save state", false)),
	saveMask_            (oceanParamList->get("Save mask", true)),
	storeEverything_     (oceanParamList->get("Store everything", false)),

	parName_             (oceanParamList->get("Continuation parameter",
											  "Combined Forcing")),

	landmaskFile_        (oceanParamList->sublist("THCM").get("Land Mask", "none"))
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
	
	// Get domain object and get the problem dimensions
	domain_ = THCM::Instance().GetDomain();
	N_ = domain_->GlobalN();
 	M_ = domain_->GlobalM();
	L_ = domain_->GlobalL();
	
	// Read starting parameters from xml
	Teuchos::ParameterList& startList =
		thcmList.sublist("Starting Parameters");	
	THCM::Instance().ReadParameters(startList);

	// If specified we load a pre-existing state and parameters (x,l)
	// thereby overwriting the starting parameters
	if (loadState_)
		loadStateFromFile(inputFile_);
	else 			// Initialize with trivial solution
		state_->PutScalar(0.0);
	
	// Now that we have the state and parameters initialize
	// the Jacobian, solution and rhs
	initializeOcean();

	// Analyze Jacobian and print/fix impossible land points
	while( analyzeJacobian() )
	{
		// If we find singular pressure rows we adjust the current landmask
		INFO(" Fixing landmask " << landmaskFile_);
		
		Teuchos::RCP<Epetra_IntVector> landmask =
			THCM::Instance().getLandMask("current", singRows_);
		
		INFO(" Putting a fixed version of " << landmaskFile_ <<
			 "  back in THCM...");
		
		THCM::Instance().setLandMask(landmask);
		THCM::Instance().evaluate(*state_, Teuchos::null, true);		
	}	

	// Initialize preconditioner
	initializePreconditioner();

	// Get current global masks to communicate with an Atmosphere
	landmask_ = THCM::Instance().getLandMask();
	surfmask_ = THCM::Instance().getSurfaceMask();

	// Inspect current state
	inspectVector(state_);

	//------------------------------------------------------------------
	// Create surface temperature restrict/import strategy
	//------------------------------------------------------------------
	// Create list of surface temp row indices
	std::vector<int> tRows;

	for (int j = 0; j != M_; ++j)
		for (int i = 0; i != N_; ++i)
			tRows.push_back(FIND_ROW2(_NUN_, N_, M_, L_,i,j,L_-1,TT));

	// Create restricted map
	Teuchos::RCP<Epetra_BlockMap> indexMap =
		Utils::CreateSubMap(state_->Map(), tRows);

	// Create the SST vector
	sst_ = Teuchos::rcp(new Epetra_Vector(*indexMap));

	// Create importer
	// Target map: indexMap
	// Source map: state_->Map()
	Teuchos::RCP<Epetra_Import> surfaceTimporter_ =
		Teuchos::rcp(new Epetra_Import(*indexMap, state_->Map()));
	
	INFO("Ocean: constructor... done");
}

//=====================================================================
// destructor
Ocean::~Ocean()
{
	INFO("Ocean destructor"); 
}

//=====================================================================
// initialize Ocean with trivial state
void Ocean::initializeOcean()
{
	// Obtain Jacobian from THCM
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	jac_ = THCM::Instance().getJacobian();
	INFO("Ocean: Obtained Jacobian from THCM");

	// Initialize solution and rhs
	sol_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	rhs_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));

	// Get the rowmap for the pressure points
	Teuchos::RCP<Epetra_Map> RowMap = domain_->GetSolveMap();
	
	mapP_     = Utils::CreateSubMap(*RowMap, _NUN_, PP);
	mapU_     = Utils::CreateSubMap(*RowMap, _NUN_, UU);
		
	singRows_ = Teuchos::rcp(new Epetra_Vector(*mapP_));
}

//====================================================================
int Ocean::analyzeJacobian()
{
	// Make preparations for extracting pressure rows
	int dim = mapP_->NumMyElements();

	assert(dim == mapU_->NumMyElements());
	
	int row, len;
	int maxlen = jac_->MaxNumEntries();
	std::vector<int> indices(maxlen);
	std::vector<double> values(maxlen);

	int singRowsFound = 0;
	double sum;
	int el;
	int pos;

	for (int i = 0; i != dim; i++)
	{
		// Check pressure rows for weird things
		row = mapP_->GID(i);
		CHECK_ZERO(jac_->ExtractGlobalRowCopy(row, maxlen,
											  len, &values[0],
											  &indices[0]));
		sum = 0.0;
		el  = 0;
		for (int p = 0; p != len; p++)
		{
			sum += values[p];
			if (std::abs(values[p]) > 1e-10)
			{
				el++;
				if (values[p] > 0)
					pos++;
				else
					pos = 0;
			}
		}
		
		if (sum == 1)
		{
			// If the sum of the elements equals 1 this is already a land cell
			// We don't need to 'fix' it
			(*singRows_)[i] = 1;
		}
		else if (el <= 2)
		{
			// If the number of contributions is less than or equal to 2,
			// it is likely to be one of our 'problem rows'.
			// These rows only depend on vertical velocity W.
			(*singRows_)[i] = 2;
			singRowsFound++;
		}
	}
	
	int maxFound;
	comm_->MaxAll(&singRowsFound, &maxFound, 1);

	INFO("  analyzeJacobian(): problem rows found: " << maxFound);
	
	return maxFound;
}

//==================================================================
void Ocean::inspectVector(Teuchos::RCP<Epetra_Vector> x)
{
	INFO("Ocean: inspect vector...");

	// for now we just check whether the surface w-values are zero
	// if not we put them to zero
	
	// this setup lets us choose more rows in a layer to reset
	// so now we are going to put the surface u,v,p anomalies to zero as well
	// this is evil
	
	std::vector<int> rows;	
	int lid = 0;
	std::vector<int> badRows;
	int flag = 0;
	int globFlag = 0;
	
	for (int j = 0; j != M_; ++j)
		for (int i = 0; i != N_; ++i)
		{
			flag = 0;
			globFlag = 0;
			
			// surface w row
			rows.push_back((L_-1)*M_*N_*_NUN_ + j*N_*_NUN_ + i*_NUN_ + 2);

			// // surface u row
			// rows.push_back((L_-1)*M_*N_*_NUN_ + j*N_*_NUN_ + i*_NUN_ + 0);

			// // surface v row
			// rows.push_back((L_-1)*M_*N_*_NUN_ + j*N_*_NUN_ + i*_NUN_ + 1);

			// // surface p row
			// rows.push_back((L_-1)*M_*N_*_NUN_ + j*N_*_NUN_ + i*_NUN_ + 3);

			for (auto &row: rows)
			{
				lid = x->Map().LID(row);
				if (lid >= 0)
				{
					if (std::abs((*x)[lid]) > 1e-10)
					{
						(*x)[lid] = 0.0;
						flag = true;
					}
				}
				comm_->SumAll(&flag, &globFlag, 1);
				
				if (globFlag)
					badRows.push_back(row);
			}
			rows.clear();
		}
	if (badRows.size() > 0)
	{
		INFO("   fixed bad w points in surface rows: ");
		for (auto &el: badRows)
			INFO(el);		
	}
	INFO("Ocean: inspect vector... done");
}

//====================================================================
Ocean::LandMask Ocean::getLandMask()
{
	return getLandMask("current");
}	

//====================================================================
Ocean::LandMask Ocean::getLandMask(std::string const &fname)
{
	LandMask mask;
	
	mask.local = THCM::Instance().getLandMask(fname);
	THCM::Instance().setLandMask(mask.local);
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	
	// zero the singRows array
	singRows_->PutScalar(0.0);
	
	while( analyzeJacobian() > 0 )
	{
		// If we find singular pressure rows we adjust the current landmask
		mask.local = THCM::Instance().getLandMask("current", singRows_);
		
		//  Putting a fixed version of the landmask back in THCM
		THCM::Instance().setLandMask(mask.local);
		THCM::Instance().evaluate(*state_, Teuchos::null, true);
	}
	
	// Get the current global landmask from THCM.
	mask.global = THCM::Instance().getLandMask();
	mask.label = fname;
	return mask;
}

//===================================================================
void Ocean::setLandMask(LandMask mask, bool global)
{
	INFO("Ocean: set landmask " << mask.label << "...");
	THCM::Instance().setLandMask(mask.local);
	
	if (global)
		THCM::Instance().setLandMask(mask.global);
	
	currentMask_ = mask.label;
	INFO("Ocean: set landmask " << mask.label << "... done");
}

//==================================================================
void Ocean::applyLandMask(LandMask mask, double factor)
{
	
	int nmask = mask.global->size();
	assert(nmask == (M_+2)*(N_+2)*(L_+2));
	
	int idx1, idx2;
	int lid;
	int ii;

	for (int k = 1; k != L_+1; ++k)
		for (int j = 1; j != M_+1; ++j)
			for (int i = 1; i != N_+1; ++i)
			{
				idx1 = k*(M_+2)*(N_+2) + j*(N_+2) + i;
				idx2 = (k-1)*M_*N_ + (j-1)*N_ + i-1;

				if ((*mask.global)[idx1])  // On land
				{
					for (ii = idx2*_NUN_; ii != (idx2+1)*_NUN_; ++ii)
					{
						lid = state_->Map().LID(ii); // Local ID of this point			
						if (lid >= 0)
							(*state_)[lid] *= factor;
					}
				}
			}
}

//==================================================================
void Ocean::applyLandMask(Teuchos::RCP<Epetra_Vector> x,
						  LandMask mask, double factor)
{
	int nmask = mask.global->size();
	assert(nmask == (M_+2)*(N_+2)*(L_+2));
	
	int idx1, idx2;
	int lid;
	int ii;

	for (int k = 1; k != L_+1; ++k)
		for (int j = 1; j != M_+1; ++j)
			for (int i = 1; i != N_+1; ++i)
			{
				idx1 = k*(M_+2)*(N_+2) + j*(N_+2) + i;
				idx2 = (k-1)*M_*N_ + (j-1)*N_ + i-1;

				if ((*mask.global)[idx1])  // On land
				{
					for (ii = idx2*_NUN_; ii != (idx2+1)*_NUN_; ++ii)
					{
						lid = x->Map().LID(ii); // Local ID of this point			
						if (lid >= 0)
							(*x)[lid] *= factor;
					}
				}
			}
}

//==================================================================
void Ocean::applyLandMask(Teuchos::RCP<Epetra_Vector> x,
						  LandMask maskA, LandMask maskB)
{
	INFO("Ocean: applyLandMask...");
	int nmask = maskA.global->size();
	assert(nmask == (M_+2)*(N_+2)*(L_+2));
	
	int idx1, idx2;
	int ii, lid;
	int iil, iir;
	int dir;
	int nnz;
	double value, globValue;
	double avg = 0.0;
	int newOcean = 0;
	int newLand  = 0;
	std::vector<int> unknowns = {4, 5}; // Only adjust these unknowns (T,S)

	for (int q = 0; q != 1; ++q)
		for (int k = 1; k != L_+1; ++k)
			for (int j = 1; j != M_+1; ++j)
				for (int i = 1; i != N_+1; ++i)
				{
					idx1 = k*(M_+2)*(N_+2) + j*(N_+2) + i;
					idx2 = (k-1)*M_*N_ + (j-1)*N_ + i-1;
				
					dir = (*maskA.global)[idx1] - (*maskB.global)[idx1];
				
					if (dir > 0)  // New ocean
					{
						for (auto &var: unknowns)
						{
							// row in vector
							ii  = idx2*_NUN_ + var;
							
							// left-most point at this latitude in this layer
							iil = ((k-1)*M_*N_ + (j-1)*N_)*_NUN_ + var;

							// right-most point at this latitude in this layer
							iir = ((k-1)*M_*N_ + (j-1)*N_ + N_-1)*_NUN_ + var;

							// get zonal average
							avg = 0.0;
							nnz = 0;
							for (; iil <= iir; iil += _NUN_)
							{
								value     = 0.0;
								globValue = 0.0;

								lid = x->Map().LID(iil);

								if (lid >= 0)
									value = (*x)[lid];
								
								comm_->SumAll(&value, &globValue, 1);
								avg += globValue;
								nnz += (std::abs(globValue) < 1e-8) ? 0 : 1;
							}

							// put zonal average in ii
							if (nnz > 0)
							{
								avg /= nnz;
								lid = x->Map().LID(ii);
								if (lid >= 0 && std::abs((*x)[lid]) < 1e-8)
								{
									newOcean++;
									(*x)[lid] = avg;
								}
								
								// synchronize counter
								int tmp;
								comm_->MaxAll(&newOcean, &tmp, 1);
								newOcean = tmp;									
							}
						}
					}
					else if (dir < 0) // New land
					{
						for (ii = idx2*_NUN_; ii != (idx2+1)*_NUN_; ++ii)
						{
							newLand++;
							lid = x->Map().LID(ii); // Local ID of this point			
							if (lid >= 0)
								(*x)[lid] = 0.0;
						}
					}
				}

	// To be sure we check whether we have introduced some unwanted values
	inspectVector(x);
	
	INFO("Ocean: applyLandmask, adjusted " << newOcean << " new ocean entries.");
	INFO("Ocean: applyLandmask, adjusted " << newLand  << " new land entries.");
	INFO("Ocean: applyLandMask... done");
}

//====================================================================
void Ocean::preProcess()
{
	// Enable computation of preconditioner
	recompPreconditioner_ = true;
	INFO("Ocean pre-processing: enabling computation of preconditioner.");
	
	// Output datafiles (fort.3 fort.44)
	printFiles(); // not sure if this is not too much...
	
}

//====================================================================
void Ocean::postProcess()
{
	if (saveState_)
		saveStateToFile(outputFile_); // Save to hdf5
	
	printFiles(); // Print in standard fortran format 
	
	if (storeEverything_)
		copyFiles(); // Copy fortran and hdf5 files
}

//=====================================================================
// Setup block preconditioner parameters
// --> xml files should have a better home
void Ocean::initializePreconditioner()
{
 	INFO("Ocean: initialize preconditioner...");
	TIMER_START("Ocean: initialize preconditioner");

	Teuchos::RCP<Teuchos::ParameterList> precParams =
		Teuchos::rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("ocean_preconditioner_params.xml",
								precParams.ptr());	

	// Create and initialize block preconditioner
	precPtr_ = Teuchos::rcp(new TRIOS::BlockPreconditioner
							(jac_, domain_, *precParams));
	
	precPtr_->Initialize();  // Initialize

	precInitialized_ = true;

	// Enable computation of preconditioner
	recompPreconditioner_ = true;

	TIMER_STOP("Ocean: initialize preconditioner");
	INFO("Ocean: initialize preconditioner done...");
}	

//=====================================================================
void Ocean::initializeSolver()
{
	INFO("Ocean: initialize solver...");
	
	solverParams_ = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams_.ptr());
	
	// Get the requested solver type
	solverType_ = solverParams_->get("Ocean solver type", 'F');
	recompTol_  = solverParams_->get("Tolerance recompute preconditioner", 0.9);
	
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
	int maxrestarts = solverParams_->get("FGMRES restarts", 0);
	int output      = solverParams_->get("FGMRES output", 100);
	
	int NumGlobalElements = state_->GlobalLength();
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
	belosParamList_->set("Output Frequency", output);
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

	// Get new preconditioner
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

		if (tol > recompTol_) // stagnation, maybe a new precon helps
		{
			INFO("Ocean: FGMRES, stagnation: " << recompTol_);
			recompPreconditioner_ = true;
		}
					 
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

//==================================================================
void Ocean::printResidual(VectorPtr rhs)
{
	VectorPtr Ax = getSolution('C');
	VectorPtr  x = getSolution('C');
	applyMatrix(*x, *Ax);
	Ax->update(1.0, *rhs, -1.0);
	Ax->norm('E',"ocean ||b-Ax||");
	Ax->print("residual");
}

//=====================================================================
void Ocean::scaleProblem(VectorPtr rhs)
{
	INFO("Ocean: scale problem...");
	if (rhs == Teuchos::null)
		ERROR("DEPRECATED FUNCTIONALITY", __LINE__, __FILE__);
	
	// Not sure if this is the right approach and/or implemented correctly.
	// Scaling is obtained from THCM and then applied to the problem.

	rowScaling_ = THCM::Instance().getRowScaling();
	// colScaling_ = THCM::Instance().getColScaling();
	
	//------------------------------------------------------
	if (rowScalingRecipr_ == Teuchos::null or
		!rowScaling_->Map().SameAs(rowScalingRecipr_->Map()))
	{
		rowScalingRecipr_ =
			rcp(new Epetra_Vector(rowScaling_->Map()));
	}

	// //------------------------------------------------------
	// if (colScalingRecipr_ == Teuchos::null or
	// 	!colScaling_->Map().SameAs(colScalingRecipr_->Map()))
	// {
	// 	colScalingRecipr_ =
	// 		rcp(new Epetra_Vector(colScaling_->Map()));
	// }
	
	rowScaling_ = THCM::Instance().getRowScaling();
	// jac_->InvRowSums(*rowScaling_);
	*rowScalingRecipr_ = *rowScaling_;
	rowScalingRecipr_->Reciprocal(*rowScaling_);

	// *colScalingRecipr_ = *colScaling_;
	// colScalingRecipr_->Reciprocal(*colScaling_);

	//------------------------------------------------------
	// double nrm;
	// sol_->Norm2(&nrm);
	// DEBUG("Ocean::scaleProblem() sol (before scaling): "  << nrm);
	
	//------------------------------------------------------
	// jac_->LeftScale(*rowScalingRecipr_);
	jac_->LeftScale(*rowScaling_);
	// jac_->RightScale(*colScaling_);

 	// (rhs->getOceanVector())->Multiply(1.0, *rowScalingRecipr_,
	// 								  *(rhs->getOceanVector()), 0.0);
	
	(rhs->getOceanVector())->Multiply(1.0, *rowScaling_,
									  *(rhs->getOceanVector()), 0.0);

	
	// recompPreconditioner_ = true;
	// sol_->ReciprocalMultiply(1.0, *colScaling, *sol_, 0.0);

	//------------------------------------------------------
	//sol_->Norm2(&nrm);
	//DEBUG("Ocean::scaleProblem() sol (after scaling): "  << nrm);
	INFO("Ocean: scale problem... done");
}

//=====================================================================
void Ocean::unscaleProblem(VectorPtr rhs)
{
	INFO("Ocean: unscale problem...");
	if (rhs == Teuchos::null)
		ERROR("DEPRECATED FUNCTIONALITY", __LINE__, __FILE__);
		
	//------------------------------------------------------
	//double nrm;
	//sol_->Norm2(&nrm);
	//DEBUG("Ocean::unscaleProblem() sol (before unscaling): " << nrm);
	
	//------------------------------------------------------
	// jac_->RightScale(*colScalingRecipr_);
	jac_->LeftScale(*rowScalingRecipr_);
	// jac_->LeftScale(*rowScaling_);
	
	(rhs->getOceanVector())->Multiply(1.0, *rowScalingRecipr_,
									  *(rhs->getOceanVector()), 0.0);
	// (rhs->getOceanVector())->Multiply(1.0, *rowScaling_,
	// 								  *(rhs->getOceanVector()), 0.0);

	
	//sol_->Multiply(1.0, *rowScaling_, *sol_, 0.0);

	//------------------------------------------------------
	//sol_->Norm2(&nrm);
	//DEBUG("Ocean::unscaleProblem() sol (after unscaling): " << nrm);
	INFO("Ocean: unscale problem... done");
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getRowScaling()
{
	return THCM::Instance().getRowScaling();
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getColScaling()
{
	return THCM::Instance().getColScaling();
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
Teuchos::RCP<Epetra_Vector> Ocean::getVector(char mode, RCP<Epetra_Vector> vec)
{
	if (mode == 'C') // copy
	{
		RCP<Epetra_Vector> copy = rcp(new Epetra_Vector(*vec));
		return copy;
	}
	else if (mode == 'V') // view
	{
		return vec;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getSolution(char mode)
{
	return getVector(mode, sol_);
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getState(char mode)
{
	return getVector(mode, state_);
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getRHS(char mode)
{
	return getVector(mode, rhs_);
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getM(char mode)
{
	RCP<Epetra_Vector> vecM =
		rcp(new Epetra_Vector(*sol_));
	
	vecM->PutScalar(1.0);
	
	int NumMyElements     = vecM->Map().NumMyElements();
	int *MyGlobalElements = vecM->Map().MyGlobalElements();
	
	for (int i = 0; i != NumMyElements; ++i)
	{
		// Ignoring w and p equations
		if ( ((MyGlobalElements[i] % _NUN_) == (WW-1)) ||
			 ((MyGlobalElements[i] % _NUN_) == (PP-1)) )
			(*vecM)[i] = 0;
	}

	return getVector(mode, vecM);
}

//====================================================================
void Ocean::applyMatrix(Epetra_Vector const &v, Epetra_Vector &out)
{
	TIMER_START("Ocean: apply matrix...");
	jac_->Apply(v, out);
	TIMER_STOP("Ocean: apply matrix...");
}

//====================================================================
void Ocean::buildPreconditioner(bool forceInit)
{
	if (!precInitialized_ || forceInit) // Initialize preconditioner
		initializePreconditioner();

	if (recompPreconditioner_)
	{
		TIMER_START("Ocean: compute preconditioner");
		INFO("Ocean: compute preconditioner...");
		precPtr_->Compute();
		INFO("Ocean: compute preconditioner... done");
		TIMER_STOP("Ocean: compute preconditioner");
		recompPreconditioner_ = false;  // Disable subsequent recomputes
	}
}

//====================================================================
void Ocean::applyPrecon(Epetra_Vector const &v, Epetra_Vector &out)
{
	if (!precInitialized_) // Initialize preconditioner
		initializePreconditioner();
	
	// Compute preconditioner
	buildPreconditioner();	

	TIMER_START("Ocean: apply preconditioning...");
	precPtr_->ApplyInverse(v, out);
	TIMER_STOP("Ocean: apply preconditioning...");
}

//====================================================================
void Ocean::setAtmosphere(Teuchos::RCP<Epetra_Vector> const &atmos)
{
	TIMER_START("Ocean: set atmosphere...");
	// This is a job for THCM
	THCM::Instance().setAtmosphere(atmos);
	TIMER_STOP("Ocean: set atmosphere...");
}

//====================================================================
// Return the coupling block containing the contribution of the
// atmosphere to the ocean
// --> Parallelize
void Ocean::getAtmosBlock(std::vector<double> &values,
						  std::vector<int> &row_inds)
{
	double Ooa, Os;
	FNAME(getooa)(&Ooa, &Os);
	values = std::vector<double>(N_ * M_, -Ooa);
	INFO("Ocean: Ooa = " << Ooa);
	
	// Apply Surface mask to values
	int ctr  = 0;
	int lctr = 0;
	for (int j = 0; j != M_; ++j)
		for (int i = 0; i != N_; ++i)
		{			
			if ((*surfmask_)[j*N_+i] == 1)
			{
				values[ctr] = 0.0;
				lctr++;
			}
			ctr++;
		}

	// clear row_inds
	row_inds = std::vector<int>(N_ * M_, 0);
	
	// fill with surface temperature values
	int idx = 0;
	for (int j = 0; j != M_; ++j)
		for (int i = 0; i != N_; ++i)
		{
			row_inds[idx] = (FIND_ROW2(_NUN_, N_, M_, L_,i,j,L_-1,TT));
			idx++;
		}
	INFO("  A->O block, zeros due to surfacemask --> " << lctr);
}

//====================================================================
// Fill and return a copy of the surface temperature
Teuchos::RCP<Epetra_Vetor> Ocean::getSurfaceT()
{
	TIMER_START("Ocean: get surface temperature...");
	
	sst_->Import(*state_, *surfaceTimporter_, Insert);
	
	TIMER_STOP("Ocean: get surface temperature...");	
	return sst_;
}

//=====================================================================
void Ocean::printFiles()
{	
	Teuchos::RCP<Epetra_MultiVector> solution = Utils::Gather(*state_, 0);
	Teuchos::RCP<Epetra_MultiVector> rhs      = Utils::Gather(*rhs_, 0);
	int filename = 3;
	int label    = 2;
	int length   = solution->GlobalLength();
	double *solutionArray = new double[length];
	double *rhsArray = new double[length]; 
	if (comm_->MyPID() == 0)
	{
		INFO( "Writing to fort." << filename 
			  << " at label " << label);

		(*solution)(0)->ExtractCopy(solutionArray); 
		(*rhs)(0)->ExtractCopy(rhsArray);
		
		FNAME(write_data)(solutionArray, &filename, &label);
	}
	
	delete [] solutionArray;
	delete [] rhsArray;
}

//=====================================================================
void Ocean::copyFiles()
{
	if (comm_->MyPID() == 0)
	{
		// Copy fort.3
		std::stringstream ss;
		ss << "ocean_state_par" << std::setprecision(4) << std::setfill('_')
		   << std::setw(2) << THCM::Instance().par2int(parName_) << "_"
		   << std::setw(6) << getPar(parName_);

		INFO("copying fort.3 to " << ss.str());
		std::ifstream src1("fort.3", std::ios::binary);
		std::ofstream dst1(ss.str(), std::ios::binary);
		dst1 << src1.rdbuf();

		if (saveState_) // Copy hdf5
		{
			ss << ".h5";
			INFO("copying " << outputFile_ << " to " << ss.str());
			std::ifstream src2(outputFile_.c_str(), std::ios::binary);
			std::ofstream dst2(ss.str(), std::ios::binary);
			dst2 << src2.rdbuf();
		}
	}
}

//=====================================================================
void Ocean::copyFiles(std::string const &filename)
{
	if (comm_->MyPID() == 0)
	{		
		INFO("copying fort.3 to " << filename);
		std::ifstream src1("fort.3", std::ios::binary);
		std::ofstream dst1(filename.c_str(), std::ios::binary);
		dst1 << src1.rdbuf();
		
		if (saveState_) // Copy hdf5
		{
			std::string fnameCpy(filename);
			fnameCpy.append(".h5");
			INFO("copying " << outputFile_ << " to " << fnameCpy);
			std::ifstream src2(outputFile_.c_str(), std::ios::binary);
			std::ofstream dst2(fnameCpy.c_str(), std::ios::binary);
			dst2 << src2.rdbuf();
		}
	}
}

//==================================================================
void Ocean::copyMask(std::string const &filename)
{
	if (comm_->MyPID() == 0)
	{		
		if (saveMask_) // Copy fort.44
		{
			std::string fnameCpy(filename);
			fnameCpy.append(".mask");
			INFO("copying " << "fort.44" << " to " << fnameCpy);
			std::ifstream src2("fort.44", std::ios::binary);
			std::ofstream dst2(fnameCpy.c_str(), std::ios::binary);
			dst2 << src2.rdbuf();
		}
		else
			WARNING("Saving mask not enabled in Ocean.", __FILE__, __LINE__);
	}
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

//====================================================================
void Ocean::setParameters(Teuchos::RCP<Teuchos::ParameterList> pars)
{
	std::string parName;
	double parValue;
	// This is similar to reading from HDF5
	for (int par = 1; par <= _NPAR_; ++par)
	{
		parName  = THCM::Instance().int2par(par);
		parValue = getPar(parName);
		
		// Overwrite continuation parameter and put it in THCM
		try
		{
			parValue = pars->get(parName, parValue);
		}
		catch (EpetraExt::Exception &e)
		{
			e.Print();
			continue;
		}
		
		setPar(parName, parValue);
	}	
}

