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

//=====================================================================
#include <math.h>

//=====================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//=====================================================================
// Get access to a few THCM functions
extern "C" _SUBROUTINE_(write_data)(double*, int*, int*);
extern "C" _SUBROUTINE_(getparcs)(int*, double*);
extern "C" _SUBROUTINE_(setparcs)(int*,double*);
extern "C" _SUBROUTINE_(getdeps)(double*, double*, double*, double*, double *);
extern "C" _SUBROUTINE_(get_constants)(double*, double*, double*);
extern "C" _SUBROUTINE_(set_ep_constants)(double*, double*, double*, double*, double*);

//=====================================================================
// Constructor:
Ocean::Ocean(RCP<Epetra_Comm> Comm, RCP<Teuchos::ParameterList> oceanParamList)
    :
    comm_                (Comm),   // Setting the communication object
    solverInitialized_   (false),  // Solver needs initialization
    precInitialized_     (false),  // Preconditioner needs initialization
    recompPreconditioner_(true),   // We need a preconditioner to start with
    recompMassMat_       (true),   // We need a mass matrix to start with

    inputFile_           (oceanParamList->get("Input file",  "ocean_input.h5")),
    outputFile_          (oceanParamList->get("Output file", "ocean_output.h5")),
    
    loadState_           (oceanParamList->get("Load state", false)),
    saveState_           (oceanParamList->get("Save state", true)),
    saveMask_            (oceanParamList->get("Save mask", true)),
    loadSalinityFlux_    (oceanParamList->get("Load salinity flux", false)),
    saveSalinityFlux_    (oceanParamList->get("Save salinity flux", true)),
    loadTemperatureFlux_ (oceanParamList->get("Load temperature flux", false)),
    saveTemperatureFlux_ (oceanParamList->get("Save temperature flux", true)),
    
    saveEveryStep_       (oceanParamList->get("Save every step", false)),
    useFort3_            (oceanParamList->get("Use legacy fortran output", false)),

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

    // Throw a few errors if the parameters are odd
    if ((thcm_->getSRES() || thcm_->getITS()) && loadSalinityFlux_)
    {
        ERROR(" (SRES = 1 || ite = 1 ) => loadSalinityFlux_ = false",
              __FILE__, __LINE__);
    }
    
    if ((thcm_->getTRES() || thcm_->getITE()) && loadTemperatureFlux_)
    {
        ERROR(" (TRES = 1 || its = 1) => loadTemperatureFlux_ = false",
              __FILE__, __LINE__);
    }

    // Obtain solution vector from THCM
    state_ = THCM::Instance().getSolution();
    INFO("Ocean: Solution obtained from THCM");

    // Get domain object and get the problem dimensions
    domain_ = THCM::Instance().GetDomain();

    N_ = domain_->GlobalN();
    M_ = domain_->GlobalM();
    L_ = domain_->GlobalL();

    // grid representation of te state
    grid_   = rcp(new OceanGrid(domain_));

    // Read starting parameters from xml
    Teuchos::ParameterList& startList =
        thcmList.sublist("Starting Parameters");
    THCM::Instance().ReadParameters(startList);

    // If specified we load a pre-existing state and parameters (x,l)
    // thereby overwriting the starting parameters
    // This will be able to load salinity and temperature fluxes as well.
    if (loadState_ || loadSalinityFlux_)
        loadStateFromFile(inputFile_);
    else            // Initialize with trivial solution
        state_->PutScalar(0.0);

    // Now that we have the state and parameters initialize
    // the Jacobian, solution and rhs
    initializeOcean();

    // Obtain adjusted landmask
    landmask_ = getLandMask("current");

    // Initialize preconditioner
    initializePreconditioner();

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
    tIndexMap_ =
        Utils::CreateSubMap(state_->Map(), tRows);

    // Create the SST vector
    sst_ = Teuchos::rcp(new Epetra_Vector(*tIndexMap_));

    // Create importer
    // Target map: tIndexMap
    // Source map: state_->Map()
    surfaceTimporter_ =
        Teuchos::rcp(new Epetra_Import(*tIndexMap_, state_->Map()));

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

    // Obtain mass matrix B from THCM. Note that we assume the mass
    // matrix is independent of the state and parameters so we only
    // compute it when asked for through recompMassMat_ flag.
    buildMassMat();

    // Initialize solution and rhs
    sol_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
    rhs_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));

    // Compute right hand side and print its norm
    computeRHS();
    
#ifdef DEBUGGING_NEW
    INFO("Ocean: initialization: ||F|| = " << Utils::norm(rhs_));
    Utils::save(state_, "initialstate");
    Utils::save(rhs_, "initialrhs");
#endif    

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
// --> This can be more efficient
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

    // Copy to tmp
    std::vector<int> tmp(*mask.global);

    // Erase everything but upper 2 layers
    tmp.erase(tmp.begin(), tmp.begin() + tmp.size() -
              ( 2 * (M_+2) * (N_+2) ));

    // Create global surface mask rcp
    mask.global_surface = std::make_shared<std::vector<int> >();

    // Put the first layer of the remaining tmp mask
    // in mask.global_surface, without borders.
    for (int j = 1; j != M_+1; ++j)
        for (int i = 1; i != N_+1; ++i)
        {
            mask.global_surface->push_back(tmp[j*(N_+2) + i]);
        }

    assert( (int) mask.global_surface->size() == N_*M_ );


    // Set label
    mask.label = fname;

    // Return the struct
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
    recompMassMat_        = true;
    INFO("Ocean pre-processing:");
    INFO("                      enabling computation of preconditioner.");
    INFO("                      enabling computation of mass matrix.");

    // Output datafiles (hdf5, fort.44)
    printFiles(); // not sure if this is not too much...

}

//====================================================================
void Ocean::postProcess()
{
    if (saveState_)
        saveStateToFile(outputFile_); // Save to hdf5

    printFiles(); // Print in standard fortran format

    if (saveEveryStep_)
        copyFiles(); // Copy fortran and hdf5 files
}

//=====================================================================
std::string const Ocean::writeData(bool describe)
{
    std::ostringstream datastring;
    if (describe)
    {
        if (solverInitialized_)
        {
            datastring << std::setw(_FIELDWIDTH_/2)
                       << "MV";
            datastring << std::setw(_FIELDWIDTH_)
                       << "Tol";

        }
            
        datastring << std::setw(_FIELDWIDTH_)
                   << "max(Psi)"
                   << std::setw(_FIELDWIDTH_) 
                   << "min(Psi)";
        
        return datastring.str();
    }
    else
    {
        datastring.precision(_PRECISION_);
    
        // compute streamfunctions and output data
        grid_->ImportData(*state_);
        double psiMax = grid_->psimMax();
        double psiMin = grid_->psimMin();

        double r0dim, udim, hdim;
        FNAME(get_constants)(&r0dim, &udim, &hdim);

        const double transc = r0dim * hdim * udim;
        psiMax = psiMax * transc * 1e-6; // conversion to Sv
        psiMin = psiMin * transc * 1e-6; //

        if (solverInitialized_)
        {
            datastring << std::scientific << std::setw(_FIELDWIDTH_/2)
                       << belosSolver_->getNumIters();
            datastring << std::scientific << std::setw(_FIELDWIDTH_)
                       << belosSolver_->achievedTol();
        }

        datastring << std::scientific << std::setw(_FIELDWIDTH_)
                   << psiMax
                   << std::setw(_FIELDWIDTH_)
                   << psiMin;
        
        return datastring.str();    
    }                           
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
    recompTol_  = solverParams_->get("Tolerance recompute preconditioner", 0.999);

    // Initialize the preconditioner
    if (!precInitialized_)
        initializePreconditioner();

    // Initialize the requested solver
    if (solverType_ == 'F')
        initializeBelos();
    else
        ERROR("No solver specified", __FILE__, __LINE__);

    solverInitialized_ = true;

    // Now that the solver and preconditioner are initialized we are allowed to
    // perform a solve.
    INFO("Ocean: initialize solver... done");
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
void Ocean::solve(Teuchos::RCP<Epetra_MultiVector> rhs)
{
    // Check whether solver is initialized, if not perform the
    // initialization here
    if (!solverInitialized_)
        initializeSolver();

    // Get new preconditioner
    buildPreconditioner();

    // Use trivial initial solution
    sol_->PutScalar(0.0);

    // Set the problem, rhs may be given as an argument to solve().
    if (solverType_ == 'F')
    {
        bool set;
        if (rhs == Teuchos::null)
            set = problem_->setProblem(sol_, rhs_);
        else
            set = problem_->setProblem(sol_, rhs);
        TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
                                   "*** Belos::LinearProblem failed to setup");
    }
    else
        ERROR("No solver specified", __FILE__, __LINE__);

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
            belosSolver_->solve();      // Solve

        }
        catch (std::exception const &e)
        {
            ERROR("Ocean: exception caught: " << e.what(), __FILE__, __LINE__);
        }
    }
    else
    {
        ERROR("No solve specified", __FILE__, __LINE__);
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
    else
    {
        ERROR("No solve specified", __FILE__, __LINE__);
    }
}

//=====================================================================
double Ocean::explicitResNorm(VectorPtr rhs)
{
    RCP<Epetra_Vector> Ax =
        rcp(new Epetra_Vector(*(domain_->GetSolveMap())));
    jac_->Apply(*sol_, *Ax);
    Ax->Update(1.0, *rhs, -1.0);
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
    Ax->Update(1.0, *rhs, -1.0);
    INFO("Ocean: ||b - A*x|| = " << Utils::norm(Ax));
    Utils::print(Ax, "residual");
}

//=====================================================================
void Ocean::scaleProblem(VectorPtr rhs)
{
    INFO("Ocean: scale problem...");
    if (rhs == Teuchos::null)
        ERROR("DEPRECATED FUNCTIONALITY", __FILE__, __LINE__);

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
    //      !colScaling_->Map().SameAs(colScalingRecipr_->Map()))
    // {
    //      colScalingRecipr_ =
    //          rcp(new Epetra_Vector(colScaling_->Map()));
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
    //                                    *(rhs->getOceanVector()), 0.0);

    rhs->Multiply(1.0, *rowScaling_, *rhs, 0.0);

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
        ERROR("DEPRECATED FUNCTIONALITY", __FILE__, __LINE__);

    //------------------------------------------------------
    //double nrm;
    //sol_->Norm2(&nrm);
    //DEBUG("Ocean::unscaleProblem() sol (before unscaling): " << nrm);

    //------------------------------------------------------
    // jac_->RightScale(*colScalingRecipr_);
    jac_->LeftScale(*rowScalingRecipr_);
    // jac_->LeftScale(*rowScaling_);

    rhs->Multiply(1.0, *rowScalingRecipr_, *rhs, 0.0);
    // (rhs->getOceanVector())->Multiply(1.0, *rowScaling_,
    //                                    *(rhs->getOceanVector()), 0.0);


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
    THCM::Instance().fixMixing(0);
    THCM::Instance().evaluate(*state_, rhs_, false);
    TIMER_STOP("Ocean: compute RHS...");
}

//=====================================================================
void Ocean::computeJacobian()
{
    TIMER_START("Ocean: compute Jacobian...");

    // Compute the Jacobian in THCM using the current state
    THCM::Instance().fixMixing(0);
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
Teuchos::RCP<Epetra_Vector> Ocean::getDiagB(char mode)
{
    return getVector(mode, diagB_);
}

//====================================================================
// Obtain M, binary diagonal matrix to select transient parts. -->
// this does not include a zero for an integral condition in the case
// of non-restoring salinity forcing.
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
void Ocean::applyMatrix(Epetra_MultiVector const &v, Epetra_MultiVector &out)
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
void Ocean::applyPrecon(Epetra_MultiVector const &v, Epetra_MultiVector &out)
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
// Not yet implemented!
void Ocean::applyDiagInv(Epetra_MultiVector const &v, Epetra_MultiVector &out)
{
    assert(v.Map().SameAs(out.Map()));

    Epetra_Vector diag(*rhs_);
    diag.PutScalar(0.0);
    jac_->ExtractDiagonalCopy(diag);

    int numMyElements = diag.Map().NumMyElements();
    for (int i = 0; i != numMyElements; ++i)
        if (diag[i] == 0)
            std::cout << "pid: " << diag.Map().Comm().MyPID() << " "
                      << i << ": zero element in diagonal" << std::endl;

    DUMP_VECTOR("diag", diag);
}

//====================================================================
void Ocean::buildMassMat()
{
    if (recompMassMat_)
    {
        THCM::Instance().evaluateB();
        diagB_ = THCM::Instance().DiagB();
    }
    recompMassMat_ = false; // Disable subsequent recomputes
}

//====================================================================
void Ocean::applyMassMat(Epetra_MultiVector const &v, Epetra_MultiVector &out)
{
    // Compute mass matrix
    buildMassMat();

    // element-wise multiplication (out = 0.0*out + 1.0*B*v)
    out.Multiply(1.0, *diagB_, v, 0.0);
}

//====================================================================
void Ocean::synchronize(std::shared_ptr<AtmospherePar> atmos)
{
    TIMER_START("Ocean: set atmosphere...");

    // Obtain and set atmosphere T at the interface
    Teuchos::RCP<Epetra_Vector> atmosT  = atmos->interfaceT();
    THCM::Instance().setAtmosphereT(atmosT);

    // Obtain and set humidity field at the interface
    Teuchos::RCP<Epetra_Vector> atmosQ  = atmos->interfaceQ();
    THCM::Instance().setAtmosphereQ(atmosQ);

    // Obtain and set precipitation field at the interface
    Teuchos::RCP<Epetra_Vector> atmosP  = atmos->interfaceP();
    THCM::Instance().setAtmosphereP(atmosP);

    // We also need to know a few atmospheric constants to compute E,
    // P and their derivatives w.r.t. SST (To) and humidity (q) These
    // could be obtained at construction but, as they may depend on
    // continuation parameters, the call belongs here.
    double qdim, nuq, eta, dqso, dqdt, Eo0;
    atmos->getConstants( qdim, nuq, eta, dqso, dqdt, Eo0);

    FNAME(set_ep_constants)( &qdim, &nuq, &eta, &dqso, &Eo0 );

    TIMER_STOP("Ocean: set atmosphere...");
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getLocalAtmosT()
{
    return THCM::Instance().getLocalAtmosT();
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getLocalAtmosQ()
{
    return THCM::Instance().getLocalAtmosQ();
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getLocalAtmosP()
{
    return THCM::Instance().getLocalAtmosP();
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getLocalOceanE()
{
    return THCM::Instance().getLocalOceanE();
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getE()
{
    return THCM::Instance().getOceanE();
}

//==================================================================
std::shared_ptr<Utils::CRSMat> Ocean::getBlock(std::shared_ptr<AtmospherePar> atmos)
{
    // initialize empty CRS matrix
    std::shared_ptr<Utils::CRSMat> block = std::make_shared<Utils::CRSMat>();


    // this block has values -Ooa on the surface temperature points
    double Ooa, Os, gamma, eta, lvscq;
    FNAME(getdeps)(&Ooa, &Os, &gamma, &eta, &lvscq);

    INFO("Ocean: getBlock() atmos specialization: Ooa = " << Ooa
         << ", gamma = " << gamma << ", eta = " << eta << ", lvscq = " << lvscq);

    int T = 1; // (1-based) in the Atmosphere, temperature is the first unknown
    int Q = 2; // (1-based) in the Atmosphere, humidity is the second unknown
    int P = 3; // (1-based) in the Atmosphere, precipitation is the final unknown

    double rowIntCon = THCM::Instance().getRowIntCon();

    // fill CRS struct
    int el_ctr = 0;
    int col;
    for (int k = 0; k != L_; ++k)
        for (int j = 0; j != M_; ++j)
            for (int i = 0; i != N_; ++i)
                for (int xx = UU; xx <= SS; ++xx)
                {
                    block->beg.push_back(el_ctr);
                    if ( (k == L_-1) && (xx == TT) ) // surface T row
                    {
                        if ((*landmask_.global_surface)[j*N_+i] == 0) // non-land
                        {
                            // tatm dependency
                            block->co.push_back(-Ooa);
                            block->jco.push_back(atmos->interface_row(i,j,T) );
                            el_ctr++;

                            // // qatm dependency
                            // block->co.push_back(-lvscq*eta);
                            // block->jco.push_back(atmos->interface_row(i,j,Q) );
                            // el_ctr++;
                        }
                    }
                    else if ( (k == L_-1) && (xx == SS) && // surface S row
                              FIND_ROW2(_NUN_, N_, M_, L_, i, j, k, xx) != rowIntCon)
                    {
                        if ((*landmask_.global_surface)[j*N_+i] == 0) // non-land
                        {
                            // humidity dependency
                            block->co.push_back(gamma*eta);
                            block->jco.push_back(atmos->interface_row(i,j,Q) );
                            el_ctr++;

                            // precipitation dependency
                            col = atmos->interface_row(i,j,P);
                            if (col >= 0)
                            {
                                block->co.push_back(gamma);
                                block->jco.push_back(col);
                                el_ctr++;
                            }
                        }
                    }
                }

    // final entry in beg ( == nnz)
    block->beg.push_back(el_ctr);

    assert( (int) block->co.size() == block->beg.back());

    return block;
}

//====================================================================
// Fill and return a copy of the surface temperature
Teuchos::RCP<Epetra_Vector> Ocean::interfaceT()
{
    TIMER_START("Ocean: get surface temperature...");

    if (!(sst_->Map().SameAs(*tIndexMap_)))
    {
        CHECK_ZERO(sst_->ReplaceMap(*tIndexMap_));
    }
    CHECK_ZERO(sst_->Import(*state_, *surfaceTimporter_, Insert));
    TIMER_STOP("Ocean: get surface temperature...");
    return sst_;
}

//=====================================================================
void Ocean::printFiles()
{
    // dummy set up, exit write function after creating fort.44
    int filename = 0;
    int label    = 0;
    int length   = 1;
    
    Teuchos::RCP<Epetra_MultiVector> solution;
    Teuchos::RCP<Epetra_MultiVector> rhs;
    
    if (useFort3_)
    {
        solution = Utils::Gather(*state_, 0);
        rhs      = Utils::Gather(*rhs_, 0);

        filename = 3;
        label    = 2;
        length   = solution->GlobalLength();
    }

    double *solutionArray = new double[length];
    double *rhsArray = new double[length];

    if (comm_->MyPID() == 0)
    {
        INFO( "Writing to fort." << filename
              << " at label " << label);

        if (useFort3_)
        {
            (*solution)(0)->ExtractCopy(solutionArray);
            (*rhs)(0)->ExtractCopy(rhsArray);
        }
        
        FNAME(write_data)(solutionArray, &filename, &label);
    }

    delete [] solutionArray;
    delete [] rhsArray;
}

//=====================================================================
void Ocean::copyFiles(std::string const &filename)
{
    if ( (comm_->MyPID() == 0) && saveState_)
    {
        std::stringstream ss;
        if (filename.length() == 0)
        {
            ss << "ocean_state_par" << std::setprecision(4) << std::setfill('_')
               << std::setw(2) << THCM::Instance().par2int(parName_) << "_"
               << std::setw(6) << getPar(parName_);

        }
        else
        {
            ss << filename;
        }
        
        ss << ".h5";
        INFO("copying " << outputFile_ << " to " << ss.str());
        std::ifstream src2(outputFile_.c_str(), std::ios::binary);
        std::ofstream dst2(ss.str(), std::ios::binary);
        dst2 << src2.rdbuf();
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
    INFO("_________________________________________________________");
    INFO("Writing ocean state and parameters to " << filename);

    INFO("   ocean state: ||x|| = " << Utils::norm(state_));

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
    
    if (saveSalinityFlux_)
    {
        // Write emip to ocean output file
        INFO("Writing salinity flux to " << filename);
        Teuchos::RCP<Epetra_Vector> salflux = THCM::Instance().getSalinityFlux();
        HDF5.Write("SalinityFlux", *salflux);        
    }

    if (saveTemperatureFlux_)
    {
        // Write emip to ocean output file
        INFO("Writing temperature flux to " << filename);
        Teuchos::RCP<Epetra_Vector> temflux = THCM::Instance().getTemperatureFlux();
        HDF5.Write("TemperatureFlux", *temflux);        
    }
    
    INFO("_________________________________________________________");
    
    return 0;
}

// =====================================================================
int Ocean::loadStateFromFile(std::string const &filename)
{
    // Create HDF5 object
    EpetraExt::HDF5 HDF5(*comm_);
    Epetra_MultiVector *readState;

    // Read state
    HDF5.Open(filename);

    if (loadState_)
    {
        INFO("Loading state from " << filename);

        // To be sure that the state is properly initialized,
        // we obtain the state vector from THCM.
        state_ = THCM::Instance().getSolution();

        // Check whether file exists
        std::ifstream file(filename);
        if (!file)
        {
            WARNING("Can't open " << filename
                    << ", continue with trivial state", __FILE__, __LINE__);

            // initialize trivial ocean
            state_->PutScalar(0.0);
            return 1;
        }
        else
            file.close();

        if (!HDF5.IsContained("State"))
        {
            ERROR("The group <State> is not contained in hdf5 " << filename,
                  __FILE__, __LINE__);
        }
        
        HDF5.Read("State", readState);      

        // Create import strategies
        // target map: thcm domain SolveMap
        // source map: state with linear map as read by HDF5.Read
        Teuchos::RCP<Epetra_Import> lin2solve =
            Teuchos::rcp(new Epetra_Import(*(domain_->GetSolveMap()),
                                           readState->Map() ));

        // Import state from HDF5 into state_ datamember
        state_->Import(*((*readState)(0)), *lin2solve, Insert);

        INFO("   ocean state: ||x|| = " << Utils::norm(state_));

        if (THCM::Instance().getSRES() == 0)
            THCM::Instance().setIntCondCorrection(state_);
    
        INFO("Loading state from " << filename << " done");

        INFO("Loading parameters from " << filename);
        // Interface between HDF5 and the THCM parameters,
        // put all the (_NPAR_ = 30) THCM parameters back in THCM.
        std::string parName;
        double parValue;

        if (!HDF5.IsContained("Parameters"))
        {
            ERROR("The group <Parameters> is not contained in hdf5 " << filename,
                  __FILE__, __LINE__);
        }

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
          
        INFO("Loading parameters from " << filename << " done");
    }
    if (loadSalinityFlux_)
    {
        INFO("Loading salinity flux from " << filename);

        Epetra_MultiVector *readSalFlux;

        if (!HDF5.IsContained("SalinityFlux"))
        {
            ERROR("The group <SalinityFlux> is not contained in hdf5 " << filename,
                  __FILE__, __LINE__);
        }

        HDF5.Read("SalinityFlux", readSalFlux);

        // This should not be factorized as we cannot be sure what Map
        // is going to come out of the HDF5.Read call.

        Teuchos::RCP<Epetra_Vector> salflux = THCM::Instance().getSalinityFlux();
        
        Teuchos::RCP<Epetra_Import> lin2solve_surf =
            Teuchos::rcp(new Epetra_Import( salflux->Map(),
                                            readSalFlux->Map() ));

        salflux->Import(*((*readSalFlux)(0)), *lin2solve_surf, Insert);

        // Instruct THCM to set/insert this as the emip in the local model
        THCM::Instance().setEmip(salflux);
        

        if (HDF5.IsContained("AdaptedSalinityFlux"))
        {
            INFO(" detected AdaptedSalinityFlux in " << filename);
            Epetra_MultiVector *readAdaptedSalFlux;
            HDF5.Read("AdaptedSalinityFlux", readAdaptedSalFlux);
            
            assert(readAdaptedSalFlux->Map().SameAs(readSalFlux->Map()));
            
            Teuchos::RCP<Epetra_Vector> adaptedSalFlux =
                Teuchos::rcp(new Epetra_Vector( salflux->Map() ) );
            
            adaptedSalFlux->Import( *((*readAdaptedSalFlux)(0)), *lin2solve_surf, Insert);
            
            // Let THCM insert the adapted salinity flux
            THCM::Instance().setEmip(adaptedSalFlux, 'A');
        }

        if (HDF5.IsContained("AdaptedSalinityFlux_Mask"))
        {
            INFO(" detected AdaptedSalinityFlux_Mask in " << filename);
            Epetra_MultiVector *readSalFluxPert;
            HDF5.Read("AdaptedSalinityFlux_Mask", readSalFluxPert);
            
            assert(readSalFluxPert->Map().SameAs(readSalFlux->Map()));
            
            Teuchos::RCP<Epetra_Vector> salFluxPert =
                Teuchos::rcp(new Epetra_Vector( salflux->Map() ) );
            
            salFluxPert->Import( *((*readSalFluxPert)(0)), *lin2solve_surf, Insert);
            
            // Let THCM insert the salinity flux perturbation mask
            THCM::Instance().setEmip(salFluxPert, 'P');            
        }

        INFO("Loading salinity flux from " << filename << " done");
    }

    if (loadTemperatureFlux_)
    {
        INFO("Loading temperature flux from " << filename);

        Epetra_MultiVector *readTemFlux;
        HDF5.Read("TemperatureFlux", readTemFlux);

        // This should not be factorized as we cannot be sure what Map
        // is going to come out of the HDF5.Read call.

        Teuchos::RCP<Epetra_Vector> temflux = THCM::Instance().getTemperatureFlux();
        
        Teuchos::RCP<Epetra_Import> lin2solve_surf =
            Teuchos::rcp(new Epetra_Import( temflux->Map(),
                                            readTemFlux->Map() ));

        temflux->Import(*((*readTemFlux)(0)), *lin2solve_surf, Insert);

        // Instruct THCM to set/insert this as tatm in the local model
        THCM::Instance().setTatm(temflux);


        INFO("Loading temperature flux from " << filename << " done");
    }

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
