//=====================================================================
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Time.h>
#include <Epetra_Import.h>

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
#include "OceanGrid.H"
#include "THCM.H"
#include "Atmosphere.H"
#include "SeaIce.H"
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
extern "C" _SUBROUTINE_(getdeps)(double*, double*, double*,
                                 double*, double*, double*,
                                 double*);
extern "C" _SUBROUTINE_(get_parameters)(double*, double*, double*);
extern "C" _SUBROUTINE_(set_atmos_parameters)(Atmosphere::CommPars*);
extern "C" _SUBROUTINE_(set_seaice_parameters)(SeaIce::CommPars*);

// extern "C" double FNAME(qtoafun)(double*, double*, double*);

//=====================================================================
// Constructor:
Ocean::Ocean(RCP<Epetra_Comm> Comm)
    : Ocean(Comm, Teuchos::rcp(new Teuchos::ParameterList()))
{}

Ocean::Ocean(RCP<Epetra_Comm> Comm, Teuchos::RCP<Teuchos::ParameterList> oceanParamList)
    : Ocean(Comm, *oceanParamList)
{}

Ocean::Ocean(RCP<Epetra_Comm> Comm, Teuchos::ParameterList& oceanParamList)
    :
    params_                ("Ocean Configuration"),
    // Create THCM object
    //  THCM is implemented as a Singleton, which allows only a single
    //  instance at a time. The Ocean class can access THCM with a call
    //  to THCM::Instance()
    thcm_                  (new THCM(oceanParamList.sublist("THCM"), Comm)),
    solverInitialized_     (false),  // Solver needs initialization
    precInitialized_       (false),  // Preconditioner needs initialization
    recompPreconditioner_  (true),   // We need a preconditioner to start with
    recompMassMat_         (true)    // We need a mass matrix to start with
{
    INFO("Ocean: constructor...");

    oceanParamList.validateParametersAndSetDefaults(getDefaultInitParameters());
    params_.setParameters(oceanParamList);

    // inherited input/output datamembers
    inputFile_   = params_.get<std::string>("Input file");
    outputFile_  = params_.get<std::string>("Output file");
    saveMask_    = params_.get<bool>("Save mask");
    loadMask_    = params_.get<bool>("Load mask");

    loadState_   = params_.get<bool>("Load state");
    saveState_   = params_.get<bool>("Save state");
    saveEvery_   = params_.get<int>("Save frequency");

    loadSalinityFlux_    = params_.get<bool>("Load salinity flux");
    saveSalinityFlux_    = params_.get<bool>("Save salinity flux");
    loadTemperatureFlux_ = params_.get<bool>("Load temperature flux");
    saveTemperatureFlux_ = params_.get<bool>("Save temperature flux");

    useFort3_            = params_.get<bool>("Use legacy fort.3 output");
    useFort44_           = params_.get<bool>("Use legacy fort.44 output");
    saveColumnIntegral_  = params_.get<bool>("Save column integral");
    maxMaskFixes_        = params_.get<int>("Max mask fixes");

    analyzeJacobian_     = params_.get<bool>("Analyze Jacobian");

    // initialize postprocessing counter
    ppCtr_ = 0;

    // set the communicator object
    comm_ = Comm;

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

    // If specified we load a pre-existing state and parameters (x,l)
    // thereby overwriting the starting parameters
    // This will be able to load salinity and temperature fluxes as well.

    state_->PutScalar(0.0);
    if (loadState_ || loadSalinityFlux_ || loadTemperatureFlux_)
        loadStateFromFile(inputFile_);

    // make sure initial state satisfies integral condition
    if (THCM::Instance().getSRES() == 0)
    {
        THCM::Instance().setIntCondCorrection(state_);
    }

    // Now that we have the state and parameters initialize
    // the Jacobian, solution and rhs
    initializeOcean();

    // Obtain landmask, in case this is a fresh start we probably need
    // to make some adjustments in the landmask
    bool adjustMask = (loadState_ && loadMask_) ? false : true;
    landmask_ = getLandMask("current", adjustMask);

    // Initialize preconditioner
    initializePreconditioner();

    // Inspect current state
    inspectVector(state_);

    //------------------------------------------------------------------
    // Create surface temperature and salinity restrict/import strategies
    //------------------------------------------------------------------
    // Create lists of surface indices
    std::vector<int> tRows, sRows;

    for (int j = 0; j != M_; ++j)
        for (int i = 0; i != N_; ++i)
        {
            tRows.push_back(FIND_ROW2(_NUN_, N_, M_, L_,i,j,L_-1,TT));
            sRows.push_back(FIND_ROW2(_NUN_, N_, M_, L_,i,j,L_-1,SS));
        }

    // Create restricted maps
    tIndexMap_ = Utils::CreateSubMap(state_->Map(), tRows);
    sIndexMap_ = Utils::CreateSubMap(state_->Map(), sRows);

    // Create SST vector
    sst_  = Teuchos::rcp(new Epetra_Vector(*tIndexMap_));

    // Create SSS vector
    sss_  = Teuchos::rcp(new Epetra_Vector(*sIndexMap_));

    // Create Qsi vector
    Qsi_  = Teuchos::rcp(new Epetra_Vector(*domain_->GetStandardSurfaceMap()));

    // Create Msi vector
    Msi_  = Teuchos::rcp(new Epetra_Vector(*domain_->GetStandardSurfaceMap()));

    // Create Gsi vector
    Gsi_  = Teuchos::rcp(new Epetra_Vector(*domain_->GetStandardSurfaceMap()));

    // Create import strategies
    // Target map: IndexMap
    // Source map: state_->Map()
    surfaceTimporter_ =
        Teuchos::rcp(new Epetra_Import(*tIndexMap_, state_->Map()));

    surfaceSimporter_ =
        Teuchos::rcp(new Epetra_Import(*sIndexMap_, state_->Map()));

    oceanParamList = params_;
    INFO(oceanParamList);
    INFO("\n");
    INFO("Ocean couplings: coupled_T = " << getCoupledT() );
    INFO("                 coupled_S = " << getCoupledS() );
    INFO("\n");
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
    // Initialize solution and rhs
    sol_ = rcp(new Epetra_Vector(*domain_->GetStandardMap(), true));
    rhs_ = rcp(new Epetra_Vector(*domain_->GetStandardMap(), true));

    // Obtain Jacobian from THCM
    THCM::Instance().evaluate(*state_, Teuchos::null, true);
    jac_ = THCM::Instance().getJacobian();

    INFO("Ocean: Obtained Jacobian from THCM");

    // Obtain mass matrix B from THCM. Note that we assume the mass
    // matrix is independent of the state and parameters so we only
    // compute it when asked for through recompMassMat_ flag.
    computeMassMat();

    // Compute right hand side and print its norm
    computeRHS();

    INFO("Ocean: initialization: ||F|| = " << Utils::norm(rhs_));
    // Export some diagnostics
    // Utils::save(state_, "initialstate");
    // Utils::save(rhs_,   "initialrhs");

    // // Perturb model
    // *sol_ = *state_;
    // double par = getPar();
    // setPar(par + 1e-8);
    // state_->PutScalar(1.0);

    // computeRHS();
    // Utils::save(rhs_, "perturbedrhs");

    // // Restore model
    // setPar(par);
    // *state_ = *sol_;
    // sol_->PutScalar(0.0);

    // Get the rowmap for the pressure points
    Teuchos::RCP<Epetra_Map> RowMap = domain_->GetSolveMap();

    mapP_     = Utils::CreateSubMap(*RowMap, _NUN_, PP);
    mapU_     = Utils::CreateSubMap(*RowMap, _NUN_, UU);

    singRows_ = Teuchos::rcp(new Epetra_Vector( *mapP_ ) );
}

//====================================================================
int Ocean::analyzeJacobian1()
{
    if (!analyzeJacobian_)
        return 0;

    INFO("\n  <><>  Analyze Jacobian P rows...");

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
            std::cout << "  <><>  problem P row: " << row << std::endl;
        }
    }

    int sumFoundP = 0;
    comm_->SumAll(&singRowsFound, &sumFoundP, 1);
    INFO("  <><>  problem P rows found  (local): " << singRowsFound);
    INFO("  <><>  problem P rows found (global): " << sumFoundP);

    return sumFoundP;
}

//==================================================================
int Ocean::analyzeJacobian2()
{
    if (!analyzeJacobian_)
        return 0;

    INFO("\n  <><>  Analyze Jacobian S column integrals...\n");
    // Another approach to analyze the Jacobian is through the volume
    // integrals of its columns. This only makes sense when the volume
    // integral of the forcing is zero and we temporarily disable the
    // integral condition. To get a meaningful test state we already
    // use the solver so the mask-fix cycle using analyzeJacobian1
    // should have been performed.

    // Copy the original Jacobian and mass matrix from THCM
    Teuchos::RCP<Epetra_CrsMatrix> tmpJac =
        Teuchos::rcp(new Epetra_CrsMatrix(*THCM::Instance().getJacobian()));
    Teuchos::RCP<Epetra_Vector> tmpB =
        Teuchos::rcp(new Epetra_Vector(*THCM::Instance().DiagB()));

    // Create converged test vector such that it satisfies boundary
    // conditions.
    Teuchos::RCP<Epetra_Vector> testvec = initialState();

    // // exaggerate a little
    // testvec->Scale(1e2);

    // Compute test Jacobian and mass matrix
    THCM::Instance().evaluate(*testvec, Teuchos::null, true, true);

    // Copy the test Jacobian from THCM
    Teuchos::RCP<Epetra_CrsMatrix> mat =
        Teuchos::rcp(new Epetra_CrsMatrix(*THCM::Instance().getJacobian()));

    // DUMPMATLAB("ocean_jac", *mat);
    // DUMP_VECTOR("intcond_coeff", *getIntCondCoeff());
    // DUMP_VECTOR("testvec", *testvec);

    // Restore the original Jacobian and mass matrix in THCM
    Teuchos::RCP<Epetra_CrsMatrix> jac = THCM::Instance().getJacobian();
    *jac = *tmpJac;

    Teuchos::RCP<Epetra_Vector> diagB = THCM::Instance().DiagB();
    *diagB = *tmpB;

    // Compute column integrals for the salinity block
    Teuchos::RCP<Epetra_Vector> colInts = getColumnIntegral(mat, false);

    Teuchos::RCP<Epetra_Map> rowMap = domain_->GetSolveMap();
    Teuchos::RCP<Epetra_Map> mapS =
        Utils::CreateSubMap(*rowMap, _NUN_, SS);

    Teuchos::RCP<Epetra_Import> importS =
        Teuchos::rcp(new Epetra_Import(*rowMap, *mapS) );

    Epetra_Vector ints(*mapS);
    ints.Export(*colInts, *importS, Zero);

    int dim = mapS->NumMyElements();
    assert(dim == ints.MyLength());

    // At this point we use values in singrows to identify bad points
    int badSintsfound = 0;
    for (int i = 0; i != dim; ++i)
        if (std::abs(ints[i]) > 1e-6)
        {
            (*singRows_)[i] = 2;
            badSintsfound++;
            std::cout << "  <><>  nonzero column integral, "
                      << "  GID = " << ints.Map().GID(i)
                      << " value = " << std::abs(ints[i]) << std::endl;
        }

    int sumFoundS = 0;
    comm_->SumAll(&badSintsfound, &sumFoundS, 1);

    INFO("\n  <><>  nonzero column integrals found: " << sumFoundS << '\n');

    //  DUMP_VECTOR("intcond_coeff", ints);

    return sumFoundS;
}

//==================================================================
// --> This can be more efficient
void Ocean::inspectVector(Teuchos::RCP<Epetra_Vector> x)
{
    INFO("Ocean: inspect vector...");

    // for now we just check whether the surface w-values are zero, if
    // not we put them to zero

    // this setup lets us choose more rows in a layer to reset so now
    // we could put the surface u,v,p anomalies to zero as well, but
    // that is evil

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
Utils::MaskStruct Ocean::getLandMask(bool adjustMask)
{
    return getLandMask("current", adjustMask);
}

//====================================================================
Utils::MaskStruct Ocean::getLandMask(std::string const &fname, bool adjustMask)
{
    Utils::MaskStruct mask;

    // Load the landmask fname
    mask.local = THCM::Instance().getLandMask(fname);
    THCM::Instance().setLandMask(mask.local);
    THCM::Instance().evaluate(*state_, Teuchos::null, true);

    if (adjustMask) // FIXME the whole analyzeJacobian stuff should be
                    // part of THCM such that we can postpone the
                    // creation of the matrix graph. -> The integral
                    // condition row needs to be chosen at an ocean
                    // point.
    {
        // zero the singRows array
        singRows_->PutScalar(0.0);

        int badProws = 1, badSints = 1;
        for (int i = 0; i != maxMaskFixes_; ++i)
        {
            for (int j = 0; j != maxMaskFixes_; ++j)
            {
                badProws = analyzeJacobian1();
                if (badProws == 0)
                    break;

                // If we find singular pressure rows we adjust the current landmask
                mask.local = THCM::Instance().getLandMask("current", singRows_);

                //  Putting a fixed version of the landmask back in THCM
                THCM::Instance().setLandMask(mask.local);

                // Perform a Newton iteration to get a physical state before
                // repeating the analysis.
                THCM::Instance().evaluate(*state_, Teuchos::null, true);

                // This adds the possibility of bad S integrals so we
                // increase this counter
                badSints++;
            }

            if ( (badSints + badProws) == 0 )
                break;

            for (int j = 0; j != maxMaskFixes_; ++j)
            {
                badSints = analyzeJacobian2();
                if ( badSints == 0 )
                    break;

                // If we find singular pressure rows we adjust the current landmask
                mask.local = THCM::Instance().getLandMask("current", singRows_);

                //  Putting a fixed version of the landmask back in THCM
                THCM::Instance().setLandMask(mask.local);

                // Perform a Newton iteration to get a physical state before
                // repeating the analysis.
                THCM::Instance().evaluate(*state_, Teuchos::null, true);

                // This adds the possibility of bad P rows so we
                // increase this counter
                badProws++;
            }

            if ( (badSints + badProws) == 0 )
                break;
        }

        // The solver of the ocean model should be reinitialized when
        // analyzeJacobian2() is done.
        solverInitialized_ = false;
    }

    // Get the current global landmask from THCM.
    mask.global = THCM::Instance().getLandMask();

    // Copy to full global mask tmp
    std::vector<int> tmp(*mask.global);

    // Create global borderless mask rcp
    mask.global_borderless = std::make_shared<std::vector<int> >();

    for (int k = 1; k != L_+1; ++k)
        for (int j = 1; j != M_+1; ++j)
            for (int i = 1; i != N_+1; ++i)
            {
                mask.global_borderless->push_back(tmp[k*(M_+2)*(N_+2) + j*(N_+2) + i]);
            }

    assert( (int) mask.global_borderless->size() == N_*M_*L_ );


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

    // Set global dimensions
    mask.n = N_;
    mask.m = M_;
    mask.l = L_;

    // Return the struct
    return mask;
}

//===================================================================
void Ocean::setLandMask(Utils::MaskStruct const &mask, bool global)
{
    INFO("Ocean: set landmask " << mask.label << "...");
    THCM::Instance().setLandMask(mask.local);

    if (global)
        THCM::Instance().setLandMask(mask.global);

    currentMask_ = mask.label;
    INFO("Ocean: set landmask " << mask.label << "... done");
}

//==================================================================
void Ocean::applyLandMask(Utils::MaskStruct mask, double factor)
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
                          Utils::MaskStruct mask, double factor)
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
                          Utils::MaskStruct maskA,
                          Utils::MaskStruct maskB)
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

    // Output legacy datafiles
    printLegacyFiles();
}

//====================================================================
void Ocean::postProcess()
{
    // increase postprocessing counter
    ppCtr_++;

    TIMER_START("Ocean: saveStateToFile");
    if (saveState_)
        saveStateToFile(outputFile_); // Save to hdf5
    TIMER_STOP("Ocean: saveStateToFile");

    if ((saveEvery_ > 0) && (ppCtr_ % saveEvery_) == 0)
    {
        std::stringstream append;
        append << "." << ppCtr_;
        copyState(append.str());
    }

    printLegacyFiles(); // Print in standard fortran format

    // Column integral can be used to check discretization: should be
    // zero, excluding integral condition and its dependencies.
    if (saveColumnIntegral_) // Compute and save column integral
        Utils::save(getColumnIntegral(jac_), "columnIntegral");
}

//=====================================================================
std::string Ocean::writeData(bool describe) const
{
    std::ostringstream datastring;
    if (describe)
    {
        if (solverInitialized_)
        {
            datastring << std::setw(_FIELDWIDTH_/3)
                       << "MV";
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
        double psiMin, psiMax;
        getPsiM(psiMin, psiMax);

        if (solverInitialized_)
        {
            datastring << std::setw(_FIELDWIDTH_/3)
                       << std::round(effort_);
            effortCtr_ = 0;
        }

        datastring << std::scientific << std::setw(_FIELDWIDTH_)
                   << psiMax
                   << std::setw(_FIELDWIDTH_)
                   << psiMin;

        return datastring.str();
    }
}

int Ocean::getPsiM(double &psiMin, double &psiMax) const
{
    grid_->ImportData(*state_);
    psiMax = grid_->psimMax();
    psiMin = grid_->psimMin();

    double r0dim, udim, hdim;
    FNAME(get_parameters)(&r0dim, &udim, &hdim);

    const double transc = r0dim * hdim * udim;
    psiMax *= transc * 1e-6; // conversion to Sv
    psiMin *= transc * 1e-6; //

    return 0;
}

//==================================================================
int Ocean::getCoupledT()
{
    return THCM::Instance().getCoupledT();
}

//==================================================================
int Ocean::getCoupledS()
{
    return THCM::Instance().getCoupledS();
}

//==================================================================
double Ocean::getSCorr()
{
    return THCM::Instance().getSCorr();
}

//==================================================================
int Ocean::dof() { return _NUN_; }

//==================================================================
int Ocean::interface_row(int i, int j, int XX)
{ return FIND_ROW2(_NUN_, N_, M_, L_, i, j, L_-1, XX); }

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

    // Initialize the preconditioner
    if (!precInitialized_)
        initializePreconditioner();

    // Initialize the solver
    initializeBelos();

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

    Teuchos::ParameterList &belosParams = params_.sublist("Belos Solver");

    // A few FGMRES parameters are made available in solver_params.xml:
    int gmresIters  = belosParams.get<int>("FGMRES iterations");
    double gmresTol = belosParams.get<double>("FGMRES tolerance");
    int maxrestarts = belosParams.get<int>("FGMRES restarts");
    int output      = belosParams.get<int>("FGMRES output");
    bool testExpl   = belosParams.get<bool>("FGMRES explicit residual test");

    int NumGlobalElements = state_->GlobalLength();
    int blocksize         = 1; // number of vectors in rhs
    int maxiters          = NumGlobalElements/blocksize - 1;

    // Create Belos parameterlist
    RCP<Teuchos::ParameterList> belosParamList = rcp(new Teuchos::ParameterList("Belos List"));
    belosParamList->set("Block Size", blocksize);
    belosParamList->set("Flexible Gmres", true);
    belosParamList->set("Adaptive Block Size", true);
    belosParamList->set("Num Blocks", gmresIters);
    belosParamList->set("Maximum Restarts", maxrestarts);
    belosParamList->set("Orthogonalization","DGKS");
    belosParamList->set("Output Frequency", output);
    belosParamList->set("Verbosity", Belos::Errors + Belos::Warnings);
    belosParamList->set("Maximum Iterations", maxiters);
    belosParamList->set("Convergence Tolerance", gmresTol);
    belosParamList->set("Explicit Residual Test", testExpl);
    belosParamList->set("Implicit Residual Scaling",
                         "Norm of Preconditioned Initial Residual");

    // belosParamList->set("Implicit Residual Scaling", "Norm of RHS");
    // belosParamList->set("Implicit Residual Scaling", "Norm of Initial Residual");
    // belosParamList->set("Explicit Residual Scaling", "Norm of RHS");

    // Belos block FGMRES setup
    belosSolver_ =
        rcp(new Belos::BlockGmresSolMgr
            <double, Epetra_MultiVector, Epetra_Operator>
            (problem_, belosParamList));

    // initialize effort counter
    effortCtr_ = 0;
    effort_ = 0.0;

}

//=====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::initialState()
{
    // Initialize result
    Teuchos::RCP<Epetra_Vector> result =
        Teuchos::rcp(new Epetra_Vector(*state_));

    // Save state
    Teuchos::RCP<Epetra_Vector> tmp =
        Teuchos::rcp(new Epetra_Vector(*state_));

    // Do a few Newton steps to get a physical state
    double par = getPar("Combined Forcing");
    setPar("Combined Forcing", 1e-8); // perturb parameter
    // start from trivial solution
    state_->PutScalar(0.0);
    computeRHS();
    INFO("Initial Newton iteration, norm rhs = " << Utils::norm(rhs_));
    for (int k = 0; k != 1; ++k)
    {
        computeJacobian();
        rhs_->Scale(-1.0);
        solve(rhs_);
        state_->Update(1.0, *sol_, 1.0);
        computeRHS();
        INFO("Initial Newton iteration, norm rhs = " << Utils::norm(rhs_));
    }
    *result = *state_;

    // Restore parameter, sol and state
    setPar("Combined Forcing", par);
    *state_ = *tmp;
    sol_->PutScalar(0.0);

    return result;
}

//=====================================================================
void Ocean::solve(Teuchos::RCP<const Epetra_MultiVector> rhs)
{
    // Check whether solver is initialized, if not perform the
    // initialization here
    if (!solverInitialized_)
        initializeSolver();

    // Get new preconditioner
    buildPreconditioner();

    // Use trivial initial solution
    sol_->PutScalar(0.0);

    // Set right hand side
    Teuchos::RCP<const Epetra_MultiVector> b;
    if (rhs == Teuchos::null)
        b = rhs_;
    else
        b = rhs;

    bool set = problem_->setProblem(sol_, b);

    TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
                               "*** Belos::LinearProblem failed to setup");

    // ---------------------------------------------------------------------
    // Start solving J*x = F, where J = jac_, x = sol_ and F = rhs
    TIMER_START("Ocean: solve...");
    INFO("Ocean: solve...");
    INFO(" |x| = " << Utils::norm(sol_));
    INFO(" |b| = " << Utils::norm(rhs));

    int    iters;
    double tol;
    try
    {
        belosSolver_->solve();      // Solve
    }
    catch (std::exception const &e)
    {
        ERROR("Ocean: exception caught: " << e.what(), __FILE__, __LINE__);
    }

    INFO("Ocean: solve... done");
    TIMER_STOP("Ocean: solve...");

    // ---------------------------------------------------------------------
    // Inspect solve and update effort
    iters = belosSolver_->getNumIters();
    tol   = belosSolver_->achievedTol();
    INFO("Ocean: FGMRES, i = " << iters << ", ||r|| = " << tol);

    // keep track of effort
    if (effortCtr_ == 0)
        effort_ = 0;

    effortCtr_++;
    effort_ = (effort_ * (effortCtr_ - 1) + iters ) / effortCtr_;

    Teuchos::RCP<Epetra_Vector> bvec =
        Teuchos::rcp(new Epetra_Vector(*(*rhs)(0)));

    double normb = Utils::norm(bvec);
    double nrm   = explicitResNorm(bvec);

    INFO("           ||b||         = " << normb);
    INFO("           ||x||         = " << Utils::norm(sol_));
    INFO("        ||b-Ax|| / ||b|| = " << nrm / normb);

    if ((tol > 0) && (normb > 0) && ( (nrm / normb / tol) > 10))
    {
        WARNING("Actual residual norm at least ten times larger: "
                << (nrm / normb) << " > " << tol
                , __FILE__, __LINE__);
    }

    TRACK_ITERATIONS("Ocean: FGMRES iterations...", iters);
}

//=====================================================================
double Ocean::explicitResNorm(VectorPtr rhs)
{
    RCP<Epetra_Vector> Ax =
        rcp(new Epetra_Vector(*(domain_->GetSolveMap())));
    jac_->Apply(*sol_, *Ax);        // A*x
    Ax->Update(1.0, *rhs, -1.0);    // b - A*x
    double nrm;
    Ax->Norm2(&nrm);                // nrm = ||b-A*x||
    // Utils::save(Ax, "lsresidual");
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
void Ocean::computeForcing()
{
    // evaluate rhs in THCM with the current state
    TIMER_START("Ocean: compute Frc...");
    THCM::Instance().computeForcing();
    frc_ = THCM::Instance().getForcing();
    TIMER_STOP("Ocean: compute Frc...");
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
Teuchos::RCP<Epetra_Vector> Ocean::getSolution(char mode)
{
    return Utils::getVector(mode, sol_);
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getState(char mode)
{
    return Utils::getVector(mode, state_);
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getRHS(char mode)
{
    return Utils::getVector(mode, rhs_);
}

//====================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getMassMat(char mode)
{
    diagB_ = THCM::Instance().DiagB();
    return Utils::getVector(mode, diagB_);
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

    return Utils::getVector(mode, vecM);
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

    // check matrix residual
    // Teuchos::RCP<Epetra_MultiVector> r =
    //     Teuchos::rcp(new Epetra_MultiVector(v));;

    // applyMatrix(out, *r);
    // r->Update(1.0, v, -1.0);
    // double rnorm = Utils::norm(r);

    // INFO("Ocean: preconditioner residual: " << rnorm);

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
void Ocean::computeMassMat()
{
    if (recompMassMat_)
    {
        THCM::Instance().evaluateB();
    }
    recompMassMat_ = false; // Disable subsequent recomputes
}

//====================================================================
void Ocean::applyMassMat(Epetra_MultiVector const &v, Epetra_MultiVector &out)
{
    // Compute mass matrix
    computeMassMat();

    diagB_ = THCM::Instance().DiagB();

    // element-wise multiplication (out = 0.0*out + 1.0*B*v)
    out.Multiply(1.0, *diagB_, v, 0.0);
}

//====================================================================
void Ocean::synchronize(std::shared_ptr<Atmosphere> atmos)
{
    TIMER_START("Ocean: set atmosphere...");

    // Obtain and set atmosphere T at the interface
    Teuchos::RCP<Epetra_Vector> atmosT  = atmos->interfaceT();
    THCM::Instance().setAtmosphereT(atmosT);

    // Obtain and set humidity field at the interface
    Teuchos::RCP<Epetra_Vector> atmosQ  = atmos->interfaceQ();
    THCM::Instance().setAtmosphereQ(atmosQ);

    // Obtain and set albedo field at the interface
    Teuchos::RCP<Epetra_Vector> atmosA  = atmos->interfaceA();
    THCM::Instance().setAtmosphereA(atmosA);

    // Obtain and set precipitation field at the interface
    Teuchos::RCP<Epetra_Vector> atmosP  = atmos->interfaceP();
    THCM::Instance().setAtmosphereP(atmosP);

    // We also need to know a few atmospheric parameters to compute E,
    // P and their derivatives w.r.t. SST (To) and humidity (q) These
    // may depend on continuation parameters, so the call belongs
    // here.
    Atmosphere::CommPars atmosPars;
    atmos->getCommPars(atmosPars);
    FNAME( set_atmos_parameters )( &atmosPars );

    TIMER_STOP("Ocean: set atmosphere...");
}

//====================================================================
void Ocean::synchronize(std::shared_ptr<SeaIce> seaice)
{
    TIMER_START("Ocean: set seaice...");
    Qsi_ = seaice->interfaceQ();
    THCM::Instance().setSeaIceQ(Qsi_);

    Msi_ = seaice->interfaceM();
    THCM::Instance().setSeaIceM(Msi_);

    Gsi_ = seaice->interfaceG();
    THCM::Instance().setSeaIceG(Gsi_);

    SeaIce::CommPars seaicePars;
    seaice->getCommPars(seaicePars);
        
    FNAME( set_seaice_parameters )( &seaicePars );

    TIMER_STOP("Ocean: set seaice...");
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
Teuchos::RCP<Epetra_Vector> Ocean::interfaceE()
{
    return THCM::Instance().getOceanE();
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getSunO()
{
    return THCM::Instance().getSunO();
}

//==================================================================
int Ocean::getRowIntCon()
{
    return THCM::Instance().getRowIntCon();
}

//==================================================================
std::shared_ptr<Utils::CRSMat> Ocean::getBlock(std::shared_ptr<Atmosphere> atmos)
{
    // initialize empty CRS matrix
    std::shared_ptr<Utils::CRSMat> block = std::make_shared<Utils::CRSMat>();

    // get parameter dependencies
    double Ooa, Os, nus, eta, lvsc, qdim, pQSnd;
    FNAME(getdeps)(&Ooa, &Os, &nus, &eta, &lvsc, &qdim, &pQSnd);
    Atmosphere::CommPars atmosPars;
    atmos->getCommPars(atmosPars);
    double albed = atmosPars.da;

    int T = ATMOS_TT_; // (1-based) atmos temperature: first unknown
    int Q = ATMOS_QQ_; // (1-based) atmos humidity: second unknown
    int A = ATMOS_AA_; // (1-based) atmos albedo: third unknown
    int P = ATMOS_PP_; // (1-based) atmos global precipitation: auxiliary

    int rowIntCon = THCM::Instance().getRowIntCon();

    // FIXME if this block would be computed locally we would not need
    // an allgather
    Teuchos::RCP<Epetra_MultiVector> Msi  = Utils::AllGather(*Msi_);

    // FIXME: factorize as this is (probably) constant
    Teuchos::RCP<Epetra_MultiVector> Pdist =
        Utils::AllGather(*atmos->getPdist());

    // Obtain shortwave radiative heat (global) field --> FIXME
    // factorize as this is constant
    Teuchos::RCP<Epetra_MultiVector> suno =
        Utils::AllGather(*THCM::Instance().getSunO());

    // fill CRS struct
    int el_ctr = 0;
    int col;
    int sr;
    double M; // sea ice mask value
    double S; // shortwave radiative flux dependency
    double dTFT; // d / dtatm (F_T)
    double dQFT; // d / dqatm (F_T)
    double dQFS; // d / dqatm (F_S)
    double dPFS; // d / dpatm (F_S)
    double dAFT; // d / dalbe (F_T)

    double comb = getPar("Combined Forcing");
    double sunp = getPar("Solar Forcing");
    double Pd;

    for (int k = 0; k != L_; ++k)
        for (int j = 0; j != M_; ++j)
            for (int i = 0; i != N_; ++i)
            {
                // surface row
                sr = j*N_+i;

                // sea ice mask value
                M  = (*(*Msi)(0))[sr];

                // shortwave distribution
                S  = (*(*suno)(0))[sr];

                // precipitation distribution
                Pd = (*(*Pdist)(0))[sr];

                for (int xx = UU; xx <= SS; ++xx)
                {
                    block->beg.push_back(el_ctr);

                    if ( (k == L_-1) &&
                         ( (*landmask_.global_surface)[sr] == 0 ) )
                    {
                        // surface T row
                        if ( (xx == TT) && getCoupledT() )
                        {
                            // tatm dependency
                            dTFT = Ooa * (1.0 - M);
                            // negating as the Jacobian is taken negative
                            block->co.push_back( -dTFT );
                            block->jco.push_back(atmos->interface_row(i,j,T) );
                            el_ctr++;

                            // albe dependency
                            dAFT = -comb * sunp * S * albed * (1.0 - M);
                            // negating as the Jacobian is taken negative
                            block->co.push_back( -dAFT );
                            block->jco.push_back(atmos->interface_row(i,j,A) );
                            el_ctr++;

                            // qatm dependency
                            dQFT = lvsc * eta * qdim * (1.0 - M);
                            // negating as the Jacobian is taken negative
                            block->co.push_back(-dQFT);
                            block->jco.push_back(atmos->interface_row(i,j,Q) );
                            el_ctr++;
                        }

                        // surface S row, exclude integral condition row
                        else if ((xx == SS) && getCoupledS() &&
                                 FIND_ROW2(_NUN_, N_, M_, L_, i, j, k, xx) != rowIntCon)
                        {
                            // humidity dependency
                            dQFS = -nus * (1.0 - M);
                            block->co.push_back(-dQFS);
                            block->jco.push_back(atmos->interface_row(i,j,Q) );
                            el_ctr++;

                            // Precipitation dependency. The
                            // derivative is taken with respect to the
                            // P anomaly, not to the full dimensional
                            // P with spatial distribution
                            col = atmos->interface_row(i,j,P);
                            if (col >= 0)
                            {
                                dPFS = -nus * Pd * (1.0 - M);
                                block->co.push_back(-dPFS);
                                block->jco.push_back(col);
                                el_ctr++;
                            }
                        }
                    }
                }
            }

    // final entry in beg ( == nnz)
    block->beg.push_back(el_ctr);

    assert( (int) block->co.size() == block->beg.back());

    return block;
}

//==================================================================
std::shared_ptr<Utils::CRSMat> Ocean::getBlock(std::shared_ptr<SeaIce> seaice)
{
    // initialize empty CRS matrix
    std::shared_ptr<Utils::CRSMat> block = std::make_shared<Utils::CRSMat>();
    int rowIntCon = THCM::Instance().getRowIntCon();

    //FIXME Gathers are unnecessary if we compute this block locally,
    //preferably in the fortran code.
    THCM::Derivatives d = THCM::Instance().getDerivatives();
    Teuchos::RCP<Epetra_MultiVector> dFTdM = Utils::AllGather(*d.dFTdM);
    Teuchos::RCP<Epetra_MultiVector> dFSdQ = Utils::AllGather(*d.dFSdQ);
    Teuchos::RCP<Epetra_MultiVector> dFSdM = Utils::AllGather(*d.dFSdM);
    Teuchos::RCP<Epetra_MultiVector> dFSdG = Utils::AllGather(*d.dFSdG);

    int el_ctr = 0;
    int sr; // surface row

    int seaiceQQ = SEAICE_QQ_; // (1-based) heat flux unknown in the sea ice model
    int seaiceMM = SEAICE_MM_; // (1-based) mask unknown in the sea ice model
    int seaiceGG = SEAICE_GG_; // (1-based) auxiliary correction in the sea ice model

    double dFTdMval;
    double dFSdQval;
    double dFSdMval;
    double dFSdGval;

    for (int k = 0; k != L_; ++k)
        for (int j = 0; j != M_; ++j)
            for (int i = 0; i != N_; ++i)
            {
                sr       = j*N_+i; // surface row
                dFTdMval = (*(*dFTdM)(0))[sr];
                dFSdQval = (*(*dFSdQ)(0))[sr];
                dFSdMval = (*(*dFSdM)(0))[sr];
                dFSdGval = (*(*dFSdG)(0))[sr];

                // surface, non-land point
                for (int XX = UU; XX <= SS; ++XX)
                {
                    block->beg.push_back(el_ctr);
                    if ( ( k == L_-1 ) &&
                         ( (*landmask_.global_surface)[sr] == 0 ))
                    {
                        // surface T row
                        if ( (XX == TT) && getCoupledT() )
                        {
                            block->co.push_back( -dFTdMval );
                            block->jco.push_back(seaice->interface_row(i,j,seaiceMM));
                            el_ctr++;
                        }
                        // surface S row, exclude integral condition row
                        else if ((XX == SS) && getCoupledS() &&
                                 FIND_ROW2(_NUN_, N_, M_, L_, i, j, k, XX) != rowIntCon)
                        {
                            block->co.push_back( -dFSdQval );
                            block->jco.push_back(seaice->interface_row(i,j,seaiceQQ));
                            el_ctr++;

                            block->co.push_back( -dFSdMval );
                            block->jco.push_back(seaice->interface_row(i,j,seaiceMM));
                            el_ctr++;

                            block->co.push_back( -dFSdGval );
                            block->jco.push_back(seaice->interface_row(i,j,seaiceGG));
                            el_ctr++;
                        }
                    }
                }
            }

    block->beg.push_back(el_ctr);
    assert( (int) block->co.size() == block->beg.back());

    return block;
}

//====================================================================
// Fill and return a copy of the surface temperature
Teuchos::RCP<Epetra_Vector> Ocean::interfaceT()
{
    TIMER_START("Ocean: get surface temperature...");
    CHECK_MAP(sst_, tIndexMap_);
    CHECK_ZERO(sst_->Import(*state_, *surfaceTimporter_, Insert));
    TIMER_STOP("Ocean: get surface temperature...");
    return sst_;
}

//====================================================================
// Fill and return a copy of the surface salinity
Teuchos::RCP<Epetra_Vector> Ocean::interfaceS()
{
    TIMER_START("Ocean: get surface salinity...");
    CHECK_MAP(sss_, sIndexMap_);
    CHECK_ZERO(sss_->Import(*state_, *surfaceSimporter_, Insert));
    TIMER_STOP("Ocean: get surface salinity...");
    return sss_;
}

//=====================================================================
void Ocean::printLegacyFiles()
{
    // This writes to fort.44. If requested we can also output the
    // legacy fort.3 here. Default behaviour: exit write function
    // after creating fort.44.

    if (!useFort44_)
        return;

    TIMER_START("Ocean: printLegacyFiles");

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

        if (useFort3_)
        {
            INFO( "Writing to fort." << filename
                  << " at label " << label);

            (*solution)(0)->ExtractCopy(solutionArray);
            (*rhs)(0)->ExtractCopy(rhsArray);
        }

        FNAME(write_data)(solutionArray, &filename, &label);
    }

    delete [] solutionArray;
    delete [] rhsArray;

    TIMER_STOP("Ocean: printLegacyFiles");
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

void Ocean::dumpBlocks()
{
    DUMPMATLAB("ocean_jac", *jac_);
}

//==================================================================
void Ocean::integralChecks(Teuchos::RCP<Epetra_Vector> state,
                           double &salt_advection,
                           double &salt_diffusion)
{
    THCM::Instance().integralChecks(state,
                                    salt_advection,
                                    salt_diffusion);
}

//==================================================================
Teuchos::RCP<Epetra_Vector>
Ocean::getColumnIntegral(Teuchos::RCP<Epetra_CrsMatrix> mat,
                         bool useSRES)
{
    TIMER_START("Column integral");
    INFO("Ocean: computing column volume integrals of Jacobian");

    // Copy into tmp
    Teuchos::RCP<Epetra_CrsMatrix> tmp_mat =
        Teuchos::rcp(new Epetra_CrsMatrix(*mat));

    // Rowscaling of the matrix with integral coefficients
    Teuchos::RCP<Epetra_Vector> icCoef =
        Teuchos::rcp(new Epetra_Vector(*getIntCondCoeff()));

    // Ignore the integral condition row if needed
    int sres = THCM::Instance().getSRES();
    if ( ( sres == 0 ) && useSRES )
    {
        int rowIntCon = getRowIntCon();
        int lidIntCon = icCoef->Map().LID(rowIntCon);
        if (lidIntCon >= 0)
            (*icCoef)[lidIntCon] = 0;
    }

    tmp_mat->LeftScale(*icCoef);

    // Change integral coefficients into column selectors
    for (int i = SS-1; i < icCoef->MyLength(); i+=_NUN_)
    {
        (*icCoef)[i] = 1;
    }

    tmp_mat->RightScale(*icCoef);

    // Create vector that will contain the column integrals
    Teuchos::RCP<Epetra_Vector> sums =
        Teuchos::rcp(new Epetra_Vector(state_->Map()));

    // Integrate the columns of mat  into sums
    Utils::colSums(*tmp_mat, *sums);

    TIMER_STOP("Column integral");
    return sums;
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Ocean::getIntCondCoeff()
{
    return THCM::Instance().getIntCondCoeff();
}

//=====================================================================
void Ocean::additionalExports(EpetraExt::HDF5 &HDF5, std::string const &filename)
{
    TIMER_START("Ocean: additionalExports");
    std::vector<Teuchos::RCP<Epetra_Vector> > fluxes =
        THCM::Instance().getFluxes();

    if (saveSalinityFlux_)
    {
        // Write emip to ocean output file
        INFO("Writing salinity fluxes to " << filename);
        HDF5.Write("SalinityFlux",       *fluxes[ THCM::_Sal  ]);
        HDF5.Write("OceanAtmosSalFlux",  *fluxes[ THCM::_QSOA ]);
        HDF5.Write("OceanSeaIceSalFlux", *fluxes[ THCM::_QSOS ]);
    }

    if (saveTemperatureFlux_)
    {
        // Write heat fluxes to ocean output file. In the coupled case
        // these are dimensional.
        INFO("Writing temperature fluxes to " << filename);
        HDF5.Write("TemperatureFlux",  *fluxes[ THCM::_Temp ]);
        HDF5.Write("ShortwaveFlux",    *fluxes[ THCM::_QSW  ]);
        HDF5.Write("SensibleHeatFlux", *fluxes[ THCM::_QSH  ]);
        HDF5.Write("LatentHeatFlux",   *fluxes[ THCM::_QLH  ]);
        HDF5.Write("SeaIceHeatFlux",   *fluxes[ THCM::_QTOS ]);
    }

    if (saveMask_)
    {
        INFO("Writing distributed and global mask to "  << filename);
        HDF5.Write("MaskLocal", *landmask_.local);

        HDF5.Write("MaskGlobal", "Global", H5T_NATIVE_INT,
                   landmask_.global->size(), &(*landmask_.global)[0]);

        HDF5.Write("MaskGlobal", "GlobalSize", (int) landmask_.global->size());

        HDF5.Write("MaskGlobal", "Surface", H5T_NATIVE_INT,
                   landmask_.global_surface->size(), &(*landmask_.global_surface)[0]);

        HDF5.Write("MaskGlobal", "Label", landmask_.label);
    }
    TIMER_STOP("Ocean: additionalExports");
}

// =====================================================================
void Ocean::additionalImports(EpetraExt::HDF5 &HDF5, std::string const &filename)
{
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

        // Import HDF5 data into THCM. This should not be
        // factorized as we cannot be sure what Map is going to come
        // out of the HDF5 Read call.

        // Create empty salflux vector
        Teuchos::RCP<Epetra_Vector> salflux =
            Teuchos::rcp(new Epetra_Vector(*domain_->GetStandardSurfaceMap()));

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

            delete readAdaptedSalFlux;

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

            delete readSalFluxPert;

            // Let THCM insert the salinity flux perturbation mask
            THCM::Instance().setEmip(salFluxPert, 'P');
        }

        delete readSalFlux;

        INFO("Loading salinity flux from " << filename << " done");
    }

    if (loadTemperatureFlux_)
    {
        INFO("Loading temperature flux from " << filename);
        if (!HDF5.IsContained("TemperatureFlux"))
        {
            ERROR("The group <SalinityFlux> is not contained in hdf5 " << filename,
                  __FILE__, __LINE__);
        }

        Epetra_MultiVector *readTemFlux;
        HDF5.Read("TemperatureFlux", readTemFlux);

        // This should not be factorized as we cannot be sure what Map
        // is going to come out of the HDF5.Read call.

        // Create empty temflux vector
        Teuchos::RCP<Epetra_Vector> temflux =
            Teuchos::rcp(new Epetra_Vector(*domain_->GetStandardSurfaceMap()));

        Teuchos::RCP<Epetra_Import> lin2solve_surf =
            Teuchos::rcp(new Epetra_Import( temflux->Map(),
                                            readTemFlux->Map() ));

        temflux->Import(*((*readTemFlux)(0)), *lin2solve_surf, Insert);

        // Instruct THCM to set/insert this as tatm in the local model
        THCM::Instance().setTatm(temflux);

        delete readTemFlux;

        INFO("Loading temperature flux from " << filename << " done");
    }

    if (loadMask_)
    {
        INFO("Loading local mask from " << filename);
        if (!HDF5.IsContained("MaskLocal") ||
            !HDF5.IsContained("MaskGlobal")
            )
        {
            WARNING("The group <Mask*> is not contained in hdf5, continue with standard mask...\n  "
                    << filename, __FILE__, __LINE__);
        }
        else
        {
            //__________________________________________________
            // We begin with the local (distributed) landmask
            Epetra_IntVector *readMask;

            // Obtain current mask to get distributed map with current
            // domain decomposition.
            Teuchos::RCP<Epetra_IntVector> tmpMask =
                THCM::Instance().getLandMask("current");

            // Read mask in hdf5 with distributed map
            HDF5.Read("MaskLocal", readMask);

            Teuchos::RCP<Epetra_Import> lin2dstr =
                Teuchos::rcp(new Epetra_Import( tmpMask->Map(),
                                                readMask->Map() ));

            tmpMask->Import(*readMask, *lin2dstr, Insert);

            delete readMask;

            // Put the new mask in THCM
            THCM::Instance().setLandMask(tmpMask, true);

            //__________________________________________________
            // Get global mask
            int globMaskSize;
            INFO("Loading global mask from " << filename);
            HDF5.Read("MaskGlobal", "GlobalSize", globMaskSize);

            std::shared_ptr<std::vector<int> > globmask =
                std::make_shared<std::vector<int> >(globMaskSize, 0);

            HDF5.Read("MaskGlobal", "Global", H5T_NATIVE_INT,
                      globMaskSize, &(*globmask)[0]);

            // Put the new global mask in THCM
            THCM::Instance().setLandMask(globmask);
        }
    }
}

//====================================================================
// Perform projection to vector: v = v - (v, s1) s1 - (v, s2) s2,
// where s1 and s2 are pressure checkerboard modes.
void Ocean::pressureProjection(Teuchos::RCP<Epetra_Vector> vec)
{
    assert(vec->Map().SameAs(state_->Map()));

    int grow = 0, lrow = 0, id3;
    double sum1 = 0.0, sum2 = 0.0;
    double gsum1 = 0.0, gsum2 = 0.0;

    Teuchos::RCP<Epetra_Vector> s1 = Teuchos::rcp(new Epetra_Vector(*vec));
    Teuchos::RCP<Epetra_Vector> s2 = Teuchos::rcp(new Epetra_Vector(*vec));
    s1->PutScalar(0.0);
    s2->PutScalar(0.0);

    for (int k = 0; k != L_; ++k)
        for (int j = 0; j != M_; ++j)
            for (int i = 0; i != N_; ++i)
            {
                // 3D mask index
                id3 = FIND_ROW2(1, N_, M_, L_, i, j, k, 1);
                if ((*landmask_.global_borderless)[id3] == 0)
                {
                    // obtain pressure row
                    grow = FIND_ROW2(_NUN_, N_, M_, L_, i, j, k, PP);
                    lrow = vec->Map().LID(grow);
                    if (lrow >= 0)
                    {
                        // create pressure modes
                        if ( (i + j) % 2 == 0 )
                        {
                            (*s1)[lrow] = 1;
                            sum1 = sum1 + 1.0;
                        }
                        else
                        {
                            (*s2)[lrow] = 1;
                            sum2 = sum2 + 1.0;
                        }
                    }
                }
            }

    comm_->SumAll(&sum1, &gsum1, 1);
    comm_->SumAll(&sum2, &gsum2, 1);

    s1->Scale(1.0 / sqrt(sum1));
    s2->Scale(1.0 / sqrt(sum2));

    double dp1 = Utils::dot(vec, s1);
    double dp2 = Utils::dot(vec, s2);

    INFO("Ocean: pressure checkerboard correction...");
    // v = v - (v, s1) s1 - (v, s2) s2
    vec->Update(-dp1, *s1, -dp2, *s2, 1.0);
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
int Ocean::npar() { return _NPAR_; }

//===================================================================
std::string Ocean::int2par(int ind) const
{
    return THCM::Instance().int2par(ind+1);
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

Teuchos::ParameterList
Ocean::getDefaultInitParameters()
{
    Teuchos::ParameterList result = getDefaultParameters();
    result.setName("Default Init Ocean List");

    result.get("Input file",  "ocean_input.h5");
    result.get("Output file", "ocean_output.h5");
    result.get("Save mask", true);
    result.get("Load mask", true);

    result.get("Load state", false);
    result.get("Save state", true);
    result.get("Save frequency", 0);

    result.get("Load salinity flux", false);
    result.get("Save salinity flux", true);
    result.get("Load temperature flux", false);
    result.get("Save temperature flux", true);

    result.get("Use legacy fort.3 output", false);
    result.get("Use legacy fort.44 output", true);
    result.get("Save column integral", false);
    result.get("Max mask fixes", 5);

    result.get("Analyze Jacobian", true);

    Teuchos::ParameterList& solverParams = result.sublist("Belos Solver");
    solverParams.get("FGMRES iterations", 500);
    solverParams.get("FGMRES tolerance", 1e-8);
    solverParams.get("FGMRES restarts", 0);
    solverParams.get("FGMRES output", 100);
    solverParams.get("FGMRES explicit residual test", false);

    result.sublist("THCM") = THCM::getDefaultInitParameters();

    return result;
}

Teuchos::ParameterList
Ocean::getDefaultParameters()
{
    Teuchos::ParameterList result("Default Ocean List");
    result.sublist("THCM") = THCM::getDefaultParameters();

    return result;
}

const Teuchos::ParameterList& Ocean::getParameters()
{
  params_.sublist("THCM") = thcm_->getParameters();
  return params_;
}

void Ocean::setParameters(Teuchos::ParameterList& newParams)
{
    newParams.validateParameters(getDefaultParameters());
    thcm_->setParameters(newParams.sublist("THCM"));
    params_.setParameters(newParams);
}
