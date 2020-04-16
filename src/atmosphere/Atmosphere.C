#include "Atmosphere.H"
#include "AtmosphereDefinitions.H"
#include "Ocean.H"
#include "SeaIce.H"

#include "Epetra_Import.h"

// Import/export
#include <EpetraExt_Exception.h>
#include <EpetraExt_HDF5.h>


//==================================================================
// Constructor
Atmosphere::Atmosphere(Teuchos::RCP<Epetra_Comm> comm, ParameterList params)
    :
    params_          (params),
    n_               (params->get("Global Grid-Size n", 16)),
    m_               (params->get("Global Grid-Size m", 16)),
    l_               (params->get("Global Grid-Size l", 1)),
    periodic_        (params->get("Periodic", false)),
    aux_             (params->get("Auxiliary unknowns", 1)),
    useIntCondQ_     (params->get("Use integral condition on q", true)),
    useFixedPrecip_  (params->get("Use idealized precipitation", false)),

    precInitialized_ (false),
    recomputePrec_   (false),
    recompMassMat_   (true)
{
    INFO("Atmosphere: constructor...");

    // inherited input/output datamembers
    inputFile_  = params->get("Input file", "atmos_input.h5");
    outputFile_ = params->get("Output file", "atmos_output.h5");
    loadState_  = params->get("Load state", false);
    saveState_  = params->get("Save state", true);
    saveMask_   = params->get("Save mask", true);
    saveEvery_  = params->get("Save frequency", 0);

    // initialize postprocessing counter
    ppCtr_ = 0;

    // set comm
    comm_ = comm;

    // Define degrees of freedom
    dof_ = ATMOS_NUN_;
    dim_ = n_ * m_ * l_ * dof_ + aux_;

    // Set integral condition row
    if (useIntCondQ_)
        rowIntCon_ = FIND_ROW_ATMOS0(dof_, n_, m_, l_,
                                     n_-1, m_-1, l_-1, ATMOS_QQ_);
    else
        rowIntCon_ = -1;

    // Define domain
    xmin_ = params->get("Global Bound xmin", 286.0) * PI_ / 180.0;
    xmax_ = params->get("Global Bound xmax", 350.0) * PI_ / 180.0;
    ymin_ = params->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
    ymax_ = params->get("Global Bound ymax", 74.0)  * PI_ / 180.0;

    INFO("Parallel atmosphere model: xmin = " << xmin_);
    INFO("                           xmax = " << xmax_);
    INFO("                           ymin = " << ymin_);
    INFO("                           ymax = " << ymax_);

    // Create domain object
    domain_ = Teuchos::rcp(new TRIOS::Domain(n_, m_, l_, dof_,
                                             xmin_, xmax_, ymin_, ymax_,
                                             periodic_, 1.0, 1.0, comm_, aux_));

    // Compute 2D decomposition
    domain_->Decomp2D();

    // Obtain local dimensions
    double xminloc = domain_->XminLoc();
    double xmaxloc = domain_->XmaxLoc();
    double yminloc = domain_->YminLoc();
    double ymaxloc = domain_->YmaxLoc();

    // Obtain local grid dimensions
    int nloc = domain_->LocalN();
    int mloc = domain_->LocalM();
    int lloc = domain_->LocalL();

    INFO("   Local atmosphere model: xmin = " << xminloc);
    INFO("                           xmax = " << xmaxloc);
    INFO("                           ymin = " << yminloc);
    INFO("                           ymax = " << ymaxloc);

    // Obtain overlapping and non-overlapping maps
    assemblyMap_ = domain_->GetAssemblyMap(); // overlapping
    standardMap_ = domain_->GetStandardMap(); // non-overlapping

    // Obtain special maps
    // depth-averaged, single unknown for ocean surface temperature
    standardSurfaceMap_ = domain_->CreateStandardMap(1, true);
    assemblySurfaceMap_ = domain_->CreateAssemblyMap(1, true);

    // Create Import object for single unknown surface values
    as2std_surf_ =
        Teuchos::rcp(new Epetra_Import(*assemblySurfaceMap_, *standardSurfaceMap_));

    // Create overlapping and non-overlapping vectors
    state_      = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    rhs_        = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    diagB_      = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    sol_        = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    lst_        = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    sst_        = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    sit_        = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    Msi_        = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    E_          = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    P_          = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    Pdist_      = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

    localState_ = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localRHS_   = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localDiagB_ = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localSol_   = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localLST_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localSST_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localSIT_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localMSI_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localE_     = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localP_     = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));

    // create graph
    createMatrixGraph();

    jac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *matrixGraph_));

    // Periodicity is handled by AtmosLocal if there is a single
    // core in the x-direction.
    Teuchos::RCP<Epetra_Comm> xComm = domain_->GetProcRow(0);
    bool perio = (periodic_ && xComm->NumProc() == 1);

    // Create local AtmosLocal object
    atmos_ = std::make_shared<AtmosLocal>(nloc, mloc, lloc, perio,
                                          xminloc, xmaxloc, yminloc, ymaxloc,
                                          params_);

    surfmask_ = std::make_shared<std::vector<int> >(m_ * n_, 0);

    // Import existing state
    if (loadState_)
        loadStateFromFile(inputFile_);

    // Create restricted maps and their importers:
    // Target map: Maps_[i]
    // Source map: state_->Map()
    int XX = 0;
    for (int i = 0; i != dof_; ++i)
    {
        XX = ATMOS_TT_ + i; // unknown
        Maps_[XX] = Utils::CreateSubMap(*standardMap_, dof_, XX);
        Imps_[XX] = Teuchos::rcp(new Epetra_Import(*Maps_[XX], *standardMap_));
    }

    // Build diagonal mass matrix
    computeMassMat();

    setupIntCoeff();

    INFO("Atmosphere: constructor done");
}

//==================================================================
void Atmosphere::setupIntCoeff()
{
    //------------------------------------------------------------------
    // Create parallelized integration coef. for integral condition q
    //------------------------------------------------------------------

    if (useIntCondQ_)
    {
        INFO("Atmosphere constructor: integral condition in row "
             << rowIntCon_);
    }
    else
    {
        WARNING("AtmposherePar: integral condition on q disabled!",
                __FILE__, __LINE__);
    }

    intcondCoeff_ = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    Teuchos::RCP<Epetra_Vector> intcondLocal =
        Teuchos::rcp(new Epetra_Vector(*assemblyMap_));

    std::vector<double> vals, inds;

    atmos_->integralCoeff(vals, inds);

    for (size_t idx = 0; idx != inds.size(); ++idx)
        (*intcondLocal)[inds[idx]-1] = vals[idx];

    // Assemble distributed version into non-overlapping vector
    domain_->Assembly2Solve(*intcondLocal, *intcondCoeff_);

    // Create allgathered version
    intcondGlob_ = Utils::AllGather(*intcondCoeff_);

// #ifdef DEBUGGING_NEW
//     std::stringstream ss1, ss2;
//     ss1 << "intcondq" << comm_->MyPID() << ".txt";
//     ss2 << "intcondq" << comm_->MyPID() << "orig.txt";
//     Utils::print(intcondCoeff_, ss1.str());
//     Utils::print(vals, ss2.str());
// #endif

    //------------------------------------------------------------------
    // Create parallelized integration coefficients for precipitation
    //------------------------------------------------------------------

    pIntCoeff_ = Teuchos::rcp( new Epetra_Vector(*standardSurfaceMap_) );

    Teuchos::RCP<Epetra_Vector> pIntCoeffLocal =
        Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));

    // Obtain integration coefficients for precipitation integral
    // Use 1 dof and ignore land
    atmos_->integralCoeff(vals, inds, 1);

    // Final index should equal local length in assemblysurfacemap
    assert(inds.back()-1 < pIntCoeffLocal->MyLength());

    // Fill local precipitation integration coefficients
    for (size_t idx = 0; idx != inds.size(); ++idx)
        (*pIntCoeffLocal)[inds[idx]-1] = vals[idx];

    // Export assembly map surface integration coeffs to standard map
    CHECK_ZERO( pIntCoeff_->Export( *pIntCoeffLocal, *as2std_surf_, Zero ) );

    // Obtain total integration area (sum of absolute values)
    pIntCoeff_->Norm1(&totalArea_);

    INFO("Atmosphere: total E,P area = " << totalArea_);
    INFO("Atmosphere:    local dA[0] = " << (*pIntCoeff_)[0]);

    // Set precipitation distribution
    setPdist();

}

//==================================================================
// --> If this turns out costly we might need to optimize using stateHash
void Atmosphere::distributeState()
{
    TIMER_START("Atmosphere: distribute state...");
    // Create assembly state
    domain_->Solve2Assembly(*state_, *localState_);

    // local problem size
    int numMyElements = assemblyMap_->NumMyElements();

    std::shared_ptr<std::vector<double> > localState =
        std::make_shared<std::vector<double> >(numMyElements, 0.0);

    localState_->ExtractCopy(&(*localState)[0], numMyElements);

    atmos_->setState(localState);

    TIMER_STOP("Atmosphere: distribute state...");
}

//==================================================================
void Atmosphere::computeRHS()
{
    TIMER_START("Atmosphere: computeRHS...");

    //------------------------------------------------------------------
    // Put parallel state in serial atmosphere.
    distributeState();

    if (aux_ <= 0)    // Compute E and P, put our precipitation field
                      // in serial AtmosLocal
    {
        // compute and obtain evaporation field
        setE();

        // compute and obtain precipitation field
        setP();

        // return precipitation field to local atmosphere
        CHECK_ZERO(localP_->Import(*P_, *as2std_surf_, Insert));

        int numMySurfaceElements = assemblySurfaceMap_->NumMyElements();
        std::shared_ptr<std::vector<double> > localP =
            std::make_shared<std::vector<double> >(numMySurfaceElements, 0.0);

        localP_->ExtractCopy(&(*localP)[0], numMySurfaceElements);

        atmos_->setPrecipitation(localP);
    }

    //------------------------------------------------------------------
    // compute local rhs and check bounds
    //------------------------------------------------------------------
    atmos_->computeRHS();
    std::shared_ptr<std::vector<double> > localRHS = atmos_->getRHS('V');

    int numMyElements = assemblyMap_->NumMyElements();
    if ((int) localRHS->size() != numMyElements)
    {
        std::cout << numMyElements    << std::endl;
        std::cout << localRHS->size() << std::endl;
        ERROR("RHS incorrect size", __FILE__, __LINE__);
    }

    // obtain view of overlapping rhs points
    double *rhs_tmp;
    localRHS_->ExtractView(&rhs_tmp);

    // fill view
    for (int i = 0; i != numMyElements; ++i)
        rhs_tmp[i] = (*localRHS)[i];

    // assemble distributed rhs into global rhs
    domain_->Assembly2Solve(*localRHS_, *rhs_);

    //------------------------------------------------------------------
    // set integral condition in RHS
    //------------------------------------------------------------------

    // If we have row rowIntCon_
    int root = comm_->NumProc()-1;

    if ( (rhs_->Map().MyGID(rowIntCon_)) && (comm_->MyPID() != root) )
    {
        ERROR("Integral should be on last processor!", __FILE__, __LINE__);
    }

    double intcond = Utils::dot(intcondCoeff_, state_);

    if (rhs_->Map().MyGID(rowIntCon_) && useIntCondQ_)
        (*rhs_)[rhs_->Map().LID(rowIntCon_)] = intcond;

    //------------------------------------------------------------------
    // Specify precipitation integral in RHS if aux > 0
    //
    // Precipitation is given by the integral of the evaporation
    // anomaly E = sigma - q + Msi*Cs, scaled with total area A
    //------------------------------------------------------------------

    // Obtain parameters
    Atmosphere::CommPars pars;
    getCommPars(pars);

    // Here we integrate sigma = dqso * sst + Msi*(dqsi*sit-dqso*sst)
    // using the increments in pIntCoeff.
    assert(sst_->Map().SameAs(sit_->Map()));

    Teuchos::RCP<Epetra_Vector> tmp =
        Teuchos::rcp(new Epetra_Vector(*sit_));

    // tmp = dqsi*sit-dqso*sst
    tmp->Update(-pars.dqso, *sst_, pars.dqsi);

    Teuchos::RCP<Epetra_Vector> sigma =
        Teuchos::rcp(new Epetra_Vector(*sst_));

    // Pointwise product and update:
    // sigma = dqso * sst + Msi*(dqsi*sit-dqso*sst)
    sigma->Multiply(1.0, *Msi_, *tmp, pars.dqso);

    // sigma is integrated and adjusted with (1/A)*(tdim/qdim)
    double sstInt = Utils::dot(pIntCoeff_, sigma) *
        (1.0 / totalArea_) * ( pars.tdim / pars.qdim ) ;

    // The humidity integral is the same as computed in the integral
    // condition above, so this can probably be simplified. We can
    // substitute it with 0 but for now we leave it and test that
    // later.
    double qInt = intcond * 1.0 / totalArea_;

    // Sublimation correction integrated over the sea ice surface int(M*Cs)
    double MCsInt = Utils::dot(pIntCoeff_, Msi_) * pars.Cs / totalArea_;

    // The integrals and P are on the same processor
    int last = FIND_ROW_ATMOS0( ATMOS_NUN_, n_, m_, l_, n_-1, m_-1, l_-1, ATMOS_NUN_ );
    int lid  = -1;

    if ( rhs_->Map().MyGID(rowIntCon_) && (aux_ == 1) )
    {
        lid = rhs_->Map().LID(last + 1);

        // F = - P - \sum_i (1 / A) * (q_i * dA_i)
        //         + \sum_i (1 / A) * (tdim / qdim) * sigma * dA_i
        //         + \sum_i (1 / A) * M * Cs dA_i
        //
        (*rhs_)[lid] = -(*state_)[lid] - qInt + sstInt + MCsInt;
    }
    TIMER_STOP("Atmosphere: computeRHS...");
}

//==================================================================
void Atmosphere::idealized(double precip)
{
    // initialize local rhs with idealized values
    atmos_->idealized(precip);

    // local problem size
    int numMyElements        = assemblyMap_->NumMyElements();
    int numMySurfaceElements = assemblySurfaceMap_->NumMyElements();

    // obtain view of assembly state
    double *state_tmp;
    localState_->ExtractView(&state_tmp);

    // obtain local state and check bounds
    std::shared_ptr<std::vector<double> > state = atmos_->getState('V');
    if ((int) state->size() != numMyElements)
    {
        ERROR("state incorrect size", __FILE__, __LINE__);
    }

   // fill assembly view with local state
    for (int i = 0; i != numMyElements; ++i)
    {
        state_tmp[i] = (*state)[i];
    }
    // set solvemap state
    domain_->Assembly2Solve(*localState_, *state_);

    INFO("Idealized state norm: " << Utils::norm(state_) );

    // obtain view of sst
    double *sst_tmp;
    localSST_->ExtractView(&sst_tmp);

    // obtain local sst and check bounds
    std::shared_ptr<std::vector<double> > sst = atmos_->getSST('V');
    if ((int) sst->size() != numMySurfaceElements)
    {
        ERROR("sst incorrect size", __FILE__, __LINE__);
    }

    for (int i = 0; i != numMySurfaceElements; ++i)
    {
        sst_tmp[i] = (*sst)[i];
    }

    // Export local values to distributed non-overlapping sst_
    CHECK_ZERO( sst_->Export( *localSST_, *as2std_surf_, Zero ) );

    INFO("Idealized   sst norm: " << Utils::norm(sst_) );
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::interface(int XX) const
{
    if ( XX >= ATMOS_PP_ )
    {
        if ( aux_ <= 0 )
        {
            WARNING("Auxiliary P disabled: "
                    << "P is 'diagnosed'", __FILE__, __LINE__);
        }

        return getP();
    }
    else
    {
        Teuchos::RCP<Epetra_Vector> out =
            Teuchos::rcp(new Epetra_Vector(*Maps_.at(XX)));

        CHECK_ZERO(out->Import(*state_, *Imps_.at(XX), Insert));
        return out;
    }
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::interfaceT() const
{
    return interface(ATMOS_TT_);
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::interfaceQ() const
{
    return interface(ATMOS_QQ_);
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::interfaceA() const
{
    return interface(ATMOS_AA_);
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::interfaceP() const
{
    return interface(ATMOS_PP_);
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::interfaceE()
{
    return getE();
}

//==================================================================
std::shared_ptr<Utils::CRSMat> Atmosphere::getBlock(std::shared_ptr<Ocean> ocean)
{
    TIMER_START("Atmosphere::getBlock(ocean)...");
    // Jacobian of the atmosphere with respect to the ocean model, see
    // the forcing in AtmosLocal.C.

    // check surfmask
    assert( (int) surfmask_->size() == m_*n_ );

    // We are going to create a 0-based global CRS matrix for this block
    std::shared_ptr<Utils::CRSMat> block = std::make_shared<Utils::CRSMat>();

    int el_ctr = 0;
    int oceanTT = 5; // (1-based) in THCM temperature is the fifth unknown

    // Obtain atmosphere commpars
    AtmosLocal::CommPars pars;
    getCommPars(pars);

    // FIXME this now depends on the sea ice mask: we enther the realm
    // of non-constant coupling coefficients. For now we use an
    // allgather to get the full sea ice mask, in the future we could
    // let this block be partly computed by the local model and
    // assemble here. Need to figure that out.
    Teuchos::RCP<Epetra_MultiVector> Msi = Utils::AllGather(*Msi_);

    int sr; // surface row

    double dTFT;  // d / dT_ocean (F_T)
    double dTFQ;  // d / dT_ocean (F_Q)
    double M;     // Mask value

    bool on_ocean;
    bool no_intrw;

    // loop over our unknowns
    for (int j = 0; j != m_; ++j)
        for (int i = 0; i != n_; ++i)
        {
            sr   = j*n_+i; // set surface row

            M = (*(*Msi)(0))[sr];

            dTFT = 1.0 - M;

            dTFQ = pars.nuq * pars.tdim / pars.qdim * pars.dqso * (1.0 - M);

            on_ocean = ((*surfmask_)[sr] == 0);
            for (int xx = ATMOS_TT_; xx <= dof_; ++xx)
            {
                no_intrw = (rowIntCon_ != FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_,
                                                          i, j, l_-1, xx));

                block->beg.push_back(el_ctr);

                if ( on_ocean && no_intrw ) // ocean points and skip integral cond
                {
                    switch (xx)
                    {

                    case ATMOS_TT_:
                        block->co.push_back(dTFT);
                        block->jco.push_back(ocean->interface_row(i,j,oceanTT));
                        el_ctr++;
                        break;

                    case ATMOS_QQ_:
                        block->co.push_back( dTFQ );
                        block->jco.push_back( ocean->interface_row(i,j,oceanTT) );
                        el_ctr++;
                        break;
                    }
                }
            }
        }

    // add dependencies of precipitation row
    int qid;
    double dTFP; // d / dTo (F_P)

    if (aux_ == 1)
    {
        block->beg.push_back(el_ctr);
        for (int j = 0; j != m_; ++j)
            for (int i = 0; i != n_; ++i)
            {
                sr = j*n_+i;                  // set surface row
                M  = (*(*Msi)(0))[sr];        // sea ice mask

                if ( (*surfmask_)[sr] == 0)   // non-land
                {
                    qid = FIND_ROW_ATMOS0( ATMOS_NUN_, n_, m_, l_,
                                           i, j, l_-1, ATMOS_QQ_ );

                    dTFP = ( *intcondGlob_ )[0][qid] * ( 1.0 / totalArea_ )
                        * ( pars.tdim / pars.qdim ) * pars.dqso * ( 1.0 - M );

                    block->co.push_back( dTFP );

                    block->jco.push_back( ocean->interface_row(i,j,oceanTT) );
                    el_ctr++;
                }
            }
    }

    block->beg.push_back(el_ctr);

    assert( (int) block->co.size() == block->beg.back());

    TIMER_STOP("Atmosphere::getBlock(ocean)...");
    return block;
}

//==================================================================
std::shared_ptr<Utils::CRSMat> Atmosphere::getBlock(std::shared_ptr<SeaIce> seaice)
{
    TIMER_START("Atmosphere::getBlock(seaice)...");
    // Jacobian of the atmosphere with respect to the sea ice model,
    // see AtmosLocal::forcing()

    // initialize empty CRS matrix
    std::shared_ptr<Utils::CRSMat> block = std::make_shared<Utils::CRSMat>();

    int el_ctr = 0;

    int seaiceMM = SEAICE_MM_; // (1-based) mask unknown in the sea ice model
    int seaiceTT = SEAICE_TT_; // (1-based) temperature unknown in the sea ice model

    // Obtain atmosphere parameters
    AtmosLocal::CommPars pars; getCommPars(pars);

    double dMFT;   // d / dMsi (F_T)
    double dMFQ;   // d / dMsi (F_Q)
    double dTFT;   // d / dTsi (F_T)
    double dTFQ;   // d / dTsi (F_Q)
    double dMFA;   // d / dMsi (F_A)

    int sr;     // surface row
    double M;   // mask value
    double To;  // sst value
    double Ti;  // sit value
    double Eo;  // evaporation value
    double Ei;  // sublimation value

    double Cs = pars.Cs;  // sublimation correction

    // FIXME: We use an allgather to get the full sea ice mask, sst
    // and sit. In the future we could let this block be partly
    // computed by the local model and assemble here. Need to figure
    // that out.
    Teuchos::RCP<Epetra_MultiVector> Msi = Utils::AllGather(*Msi_);
    Teuchos::RCP<Epetra_MultiVector> sst = Utils::AllGather(*sst_);
    Teuchos::RCP<Epetra_MultiVector> sit = Utils::AllGather(*sit_);


    for (int j = 0; j != m_; ++j)
        for (int i = 0; i != n_; ++i)
        {
            sr = j*n_+i;

            M  = (*(*Msi)(0))[sr];
            To = (*(*sst)(0))[sr];
            Ti = (*(*sit)(0))[sr];

            Eo = pars.tdim / pars.qdim * pars.dqso * To;
            Ei = pars.tdim / pars.qdim * pars.dqsi * Ti;

            dMFT = Ti + pars.t0i - To - pars.t0o;
            dTFT = M;

            dMFQ = pars.nuq  * (Ei - Eo + Cs);
            dTFQ = pars.nuq  * pars.tdim / pars.qdim * pars.dqsi * M;
            dMFA = pars.comb * pars.albf / pars.tauc;

            for (int xx = ATMOS_TT_; xx <= dof_; ++xx)
            {
                block->beg.push_back(el_ctr);

                // skip land and integral condition
                if ( (*surfmask_)[sr] == 0 &&
                     (rowIntCon_ != FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_,
                                                    i, j, l_-1, xx)))
                {
                    switch (xx)
                    {

                    case ATMOS_TT_:
                        block->co.push_back(dMFT);
                        block->jco.push_back(seaice->interface_row(i,j,seaiceMM));
                        el_ctr++;

                        block->co.push_back(dTFT);
                        block->jco.push_back(seaice->interface_row(i,j,seaiceTT));
                        el_ctr++;
                        break;

                    case ATMOS_QQ_:
                        block->co.push_back(dMFQ);
                        block->jco.push_back(seaice->interface_row(i,j,seaiceMM));
                        el_ctr++;

                        block->co.push_back(dTFQ);
                        block->jco.push_back(seaice->interface_row(i,j,seaiceTT));
                        el_ctr++;
                        break;

                    case ATMOS_AA_:
                        block->co.push_back(dMFA);
                        block->jco.push_back(seaice->interface_row(i,j,seaiceMM));
                        el_ctr++;
                        break;
                    }
                }
            }
        }

    // add dependencies of precipitation row
    int qid;
    double dMFP; // d / dMsi (F_P)
    double dTFP; // d / dTsi (F_P)
    double dA;   // integral coefficient

    if (aux_ == 1)
    {
        block->beg.push_back(el_ctr);
        for (int j = 0; j != m_; ++j)
            for (int i = 0; i != n_; ++i)
            {
                sr = j*n_+i; // set surface row

                M  = (*(*Msi)(0))[sr];
                To = (*(*sst)(0))[sr];
                Ti = (*(*sit)(0))[sr];


                if ( (*surfmask_)[sr] == 0)   // non-land
                {
                    qid = FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_,
                                          i, j, l_-1, ATMOS_QQ_);

                    dA   = (*intcondGlob_)[0][qid];

                    dMFP = (dA / totalArea_)
                        * ( (pars.tdim / pars.qdim) *
                            (pars.dqsi * Ti - pars.dqso * To)
                            + Cs );

                    block->co.push_back(dMFP);
                    block->jco.push_back(seaice->interface_row(i,j,seaiceMM));
                    el_ctr++;

                    dTFP = (dA / totalArea_)
                        * ( pars.tdim / pars.qdim ) * pars.dqsi * M;

                    block->co.push_back(dTFP);
                    block->jco.push_back(seaice->interface_row(i,j,seaiceTT));
                    el_ctr++;
                }
            }
    }

    block->beg.push_back(el_ctr);
    assert( (int) block->co.size() == block->beg.back());

    TIMER_STOP("Atmosphere::getBlock(seaice)...");
    return block;
}

//==================================================================
void Atmosphere::synchronize(std::shared_ptr<Ocean> ocean)
{
    // This is a simple interface. The atmosphere only needs the
    // ocean temperature at the ocean-atmosphere interface (SST).

    // Get ocean surface temperature
    Teuchos::RCP<Epetra_Vector> sst = ocean->interfaceT();

    // Set ocean surface temperature in parallel and serial atmosphere model.
    setOceanTemperature(sst);
}

//==================================================================
void Atmosphere::synchronize(std::shared_ptr<SeaIce> seaice)
{
    // Get sea ice mask
    Teuchos::RCP<Epetra_Vector> Msi = seaice->interfaceM();
    setSeaIceMask(Msi);

    // Get sea ice temperature
    Teuchos::RCP<Epetra_Vector> sit = seaice->interfaceT();
    setSeaIceTemperature(sit);
}

//==================================================================
void Atmosphere::setOceanTemperature(Teuchos::RCP<Epetra_Vector> sst)
{
    CHECK_MAP(sst, standardSurfaceMap_);

    // assign to our own datamember
    sst_ = sst;

    // create assembly
    // domain_->Solve2Assembly(*sst_, *localSST_);
    CHECK_ZERO(localSST_->Import(*sst_, *as2std_surf_, Insert));

    // local vector size
    int numMyElements = assemblySurfaceMap_->NumMyElements();

    std::vector<double> localSST(numMyElements, 0.0);

    localSST_->ExtractCopy(&localSST[0], numMyElements);
    atmos_->setOceanTemperature(localSST);
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::getLandTemperature()
{
    double *local;
    localLST_->ExtractView(&local);
    int numMyElements = localLST_->Map().NumMyElements();
    std::shared_ptr<std::vector<double> > lst = atmos_->getLandTemperature();

    assert((int) lst->size() == numMyElements);
    for (int i = 0; i != numMyElements; ++i)
        local[i] = (*lst)[i];

    CHECK_ZERO(lst_->Export(*localLST_, *as2std_surf_, Zero));
    return lst_;
}

//==================================================================
void Atmosphere::setSeaIceTemperature(Teuchos::RCP<Epetra_Vector> sit)
{
    // Replace map if necessary
    if (!(sit->Map().SameAs(*standardSurfaceMap_)))
    {
        // INFO("Atmosphere::setOceanTemperature sit map -> standardSurfaceMap_");
        CHECK_ZERO(sit->ReplaceMap(*standardSurfaceMap_));
    }

    sit_ = sit;
    CHECK_ZERO(localSIT_->Import(*sit_, *as2std_surf_, Insert));

    // local vector size
    int numMyElements = assemblySurfaceMap_->NumMyElements();

    std::vector<double> localSIT(numMyElements, 0.0);

    localSIT_->ExtractCopy(&localSIT[0], numMyElements);
    atmos_->setSeaIceTemperature(localSIT);
}

//==================================================================
void Atmosphere::setSeaIceMask(Teuchos::RCP<Epetra_Vector> mask)
{
    if (!(mask->Map().SameAs(*standardSurfaceMap_)))
    {
        CHECK_ZERO(mask->ReplaceMap(*standardSurfaceMap_));
    }

    Msi_ = mask;
    CHECK_ZERO(localMSI_->Import(*Msi_, *as2std_surf_, Insert));
    int numMyElements = assemblySurfaceMap_->NumMyElements();

    std::vector<double> localMSI(numMyElements, 0.0);
    localMSI_->ExtractCopy( &localMSI[0], numMyElements);

    atmos_->setSeaIceMask(localMSI);
}

//==================================================================
void Atmosphere::setLandMask(Utils::MaskStruct const &mask)
{
    // create global surface mask
    // we do the same thing for the local mask in the AtmosLocal object
    surfmask_->clear();

    if ((int) mask.global_surface->size() < (m_ * n_))
    {
        ERROR("mask.global_surface->size() not ok:",  __FILE__, __LINE__);
    }
    else // we trust surfm
    {
        surfmask_ = mask.global_surface;
    }

    // create rcp
    int numMyElements = mask.local->MyLength();

    std::shared_ptr<std::vector<int> > landmask =
        std::make_shared<std::vector<int> >(numMyElements, 0);

    CHECK_ZERO(mask.local->ExtractCopy(&(*landmask)[0]));

    // local atmosphere builds its own landmask from full distributed
    // mask
    atmos_->setSurfaceMask(landmask);

    // Some of the integral coefficients depend on the mask
    // so we repeat that setup.
    setupIntCoeff();
}

//==================================================================
void Atmosphere::computeJacobian()
{
    TIMER_START("Atmosphere: compute Jacobian...");
    // set all entries to zero
    CHECK_ZERO(jac_->PutScalar(0.0));

    // compute jacobian in local atmosphere
    atmos_->computeJacobian();

    // obtain 1-based CRS matrix from local atmosphere
    std::shared_ptr<Utils::CRSMat> localJac =
        atmos_->getJacobian();

    // max nonzeros per row
    const int maxnnz = ATMOS_NUN_ * ATMOS_NP_ + 1 + aux_;

    // indices array
    int indices[maxnnz];

    // values array
    double values[maxnnz];

    // check size
    int numMyElements = assemblyMap_->NumMyElements() - aux_;

    assert(numMyElements == (int) localJac->beg.size() - 1);

    // loop over local elements
    int index, numentries;
    for (int i = 0; i < numMyElements; ++i)
    {
        // ignore ghost rows
        if (!domain_->IsGhost(i, ATMOS_NUN_))
        {
            // obtain indices and values from CRS container
            index = localJac->beg[i]; // beg contains 1-based indices!
            numentries = localJac->beg[i+1] - index;
            for (int j = 0; j < numentries; ++j)
            {
                indices[j] = assemblyMap_->GID(localJac->jco[index-1+j] - 1);
                values[j]  = localJac->co[index-1+j];
            }

            // put values in Jacobian
            int ierr = jac_->ReplaceGlobalValues(assemblyMap_->GID(i),
                                                 numentries,
                                                 values, indices);

            // debugging
            if (ierr != 0)
            {
                int grow = assemblyMap_->GID(i);
                int lrid = jac_->Graph().LRID(grow);
                std::cout << " my PID: " << comm_->MyPID()
                          << ", ierr = " << ierr
                          << ", GRID = " << grow
                          << "  LRID = " << lrid << std::endl;

                std::cout << " graph inds in LRID:   "
                          << jac_->Graph().NumMyIndices(lrid) << std::endl;

                for (int ii = 0; ii < numentries; ++ii)
                {
                    std::cout << "proc " << comm_->MyPID()
                              << " entries: (" << indices[ii]
                              << " " << values[ii] << ")" <<  std::endl;

                    std::cout << "proc " << comm_->MyPID() << " "
                              << assemblyMap_->GID(i) << std::endl;

                    INFO(" debug info: " << indices[ii] << " " << values[ii]);
                }

                INFO(" GRID: "<< assemblyMap_->GID(i));
                INFO(" number of entries: " << numentries);
                INFO(" numMyElements: " << numMyElements);
                INFO(" is ghost " << domain_->IsGhost(i, ATMOS_NUN_));
                INFO(" maxnnz: " << maxnnz);

                CHECK_ZERO(jac_->ExtractGlobalRowCopy(assemblyMap_->GID(i),
                                                      maxnnz, numentries,
                                                      values, indices));
                INFO("\noriginal row: ");
                INFO("number of entries: "<<numentries);

                for (int ii = 0; ii < numentries; ++ii)
                {
                    std::cout << "proc" << comm_->MyPID()
                              << " entries: (" << indices[ii]
                              << " " << values[ii] << ")" <<  std::endl;
                    INFO(" debug info: " << indices[ii] << " " << values[ii]);
                }

                INFO ("Error in ReplaceGlobalValues: " << ierr);
                ERROR("Error in ReplaceGlobalValues", __FILE__, __LINE__);
            }
        }
    }

    //------------------------------------------------------------------
    // Implementation of integrals:
    //  - integral condition on Q
    //  - P integral in auxiliary row
    //------------------------------------------------------------------

    int root = comm_->NumProc()-1;

    // If we have row rowIntCon_
    if ( (jac_->MyGRID(rowIntCon_)) && (comm_->MyPID() != root) )
    {
        ERROR("Integral should be on last processor!", __FILE__, __LINE__);
    }

    // length is the same for all QQ integrals
    int len = n_ * m_ * l_ + aux_;
    int icinds[len];    // integral condition indices
    int ipinds[len];    // P integral indices
    double icvals[len]; // integral condition values
    double ipvals[len]; // P integral values

    // Obtain indices and values for integrals
    int pos = 0;
    int gid;

    // Obtain some constants from local model
    Atmosphere::CommPars pars;
    getCommPars(pars);

    for (int k = 0; k != l_; ++k)
        for (int j = 0; j != m_; ++j)
            for (int i = 0; i != n_; ++i)
            {
                gid = FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_, i, j, k, ATMOS_QQ_);
                icinds[pos] = gid;
                icvals[pos] = (*intcondGlob_)[0][gid];
                ipinds[pos] = gid;
                ipvals[pos] = (-1.0 / totalArea_ ) * (*intcondGlob_)[0][gid];
                pos++;
            }

    // Add auxiliary dependencies
    int last = FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_, n_-1, m_-1, l_-1, ATMOS_NUN_);
    for (int aa = 1; aa <= aux_; ++aa)
    {
        gid         = last + aa;
        icinds[pos] = gid;
        icvals[pos] = 0.0;

        // final dependency P -> P
        if (aa == 1)
        {
            ipinds[pos] = gid;
            ipvals[pos] = -1;
        }

        pos++;
    }

    // Integrity check
    assert(pos == len);

    int icerr = 0;
    int iperr = 0;
    int icerrGlob = 0;
    int iperrGlob = 0;

    // Set indices and values for integrals. Integral rows are on final
    // processor
    if (jac_->MyGRID(rowIntCon_))
    {
        if (jac_->Filled())
        {
            // Integral condition Q
            if (useIntCondQ_)
                icerr = jac_->ReplaceGlobalValues(rowIntCon_, len, icvals, icinds);

            if (aux_ > 0)
                // Global precipiation integral
                iperr = jac_->ReplaceGlobalValues(last + 1, len, ipvals, ipinds);
        }
        else
        {
            if (useIntCondQ_)
                icerr = jac_->InsertGlobalValues(rowIntCon_, len, icvals, icinds);

            if (aux_ > 0)
                iperr = jac_->InsertGlobalValues(last + 1, len, ipvals, ipinds);
        }
    }

    comm_->SumAll(&icerr, &icerrGlob, 1);
    comm_->SumAll(&iperr, &iperrGlob, 1);

    if ((icerrGlob != 0) || (iperrGlob != 0))
    {
        INFO( "Insertion ERROR! " << icerr << " " << iperr << " filled = "
              << jac_->Filled());
        INFO( " while inserting/replacing values in local Jacobian");
        INFO( "  GRID: " << rowIntCon_);
        ERROR("Error during insertion/replacing of values in local Jacobian",
              __FILE__, __LINE__);
    }

    //------------------------------------------------------------------
    // Finalize matrix
    CHECK_ZERO(jac_->FillComplete());

// #ifdef DEBUGGING_NEW
//     DUMPMATLAB("atmos_jac", *jac_);
// #endif

    // With a new Jacobian we need to recompute the factorization
    recomputePrec_ = true;

    TIMER_STOP("Atmosphere: compute Jacobian...");
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::getE(char mode)
{
    // compute dimensional E in local AtmosLocal
    atmos_->computeEvaporation();

    // obtain view of E from serial AtmosLocal
    std::shared_ptr<std::vector<double> > localE = atmos_->interfaceE('V');

    // assign obtained E values to distributed vector (overlapping)
    double *tmpE;
    localE_->ExtractView( &tmpE );

    int numMySurfaceElements = assemblySurfaceMap_->NumMyElements();

    // fill localE_ with dimensional evaporation
    for (int i = 0; i != numMySurfaceElements; ++i)
        tmpE[i] = (*localE)[i];

    // export overlapping into non-overlapping E values
    CHECK_ZERO(E_->Export(*localE_, *as2std_surf_, Zero));

    return Utils::getVector(mode, E_);
}

//==================================================================
void Atmosphere::setE()
{
    getE('V');
}

//==================================================================
Teuchos::RCP<Epetra_Vector> Atmosphere::getP(Teuchos::RCP<Epetra_Vector> P) const
{
    if (P == Teuchos::null)
        P = Utils::getVector('C', P_);

    // P is a single, globally computed unknown
    double Pvalue = 0.0;
    int numGlobalElements = P->Map().NumGlobalElements();
    int numMyElements     = P->Map().NumMyElements();
    assert((int) surfmask_->size() == numGlobalElements);

    Atmosphere::CommPars pars;
    getCommPars(pars);

    // Without auxiliary unknowns we have to do the integral manually
    if (aux_ <= 0)
    {
        if (useFixedPrecip_)
            Pvalue = 1.0e-6;
        else
        {
            // Compute integral: we assume E_ is up to date and fully
            // dimensional
            Pvalue = Utils::dot(pIntCoeff_, E_) / totalArea_;

            // nondimensionalize Pvalue
            Pvalue = (Pvalue - pars.Eo0)/(pars.eta*pars.qdim);
        }
    }
    else if (aux_ == 1) // when aux = 1, P is an anomaly (state component)
    {
        // get last ordinary index
        int last = FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_,
                                   n_-1, m_-1, l_-1, ATMOS_NUN_);

        // get P value
        int lid;
        double PvalueLoc;
        if (state_->Map().MyGID(last+1))
        {
            lid       = state_->Map().LID(last+1);
            PvalueLoc = (*state_)[lid];
        }

        // communicate Pvalue
        comm_->SumAll(&PvalueLoc, &Pvalue, 1.0);
    }
    else
    {
        ERROR("Invalid aux", __FILE__, __LINE__);
    }

    // Fill dimensional P field. This can be homogoneous or applied
    // with any function f satisfying int f(theta) dA = A
    int gid;
    for (int i = 0; i != numMyElements; ++i)
    {
        gid = P->Map().GID(i);
        if ((*surfmask_)[gid] == 0)
            (*P)[i] = (*Pdist_)[i] *
                (pars.Eo0 + pars.eta * pars.qdim * Pvalue);
        else
            (*P)[i] = 0.0;
    }
    return P;
}

//==================================================================
void Atmosphere::setP()
{
    getP(P_);
}

//==================================================================
void Atmosphere::setPdist()
{
    INFO("Atmosphere::setPdist... ");
    Teuchos::RCP<Epetra_Vector> localPdist =
        Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));

    Pdist_->PutScalar(0.0);
    // let local model fill Pdist with some function
    domain_->Standard2AssemblySurface(*Pdist_, *localPdist);
    double *tmpPdist;
    localPdist->ExtractView(&tmpPdist);
    atmos_->fillPdist(tmpPdist);

    // correct global vector with integral
    domain_->Assembly2StandardSurface(*localPdist, *Pdist_);
    double corr = 1 - Utils::dot(pIntCoeff_, Pdist_) / totalArea_;

    // create vector of ones at non-land points
    Epetra_Vector ones(*standardSurfaceMap_);
    for (int i = 0; i != ones.Map().NumMyElements(); ++i)
        if (std::abs((*pIntCoeff_)[i]) > 1e-7)
            ones[i] = 1.0;

    CHECK_ZERO(Pdist_->Update(corr, ones, 1.0));

    // return the corrected Pdist in the local model
    domain_->Standard2AssemblySurface(*Pdist_, *localPdist);
    localPdist->ExtractView(&tmpPdist);
    atmos_->setPdist(tmpPdist);

    INFO("Atmosphere::setPdist... done");
}

//==================================================================
void Atmosphere::applyMatrix(Epetra_MultiVector const &in,
                                Epetra_MultiVector &out)
{
    TIMER_START("Atmosphere: apply matrix...");
    jac_->Apply(in, out);
    TIMER_STOP("Atmosphere: apply matrix...");
}

//====================================================================
void Atmosphere::computeMassMat()
{
    if (recompMassMat_)
    {
        INFO("Atmosphere: build mass matrix...");
        // compute mass matrix
        atmos_->computeMassMat();

        int numMyElements = localDiagB_->Map().NumMyElements();

        // obtain local diagonal
        std::shared_ptr<std::vector<double> > localDiagB =
            atmos_->getMassMat('V');

        if ((int) localDiagB->size() != numMyElements)
        {
            ERROR("Local diagB incorrect size", __FILE__, __LINE__);
        }

        // obtain view of assembly diagB
        double *tmpview;
        localDiagB_->ExtractView(&tmpview);

        // fill view
        for (int i = 0; i != numMyElements; ++i)
        {
            tmpview[i] = (*localDiagB)[i];
        }

        domain_->Assembly2Solve(*localDiagB_, *diagB_);

        // Set zero for integral condition on QQ
        if ( useIntCondQ_ && diagB_->Map().MyGID(rowIntCon_) )
        {
            (*diagB_)[diagB_->Map().LID(rowIntCon_)] = 0.0;
        }

        // Set zero for integral equation for PP
        if (aux_ == 1)
        {
            int last =
                FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_, n_-1, m_-1, l_-1, ATMOS_NUN_);

            int lid;
            if (diagB_->Map().MyGID(last+1))
            {
                lid = state_->Map().LID(last+1);
                (*diagB_)[lid] = 0.0;
            }
        }
        INFO("Atmosphere: build mass matrix... done");
    }

    recompMassMat_ = false; // Disable subsequent recomputes
}

//==================================================================
void Atmosphere::applyMassMat(Epetra_MultiVector const &v,
                                 Epetra_MultiVector &out)
{
    TIMER_START("Atmosphere: apply mass matrix...");

    // Compute mass matrix
    computeMassMat();

    // element-wise multiplication (out = 0.0*out + 1.0*B*v)
    out.Multiply(1.0, *diagB_, v, 0.0);

    TIMER_STOP("Atmosphere: apply mass matrix...");
}

//==================================================================
void Atmosphere::buildPreconditioner()
{
    INFO("Atmosphere: initialize preconditioner...");

    Ifpack Factory;
    std::string precType = "Amesos"; // direct solve on subdomains with some overlap
    int overlapLevel = params_->get("Ifpack overlap level", 2);

    // Create preconditioner
    precPtr_ = Teuchos::rcp(Factory.Create(precType, jac_.get(), overlapLevel));
    precPtr_->Initialize();
    precPtr_->Compute();

    precInitialized_ = true;

    INFO("Atmosphere: initialize preconditioner... done");
}

//==================================================================
void Atmosphere::getCommPars(Atmosphere::CommPars &parStruct) const
{
    atmos_->getCommPars(parStruct);
}

//==================================================================
void Atmosphere::preProcess()
{
    recomputePrec_ = true;
    recompMassMat_ = true;
}

//==================================================================
void Atmosphere::postProcess()
{
    // increase postprocessing counter
    ppCtr_++;

    // save state -> hdf5
    if (saveState_)
        saveStateToFile(outputFile_); // Save to hdf5

    if ((saveEvery_ > 0) && (ppCtr_ % saveEvery_) == 0)
    {
        std::stringstream append;
        append << "." << ppCtr_;
        copyState(append.str());
    }

}

//==================================================================
std::string Atmosphere::writeData(bool describe) const
{
    std::ostringstream datastring;

    if (describe)
    {
        datastring << std::setw(_FIELDWIDTH_)
                   << "max(Ta)"
                   << std::setw(_FIELDWIDTH_)
                   << "min(Ta)"
                   << std::setw(_FIELDWIDTH_)
                   << "max(Qa)"
                   << std::setw(_FIELDWIDTH_)
                   << "min(Qa)"
                   << std::setw(_FIELDWIDTH_)
                   << "||A||";

        return datastring.str();
    }
    else
    {
        datastring.precision(_PRECISION_);

        double maxT, minT;
        interfaceT()->MinValue(&minT);
        interfaceT()->MaxValue(&maxT);

        double maxQ, minQ;
        interfaceQ()->MinValue(&minQ);
        interfaceQ()->MaxValue(&maxQ);

        double normA = Utils::norm(interfaceA());

        datastring << std::scientific << std::setw(_FIELDWIDTH_)
                   << maxT << std::setw(_FIELDWIDTH_)
                   << minT << std::setw(_FIELDWIDTH_)
                   << maxQ << std::setw(_FIELDWIDTH_)
                   << minQ << std::setw(_FIELDWIDTH_)
                   << normA ;

        return datastring.str();
    }
}

//==================================================================
void Atmosphere::applyPrecon(Epetra_MultiVector const &in,
                                Epetra_MultiVector &out)
{
    TIMER_START("Atmosphere: apply preconditioner...");
    if (!precInitialized_)
    {
        buildPreconditioner();
    }
    if (recomputePrec_)
    {
        INFO("Atmosphere: recomputing prec");
        // precPtr_->Initialize();
        precPtr_->Compute();
        recomputePrec_ = false;
    }
    precPtr_->ApplyInverse(in, out);

    // check matrix residual
    // Teuchos::RCP<Epetra_MultiVector> r =
    //     Teuchos::rcp(new Epetra_MultiVector(in));;

    // applyMatrix(out, *r);
    // r->Update(1.0, in, -1.0);
    // double rnorm = Utils::norm(r);

    // INFO("Atmosphere: preconditioner residual: " << rnorm);


    TIMER_STOP("Atmosphere: apply preconditioner...");
}

//==================================================================
void Atmosphere::solve(Teuchos::RCP<Epetra_MultiVector> const &b)
{
    // when using the preconditioner as a solver make sure the overlap
    // is large enough (depending on number of cores obv).
    applyPrecon(*b, *sol_);
}

//==================================================================
// Very similar to the THCM function
// Create a graph to initialize the Jacobian
void Atmosphere::createMatrixGraph()
{
    // We know from the discretization that we have at most 5 + aux
    // dependencies in each row. This will change, obviously, when
    // extending the atmosphere model. If the number of dependencies
    // varies greatly per row we need to supply them differently
    // (variable).
    int maxDeps = 5 + aux_;
    matrixGraph_ =
        Teuchos::rcp(new Epetra_CrsGraph(Copy, *standardMap_, maxDeps, false));

    // Here we start specifying the indices
    int indices[maxDeps];

    // Get global domain size
    int N = domain_->GlobalN();
    int M = domain_->GlobalM();
    int L = domain_->GlobalL();

    // Get our local range in all directions
    // 0-based
    int I0 = domain_->FirstRealI();
    int J0 = domain_->FirstRealJ();
    int K0 = domain_->FirstRealK();
    int I1 = domain_->LastRealI();
    int J1 = domain_->LastRealJ();
    int K1 = domain_->LastRealK();

    int pos; // position in indices array
    int gidT, gid0;

    // Last ordinary row in the grid, i.e., the last non-auxiliary row.
    int last = FIND_ROW_ATMOS0(ATMOS_NUN_, N, M, L, N-1, M-1, L-1, ATMOS_NUN_);

    for (int k = K0; k <= K1; ++k)
        for (int j = J0; j <= J1; ++j)
            for (int i = I0; i <= I1; ++i)
            {
                // Obtain offset row corresponding to i,j,k,TT
                gidT = FIND_ROW_ATMOS0(ATMOS_NUN_, N, M, L, i, j, k, ATMOS_TT_);
                gid0 = gidT - 1; // used as offset
                pos  = 0;

                // T-equation

                // Specify dependencies, see Atmosphere::discretize()
                // ATMOS_TT_-ATMOS_TT_: 5-point stencil
                insert_graph_entry(indices, pos, i, j, k,   ATMOS_TT_, N, M, L);
                insert_graph_entry(indices, pos, i-1, j, k, ATMOS_TT_, N, M, L);
                insert_graph_entry(indices, pos, i+1, j, k, ATMOS_TT_, N, M, L);
                insert_graph_entry(indices, pos, i, j-1, k, ATMOS_TT_, N, M, L);
                insert_graph_entry(indices, pos, i, j+1, k, ATMOS_TT_, N, M, L);

                // ATMOS_TT_-ATMOS_AA_
                insert_graph_entry(indices, pos, i, j, k,   ATMOS_AA_, N, M, L);

                // T rows have a dependency on aux rows
                // ATMOS_TT_-ATMOS_PP_
                for (int aa = 1; aa <= aux_; ++aa)
                    indices[pos++] = last + aa;

                // Insert dependencies in matrixGraph
                CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                               gid0 + ATMOS_TT_, pos, indices));

                // Q-equation
                pos = 0;

                // ATMOS_QQ_-ATMOS_QQ_: 5-point stencil
                insert_graph_entry(indices, pos, i, j, k,   ATMOS_QQ_, N, M, L);
                insert_graph_entry(indices, pos, i-1, j, k, ATMOS_QQ_, N, M, L);
                insert_graph_entry(indices, pos, i+1, j, k, ATMOS_QQ_, N, M, L);
                insert_graph_entry(indices, pos, i, j-1, k, ATMOS_QQ_, N, M, L);
                insert_graph_entry(indices, pos, i, j+1, k, ATMOS_QQ_, N, M, L);

                // Add the dependencies on auxiliary unknowns
                // ATMOS_QQ_-ATMOS_PP_
                for (int aa = 1; aa <= aux_; ++aa)
                    indices[pos++] = last + aa;

                // Skip the final insertion when we are at rowIntCon_
                if ( (gid0 + ATMOS_QQ_) == rowIntCon_ && useIntCondQ_)
                {
                    // do nothing
                }
                else
                {
                    // Insert dependencies in matrixGraph
                    CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                                   gid0 + ATMOS_QQ_, pos, indices));
                }

                // A-equation
                pos = 0;

                // ATMOS_AA_-ATMOS_AA_: no spatial dependencies
                insert_graph_entry(indices, pos, i, j, k, ATMOS_AA_, N, M, L);

                // ATMOS_AA_-ATMOS_TT_: no spatial dependencies
                insert_graph_entry(indices, pos, i, j, k, ATMOS_TT_, N, M, L);

                // ATMOS_AA_-ATMOS_PP_
                for (int aa = 1; aa <= aux_; ++aa)
                    indices[pos++] = last + aa;

                // Insert dependencies in matrixGraph
                CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                               gid0 + ATMOS_AA_, pos, indices));
            }

    // Create graph entries for integral condition row
    if (standardMap_->MyGID(rowIntCon_) && useIntCondQ_)
    {
        int len = n_ * m_ * l_ + aux_;
        int icinds[len];
        int gcid;
        pos = 0;
        for (int k = 0; k != l_; ++k)
            for (int j = 0; j != m_; ++j)
                for (int i = 0; i != n_; ++i)
                {
                    gcid = FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_, i, j, k, ATMOS_QQ_);
                    icinds[pos] = gcid;
                    pos++;
                }

        // Although they effectively do not exist, we still need to
        // prepare dependencies on auxiliary unknowns as these are
        // needed when filling the subdomain Jacobians.
        for (int aa = 1; aa <= aux_; ++aa)
        {
            icinds[pos] = last + aa;
            pos++;
        }

        assert(len == pos);
        CHECK_NONNEG(matrixGraph_->InsertGlobalIndices(rowIntCon_, len, icinds));
    }

    // Dependencies of the auxiliary unknown, on the same proc as integral condition
    if ( standardMap_->MyGID(rowIntCon_) )
    {
        for (int aa = 1; aa <= aux_; ++aa)
        {
            int len = n_ * m_ * l_ + aux_;
            int auxinds[len];
            int gqid;
            pos = 0;
            for (int k = 0; k != l_; ++k)
                for (int j = 0; j != m_; ++j)
                    for (int i = 0; i != n_; ++i)
                    {
                        gqid = FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_, i, j, k, ATMOS_QQ_);
                        auxinds[pos] = gqid;
                        pos++;
                    }

            // Add the dependencies on auxiliary unknowns
            for (int aa = 1; aa <= aux_; ++aa)
                auxinds[pos++] = last + aa;

            assert(len == pos);

            CHECK_NONNEG(matrixGraph_->InsertGlobalIndices(last + aa, len, auxinds));
        }
    }

    // Finalize matrixgraph
    CHECK_ZERO(matrixGraph_->FillComplete() );
}

//=============================================================================
// Copied from THCM, adjusted for Atmosphere
void Atmosphere::insert_graph_entry(int* indices, int& pos,
                                       int i, int j, int k, int xx,
                                       int N, int M, int L) const
{
    int ii = i; // if x-boundary is periodic i may be out of bounds.
    // ii will be adjusted in that case:
    if (domain_->IsPeriodic())
    {
        ii = MOD((double)i, (double)N);
    }
    if ((ii>=0) && (j>=0) && (k>=0) &&
        (ii< N) && (j< M) && (k< L) )
    {
        // find index using 0-based find_row
        indices[pos++] = FIND_ROW_ATMOS0(ATMOS_NUN_, N, M, L, ii, j, k, xx);
    }
}

//=============================================================================
void Atmosphere::additionalExports(EpetraExt::HDF5 &HDF5, std::string const &filename)
{
    Atmosphere::CommPars pars;
    getCommPars(pars);

    // Get and write dimensional evaporation and precipitation fields
    Teuchos::RCP<Epetra_Vector> E = getE();
    Teuchos::RCP<Epetra_Vector> P = getP();

    HDF5.Write("E", *E);
    HDF5.Write("P", *P);

    // E->Update(-1.0, *P, 1.0); // E-P
    // std::cout << " int E - P dA = " << Utils::dot(pIntCoeff_, E) << std::endl;

    // Write surface temperatures
    getLandTemperature();
    HDF5.Write("lst", *lst_);
    HDF5.Write("sst", *sst_);
    HDF5.Write("sit", *sit_);
    HDF5.Write("Msi", *Msi_);

    // Write dimensional fluxes
    std::vector<Teuchos::RCP<Epetra_Vector> > fluxes = getFluxes();
    HDF5.Write("LongwaveFlux",     *fluxes[AtmosLocal::_QLW]);
    HDF5.Write("ShortwaveFlux",    *fluxes[AtmosLocal::_QSW]);
    HDF5.Write("SensibleHeatFlux", *fluxes[AtmosLocal::_QSH]);
    HDF5.Write("LatentHeatFlux",   *fluxes[AtmosLocal::_QLH]);

    if (saveMask_)
    {
        HDF5.Write("MaskGlobal", "Surface", H5T_NATIVE_INT,
                   surfmask_->size(), &(*surfmask_)[0]);
    }
}

//=============================================================================
std::vector<Teuchos::RCP<Epetra_Vector> > Atmosphere::getFluxes()
{
    std::vector<Teuchos::RCP<Epetra_Vector> > fluxes;
    std::vector<Epetra_Vector> localFluxes;

    int numFluxes = AtmosLocal::_QLH+1;
    for (int i = 0; i != numFluxes; ++i)
    {
        fluxes.push_back(Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_)));
        localFluxes.push_back(Epetra_Vector(*assemblySurfaceMap_));
    }

    std::vector<double*> tmpPtrs(numFluxes);

    for (int i = 0; i != numFluxes; ++i)
    {
        localFluxes[i].ExtractView(&tmpPtrs[i]);
    }

    atmos_->getFluxes(tmpPtrs[AtmosLocal::_QLW],
                      tmpPtrs[AtmosLocal::_QSW],
                      tmpPtrs[AtmosLocal::_QSH],
                      tmpPtrs[AtmosLocal::_QLH]);

    for (int i = 0; i != numFluxes; ++i)
    {
        CHECK_ZERO(fluxes[i]->Export(localFluxes[i], *as2std_surf_, Zero));
    }

    return fluxes;
}

//=============================================================================
int Atmosphere::interface_row(int i, int j, int XX)
{
    if (XX >= ATMOS_PP_)
    {
        if (aux_ <= 0)
        {
            WARNING("Atmosphere: invalid row", __FILE__, __LINE__);
            return -1;
        }
        else
        {
            // return rows beyond the regular n*m*l*dof
            return dim_ - aux_ + (XX - ATMOS_PP_);
        }
    }
    else
        return FIND_ROW_ATMOS0(dof_, n_, m_, l_, i, j, l_-1, XX);
}
