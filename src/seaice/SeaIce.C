#include "SeaIce.H"
#include "Ocean.H"
#include "Atmosphere.H"

extern "C" _SUBROUTINE_(getdeps)(double*, double*, double*,
                                 double*, double*, double*, double *);

//=============================================================================
// Constructor
SeaIce::SeaIce(Teuchos::RCP<Epetra_Comm> comm, ParameterList params)
    :
    params_          (params),
    nGlob_           (params->get("Global Grid-Size n", 16)),
    mGlob_           (params->get("Global Grid-Size m", 16)),
    periodic_        (params->get("Periodic", false)),
    
    precInitialized_ (false),
    recomputePrec_   (false),
    
    taus_         (0.01),    // threshold ice thickness
    epsilon_      (1.0e-2),  // Heavyside approximation steepness

    // background mean values
    t0o_   (params->get("background ocean temp t0o", 15)),
    t0a_   (params->get("background atmos temp t0a", 15)),
    s0_    (params->get("ocean background salinity s0", 35)),
    q0_    (params->get("atmos reference humidity",8e-3)),
    qdim_  (params->get("atmos humidity scale", 1e-3)),
    tdim_  (params->get("temperature scale", 1.0)),
    H0_    (params->get("seaice background thickness H0", taus_)),
    M0_    (params->get("seaice background mask M0", 0)),

    // ice formation constants
    ch_    (params->get("empirical constant", 0.0058)),
    utau_  (params->get("skin friction velocity, ms^{-1}", 0.02)),
    rhoo_  (params->get("sea water density, kg m^{-3}", 1.024e3)),
    rhoi_  (params->get("ice density, kg m^{-3}", 0.913e3)),
    rhoa_  (params->get("atmospheric density, kg m^{-3}", 1.25)),
    cpo_   (params->get("sea water heat capacity, W s kg^{-1} K^{-1]", 4.2e3)),
    Lf_    (params->get("latent heat of fusion of ice, J kg^{-1}", 3.347e5)),
    Ls_    (params->get("latent heat of sublimation of ice, J kg^{-1}", 2.835e6)),
    Ic_    (params->get("constant ice conductivity, W m^{-1} K^{-1}", 2.166)),
    a0_    (params->get("freezing temperature sensitivity", -0.0575)),

    // combined parameter
    zeta_  (ch_ * utau_ * rhoo_ * cpo_),

    // sublimation and evaporation constants, parameters for
    // saturation humidity over ice
    c1_    (params->get("c1", 3.8e-3)),
    c2_    (params->get("c2", 21.87)),
    c3_    (params->get("c3", 265.5)),
    c4_    (params->get("c4", 17.67)),
    c5_    (params->get("c5", 243.5)),

    ce_    (params->get("Dalton number", 1.3e-03)),
    uw_    (params->get("mean atmospheric surface wind speed, ms^{-1}", 8.5)),

// typical vertical velocity
    eta_   (( rhoa_ / rhoo_ ) * ce_ * uw_),

// Shortwave radiation constants and functions
    albe0_   (params->get("reference albedo", 0.3)),
    albed_   (params->get("albedo excursion", 0.5)),
    sun0_    (params->get("solar constant", 1360)),
    c0_      (params->get("atmospheric absorption coefficient", 0.43)),
    Ch_      (params->get("Ch", 1.22e-3)),
    cpa_     (params->get("heat capacity", 1000)),

// exchange coefficient
    muoa_    (rhoa_ * Ch_ * cpa_ * uw_)
{
    INFO("SeaIce constructor");
    // Continuation parameters
    allParameters_ = { "Combined Forcing",
                       "Solar Forcing",
                       "Latent Heat Forcing",
                       "Mask Forcing" };

    parName_ = params->get( "Continuation parameter",
                             allParameters_[0] );

    comb_   = params->get(allParameters_[0], 0.0);
    sunp_   = params->get(allParameters_[1], 1.0);
    latf_   = params->get(allParameters_[2], 1.0);
    maskf_  = params->get(allParameters_[3], 1.0);

    // inherited input/output datamembers
    inputFile_  = params->get("Input file",  "seaice_input.h5");
    outputFile_ = params->get("Output file", "seaice_output.h5");
    loadState_  = params->get("Load state", false);
    saveState_  = params->get("Save state", true);

    // set communicator 
    comm_ = comm;
    
    // Background sublimation and derivatives: calculate background
    // saturation specific humidity according to [Bolton,1980], T in
    // \deg C
    qsi_ = [&] (double t0i)
        {
            return c1_ * exp(c2_ * t0i / (t0i + c3_) );
        };

    qso_ = [&] (double t0o)
        {
            return c1_ * exp(c4_ * t0o / (t0o + c5_) );
        };

    dqsi_ = [&] (double t0i)
        {
            return (c1_ * c2_ * c3_) / pow(t0i + c3_, 2) *
            exp( (c2_ * t0i) / (t0i + c3_) );
        };

    dqso_ = [&] (double t0o)
        {
            return (c1_ * c4_ * c5_) / pow(t0o + c5_, 2) *
            exp( (c4_ * t0o) / (t0o + c5_) );
        };

    // Background sea ice surface temperature is chosen such that
    // background evaporation and sublimation cancel:
    t0i_  =  c3_*c4_*t0o_ / (c2_*c5_+(c2_-c4_)*t0o_);

    // Background sublimation and derivatives
    E0i_   =  eta_ * ( qsi_(t0i_) - q0_ );
    E0o_   =  eta_ * ( qso_(t0o_) - q0_ );

    // Check whether background values cancel
    assert(std::abs(E0o_-E0i_) < 1e-12);

    dEdT_  =  eta_ * qdim_ * tdim_ / qdim_ * dqsi_(t0i_);
    dEdq_  =  eta_ * qdim_ * -1;

    //! Ocean nondimensionalization prefactor (set in synchronization)
    pQSnd_ = 1.0;

    // background heat flux variation
    Qvar_  = zeta_;

    // background heat flux
    // Q0_    = zeta_ * (freezingT(0) - t0o_) - rhoo_ * Lf_ * E0_;
    Q0_ = -100.0;

    xmin_ = params->get("Global Bound xmin", 286.0) * PI_ / 180.0;
    xmax_ = params->get("Global Bound xmax", 350.0) * PI_ / 180.0;
    ymin_ = params->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
    ymax_ = params->get("Global Bound ymax", 80.0)  * PI_ / 180.0;

    dof_      = SEAICE_NUN_;
    dimGlob_  = mGlob_ * nGlob_ * dof_;

    // auxiliary integral correction
    aux_      = 1;
    
    // Create domain object
    domain_ = Teuchos::rcp(new TRIOS::Domain(nGlob_, mGlob_, 1, dof_,
                                             xmin_, xmax_, ymin_, ymax_,
                                             periodic_, 1.0, comm_, aux_));

    // Compute 2D decomposition
    domain_->Decomp2D();

    // local dimensions
    xminLoc_ = domain_->XminLoc();
    xmaxLoc_ = domain_->XmaxLoc();
    yminLoc_ = domain_->YminLoc();
    ymaxLoc_ = domain_->YmaxLoc();

    INFO("   local sea ice model: local xmin = " << xminLoc_);
    INFO("                        local xmax = " << xmaxLoc_);
    INFO("                        local ymin = " << yminLoc_);
    INFO("                        local ymax = " << ymaxLoc_);

    // local grid dimensions
    nLoc_   =  domain_->LocalN();
    mLoc_   =  domain_->LocalM();
    dimLoc_ =  mLoc_ * nLoc_ * dof_;

    // Local flux containers
    QSos_ = std::vector<double>(mLoc_ * nLoc_);
    EmiP_ = std::vector<double>(mLoc_ * nLoc_);

    // initialize mask
    surfmask_ = std::make_shared<std::vector<int> >(mGlob_ * nGlob_);

    // Obtain overlapping and non-overlapping maps
    assemblyMap_ = domain_->GetAssemblyMap(); // overlapping
    standardMap_ = domain_->GetStandardMap(); // non-overlapping

    // Obtain special maps: depth-averaged, single unknown
    standardSurfaceMap_ = domain_->GetStandardSurfaceMap();
    assemblySurfaceMap_ = domain_->GetAssemblySurfaceMap();

    // Create vectors
    state_      = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    rhs_        = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    diagB_      = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    sol_        = Teuchos::rcp(new Epetra_Vector(*standardMap_));

    // External data
    sst_        = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    sss_        = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    tatm_       = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    qatm_       = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    patm_       = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    albe_       = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    intCoeff_   = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

    // Local (overlapping) vectors
    localState_  = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localRHS_    = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localDiagB_  = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localSol_    = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));

    // External data (overlapping)
    localSST_      = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localSSS_      = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosT_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosQ_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosP_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosA_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localIntCoeff_ = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    
    localSurfmask_ = Teuchos::rcp(new Epetra_IntVector(*assemblySurfaceMap_));

    // Create local computational grid in x_, y_
    createGrid();

    // Create integral coefficients intCoeff_ and localIntCoeff_
    createIntCoeff();

    // Construct local dependency grid:
    Al_ = std::make_shared<DependencyGrid>(nLoc_, mLoc_, 1, 1, dof_ + aux_);

    createMatrixGraph();

    // Initialize Jacobian
    jac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *matrixGraph_));

    // Import existing state
    if (loadState_)
        loadStateFromFile(inputFile_);
    else
        initializeState();

    // Create importers for communication with other models
    int XX;
    for (int i = 0; i < (dof_ + aux_); ++i)
    {
        XX = SEAICE_HH_ + i;
        Maps_[XX] = Utils::CreateSubMap(*standardMap_, dof_, XX);
        Imps_[XX] = Teuchos::rcp(new Epetra_Import(*Maps_[XX], *standardMap_));
    }
    INFO("SeaIce constructor done");    
}

//=============================================================================
void SeaIce::computeRHS()
{
    // zero rhs vector
    localRHS_->PutScalar(0.0);

    // obtain view of rhs, state and external data
    double *rhs, *state, *sst, *sss;
    double *tatm, *qatm, *patm, *albe, *flxd;
    localRHS_->ExtractView(&rhs);

    domain_->Standard2Assembly(*state_, *localState_);

    domain_->Standard2AssemblySurface(*sst_,   *localSST_);
    domain_->Standard2AssemblySurface(*sss_,   *localSSS_);
    domain_->Standard2AssemblySurface(*tatm_,  *localAtmosT_);
    domain_->Standard2AssemblySurface(*qatm_,  *localAtmosQ_);
    domain_->Standard2AssemblySurface(*patm_,  *localAtmosP_);
    domain_->Standard2AssemblySurface(*albe_,  *localAtmosA_);

    localState_->ExtractView(&state);

    localSST_->ExtractView(&sst);
    localSSS_->ExtractView(&sss);

    localAtmosT_->ExtractView(&tatm);
    localAtmosQ_->ExtractView(&qatm);
    localAtmosP_->ExtractView(&patm);
    localAtmosA_->ExtractView(&albe);

    if (aux_ == 1)
        computeLocalFluxes(state, sss, sst, qatm, patm);
            
    // Create assembly and standard vectors to hold flux
    /// difference for integral correction
    Epetra_Vector fluxDiff(*standardSurfaceMap_);
    Epetra_Vector localFluxDiff(*assemblySurfaceMap_);
    localFluxDiff.ExtractView(&flxd);

    // row indices for rhs and data
    int rr, sr;
    double Tsi, Hval, Qval, Mval, Tval, Gval, val;
    for (int j = 0; j != mLoc_; ++j)
        for (int i = 0; i != nLoc_; ++i)
        {
            sr = j*nLoc_ + i; // surface position for vectors with dof=1

            Hval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_HH_)];
            Qval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_QQ_)];
            Mval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_MM_)];
            Tval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_TT_)];

            // Although T is part of the state, we substitute the sea
            // ice temperature equation in the H equation to improve
            // diagonal dominance in these rows.
            Tsi = iceSurfT(Qval, Hval, sss[sr]);

            for (int XX = 1; XX <= dof_; ++XX)
            {
                switch (XX)
                {

                case SEAICE_HH_:   // H row (thickness)

                    val = freezingT(sss[sr]) - sst[sr] - t0o_ -
                        ( Q0_ / zeta_ +  Qvar_ / zeta_ * Qval ) -
                        ( rhoo_ * Lf_ / zeta_ ) *
                        ( E0i_ + dEdT_ * Tsi + dEdq_ * qatm[sr] );
                    
                    break;

                case SEAICE_QQ_:   // Q row (heat flux)

                    val = 1.0 / muoa_ * Q0_ + Qvar_ / muoa_ * Qval -
                        (comb_ * sunp_ * sun0_ / 4. / muoa_) * shortwaveS(y_[j]) *
                        (1. - albe0_ - albed_*albe[sr]) * c0_ +
                        (Tval + t0i_ - tatm[sr] - t0a_) +
                        (comb_ * latf_ * rhoo_ * Ls_ / muoa_) *
                        (E0i_ + dEdT_ * Tval +  dEdq_ * qatm[sr]);

                    break;

                case SEAICE_MM_:   // M row (mask)

                    val = Mval - maskFun(Hval);
                    
                    break;

                case SEAICE_TT_:   // T row (surface temperature)

                    val = Tsi - Tval;

                    break;
                }

                rr = find_row0(nLoc_, mLoc_, i, j, XX);
//                if ((*localSurfmask_)[sr] == 0)
                    rhs[rr] = val;
//                else
//                    rhs[rr] = 0.0;
            }

            if (aux_ == 1)
            {
                // Flux difference over sea ice                
                flxd[sr] = Mval * (QSos_[sr] - EmiP_[sr]);
            }
        }
    
    domain_->Assembly2Standard(*localRHS_, *rhs_);

    if (aux_ == 1)
    {
        // create non-overlapping flux difference
        domain_->Assembly2StandardSurface(localFluxDiff, fluxDiff);

        // integrate flux difference over total surface area
        double fluxInt;
        intCoeff_->Dot(fluxDiff, &fluxInt);

        // implement integral equation in RHS
        int Grow = find_row0(nGlob_, mGlob_, 0, 0, SEAICE_GG_);
        int lid  = standardMap_->LID(Grow);
        if (lid >= 0)
        {
            Gval = state[lid];
            (*rhs_)[lid] = pQSnd_ * fluxInt - Gval * totalArea_;
        }
    }
    
    INFO(" seaic F = " << Utils::norm(rhs_));
}

//=============================================================================
void SeaIce::computeLocalFluxes(double *state, double *sss, double *sst,
                                double *qatm, double *patm)

{
    assert((int) QSos_.size() == nLoc_ * mLoc_);
    int sr;
    double Qval, Tval;
    for (int j = 0; j != mLoc_; ++j)
        for (int i = 0; i != nLoc_; ++i)
        {
            sr   = j*nLoc_ + i;
            
            Qval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_QQ_)];
            Tval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_TT_)];

            QSos_[sr] = (
                zeta_ * ( freezingT(sss[sr])    // QTos component
                          - (sst[sr]+t0o_))
                - (Qvar_ * Qval +  Q0_)         // QTsa component
                ) / rhoo_ / Lf_;
            
            EmiP_[sr] = dEdT_ * Tval + dEdq_ * qatm[sr] // E component
                - eta_ * qdim_ * patm[sr];              // P component
            
        }            
}

//=============================================================================
void SeaIce::computeLocalJacobian()
{
    int H = SEAICE_HH_;
    int Q = SEAICE_QQ_;
    int M = SEAICE_MM_;
    int T = SEAICE_TT_;
    int G = SEAICE_GG_;

    // reset entries
    Al_->zero();

    // obtain local state
    domain_->Standard2Assembly(*state_, *localState_);
    
    domain_->Standard2AssemblySurface(*sst_,   *localSST_);
    domain_->Standard2AssemblySurface(*sss_,   *localSSS_);
    domain_->Standard2AssemblySurface(*qatm_,  *localAtmosQ_);
    domain_->Standard2AssemblySurface(*patm_,  *localAtmosP_);

    double *state, *sst, *sss, *qatm, *patm;
    localState_->ExtractView(&state);
    localSST_->ExtractView(&sst);
    localSSS_->ExtractView(&sss);
    localAtmosQ_->ExtractView(&qatm);
    localAtmosP_->ExtractView(&patm);
    
    if (aux_ == 1)
        computeLocalFluxes(state, sss, sst, qatm, patm);

    // our range is the entire local domain (1-based)
    int range[8] = {1, nLoc_, 1, mLoc_, 1, 1, 1, 1};

    // initialize dependencies
    double HH_HH, HH_QQ, QQ_QQ, QQ_TT, MM_MM;
    double TT_HH, TT_QQ, TT_TT, GG_GG;
    Atom MM_HH(nLoc_, mLoc_, 1, 1);
    Atom GG_QQ(nLoc_, mLoc_, 1, 1);        
    Atom GG_MM(nLoc_, mLoc_, 1, 1);
    Atom GG_TT(nLoc_, mLoc_, 1, 1);

    // dHdt equation ----------------------------------
    HH_HH = -rhoo_ * Lf_ / zeta_ * dEdT_ * Q0_ / Ic_;
    HH_QQ = -Qvar_ / zeta_ -
        rhoo_ * Lf_ / zeta_ * dEdT_ * H0_ * Qvar_ / Ic_;

    Al_->set(range, H, H, HH_HH);
    Al_->set(range, H, Q, HH_QQ);

    // Qtsa equation ---------------------------------
    QQ_QQ = Qvar_ / muoa_;
    QQ_TT = 1.0 + comb_ * latf_ * rhoo_ * Ls_ / muoa_ * dEdT_;

    Al_->set(range, Q, Q, QQ_QQ);
    Al_->set(range, Q, T, QQ_TT);

    // Msi equation -----------------------------------
    // fill atom for nonlinear contribution H
    int Hrow;     // state index
    double val;
    for (int j = 0; j < mLoc_; ++j)
        for (int i = 0; i < nLoc_; ++i)
        {
            Hrow  = find_row0(nLoc_, mLoc_, i, j, H); // H row

            val = -dMdH(state[Hrow]);

            MM_HH.set(i+1, j+1, 1, 1, val); // Atoms are 1-based
        }

    MM_MM = 1.0;

    Al_->set(range, M, H, MM_HH);
    Al_->set(range, M, M, MM_MM);

    // Tsi equation ------------------------------------
    TT_HH = Q0_ / Ic_;
    TT_QQ = H0_ * Qvar_ / Ic_;
    TT_TT = -1.0;

    Al_->set(range, T, H, TT_HH);
    Al_->set(range, T, Q, TT_QQ);
    Al_->set(range, T, T, TT_TT);

    // Gamma integral equation ----------------------------------
    if (aux_ == 1)
    {
        int sr, Mrow;
        double Mval, ICval;
        for (int j = 0; j < mLoc_; ++j)
            for (int i = 0; i < nLoc_; ++i)
            {
                sr    = j*nLoc_ + i;
                Mrow  = find_row0(nLoc_, mLoc_, i, j, M); // M row
                Mval  = state[Mrow];
                ICval = (*localIntCoeff_)[sr] * pQSnd_;

                // GG_QQ
                val  = -1.0 * ICval * Mval * Qvar_ / rhoo_ / Lf_;
                GG_QQ.set(i+1, j+1, 1, 1, val); // Atoms are 1-based
            
                // GG_MM
                val  = ICval * (QSos_[sr] - EmiP_[sr]);
                GG_MM.set(i+1, j+1, 1, 1, val);

                // GG_TT
                val  = -1.0 * ICval * Mval * dEdT_;
                GG_TT.set(i+1, j+1, 1, 1, val); // Atoms are 1-based
            }
    
        // diagonal dependence
        GG_GG = -totalArea_;

        Al_->set(range, G, Q, GG_QQ);
        Al_->set(range, G, M, GG_MM);
        Al_->set(range, G, T, GG_TT);
        Al_->set(range, G, G, GG_GG);
    }

    //---------------------------------------------------------------
    // Set land points to 1, RHS elements should be 0 at land points
    // int sr;
    // for (int j = 0; j < mLoc_; ++j)
    //     for (int i = 0; i < nLoc_; ++i)
    //         for (int A = 1; A <= (dof_+aux_); ++A)
    //             for (int B = 1; B <= (dof_+aux_); ++B)
    //                 {
    //                     sr = j*nLoc_ + i;       // local surface index
    //                     if ((*localSurfmask_)[sr])
    //                     {
    //                         if (B==A)
    //                             Al_->set(i+1,j+1,1,1,A,B,1.0);
    //                         else
    //                             Al_->set(i+1,j+1,1,1,A,B,0.0);
    //                     }
    //                 }
}

//=============================================================================
void SeaIce::computeJacobian()
{

    // set all entries to zero
    CHECK_ZERO(jac_->PutScalar(0.0));

    // Create local dependency grid
    computeLocalJacobian();

    // Create local crs arrays
    assemble();

    // Obtain crs arrays
    std::shared_ptr<Utils::CRSMat> localJac = getLocalJacobian();

    //--------------------------------------------
    // Create parallel crs matrix
    //--------------------------------------------

    // max nonzeros per row
    const int maxnnz = dof_ + 1;
    Utils::assembleCRS(jac_, *localJac, maxnnz, domain_);

    jac_->FillComplete();

    // With a new Jacobian we need to recompute the factorization
    recomputePrec_ = true;
}

//=============================================================================
void SeaIce::assemble()
{
    // Assemble local Al_ into local crs vectors
    // clear old CRS matrix
    beg_.clear();
    co_.clear();
    jco_.clear();

    // We do this 1-based
    int elm_ctr = 1, col;
    double value;
    for (int j = 1; j <= mLoc_; ++j)
        for (int i = 1; i <= nLoc_; ++i)
            for (int A = 1; A <= dof_; ++A)
            {
                // fill beg with element cntr
                beg_.push_back(elm_ctr);

                for (int B = 1; B <= (dof_ + aux_); ++B)
                {
                    value = Al_->get( i, j, 1, 1, A, B );

                    if (std::abs(value) > 0)
                    {
                        co_.push_back(value);

                        // obtain column
                        col = find_row1(nLoc_, mLoc_,  i, j, B);
                        jco_.push_back(col);
                        ++elm_ctr;
                    }
                }
            }

    // auxiliary equation
    if (aux_ == 1)
    {
        beg_.push_back(elm_ctr);
        for (int j = 1; j <= mLoc_; ++j)
            for (int i = 1; i <= nLoc_; ++i)
                for (int B = 1; B <= (dof_ + aux_); ++B)
                {
                    value = Al_->get(i, j, 1, 1, SEAICE_GG_, B);
                    if (std::abs(value) > 0)
                    {
                        co_.push_back(value);
                        
                        // obtain column
                        col = find_row1(nLoc_, mLoc_,  i, j, B);
                        jco_.push_back(col);
                        ++elm_ctr;
                    }
                }
    }
    
    // final element of beg
    beg_.push_back(elm_ctr);
}

//=============================================================================
std::shared_ptr<Utils::CRSMat> SeaIce::getLocalJacobian()
{
    std::shared_ptr<Utils::CRSMat> jac = std::make_shared<Utils::CRSMat>();
    jac->co  = co_;
    jac->jco = jco_;
    jac->beg = beg_;
    return jac;
}

//=============================================================================
void SeaIce::getCommPars(SeaIce::CommPars &parStruct)
{
    parStruct.zeta = zeta_;
    parStruct.a0   = a0_;
    parStruct.Lf   = Lf_;
    parStruct.s0   = s0_;
    parStruct.rhoo = rhoo_;
    parStruct.Qvar = Qvar_;
    parStruct.Q0   = Q0_;

}

// ---------------------------------------------------------------------------
// Adjust locally defined parameter
double SeaIce::getPar()
{
    return getPar(parName_);
}

// ---------------------------------------------------------------------------
// Get parameter value
double SeaIce::getPar(std::string const &parName)
{
    if (parName.compare(allParameters_[0]) == 0)
        return comb_;
    else if (parName.compare(allParameters_[1]) == 0)
        return sunp_;
    else if (parName.compare(allParameters_[2]) == 0)
        return latf_;
    else if (parName.compare(allParameters_[3]) == 0)
        return maskf_;
    else // If parameter not available we return 0
        return 0;
}

// ---------------------------------------------------------------------------
// Set continuation parameter
void SeaIce::setPar(double value)
{
    setPar(parName_, value);
}

//-----------------------------------------------------------------------------
int SeaIce::npar()
{
    return (int) allParameters_.size();
}

//-----------------------------------------------------------------------------
std::string const SeaIce::int2par(int ind)
{
    return allParameters_[ind];
}

//-----------------------------------------------------------------------------
void SeaIce::setLandMask(Utils::MaskStruct const &mask)
{
    // create global surface mask
    surfmask_->clear();

    if ((int) mask.global_surface->size() < (mGlob_ * nGlob_))
    {
        ERROR("mask.global_surface->size() not ok:",  __FILE__, __LINE__);
    }
    else // we trust surfm
    {
        surfmask_ = mask.global_surface;
    }

    int lsr, gsr;
    for (int j = 0; j != mLoc_; ++j)
        for (int i = 0; i != nLoc_; ++i)
        {
            lsr = j*nLoc_ + i;
            gsr = assemblySurfaceMap_->GID(lsr);
            (*localSurfmask_)[lsr] = (*surfmask_)[gsr];
        }
               
    // adjust integral coefficients
    createIntCoeff();
}

// ---------------------------------------------------------------------------
// Set specific continuation parameter
void SeaIce::setPar(std::string const &parName, double value)
{
    parName_ = parName; // Overwrite our parameter name

    if (parName.compare(allParameters_[0]) == 0)
        comb_  = value;
    else if (parName.compare(allParameters_[1]) == 0)
        sunp_  = value;
    else if (parName.compare(allParameters_[2]) == 0)
        latf_  = value;
    else if (parName.compare(allParameters_[3]) == 0)
        maskf_ = value;
}

//=============================================================================
std::shared_ptr<Utils::CRSMat> SeaIce::getBlock(std::shared_ptr<Atmosphere> atmos)
{
    // initialize empty CRS matrix
    std::shared_ptr<Utils::CRSMat> block = std::make_shared<Utils::CRSMat>();

    // construct global 0-based CRS matrix
    int el_ctr = 0;

    int T = ATMOS_TT_; // (1-based) in the Atmosphere, temperature is the first unknown
    int Q = ATMOS_QQ_; // (1-based) in the Atmosphere, humidity is the second unknown
    int A = ATMOS_AA_; // (1-based) in the Atmosphere, albedo is the third unknown
    int P = ATMOS_PP_; // (1-based) in the Atmosphere, precipitation is auxiliary

    // FIXME: build this block locally to avoid Msi gather
    Teuchos::RCP<Epetra_Vector> Msi = interfaceM();
    Teuchos::RCP<Epetra_MultiVector> MsiG = Utils::AllGather(*Msi);

    // compute a few constant derivatives (see computeRHS)
    // d / dq_atm (F_H)
    double dqatmFH = -(rhoo_ * Lf_ / zeta_) * dEdq_;

    // d / dt_atm (F_Q)
    double dtatmFQ = -1.0;

    // d / dq_atm (F_Q)
    double dqatmFQ =  (comb_ * latf_ * rhoo_ * Ls_ / muoa_) * dEdq_;

    // d / da_atm (F_Q)
    double daatmFQ;
    double tmp = 0.0;
    // int sr;
    int col;

    for (int j = 0; j != mGlob_; ++j)
    {
        int gid = j * nGlob_;
        int lid = standardSurfaceMap_->LID(gid);

        if (lid >= 0)
            tmp = (comb_ * sunp_ * sun0_ / 4. / muoa_) *
                shortwaveS(y_[lid / nLoc_]) * albed_ * c0_;

        comm_->SumAll( &tmp, &daatmFQ, 1);

        for (int i = 0; i != nGlob_; ++i)
            for (int XX = 1; XX <= dof_; ++XX)
            {
                block->beg.push_back(el_ctr);
                
//                sr = j*nGlob_ + i;            // global surface index
//                if ((*surfmask_)[sr] == 0)
                {
                    switch (XX)
                    {
                    case SEAICE_HH_:
                        block->co.push_back(dqatmFH);
                        block->jco.push_back(atmos->interface_row(i,j,Q));
                        el_ctr++;
                        break;

                    case SEAICE_QQ_:
                        block->co.push_back(dtatmFQ);
                        block->jco.push_back(atmos->interface_row(i,j,T));
                        el_ctr++;

                        block->co.push_back(dqatmFQ);
                        block->jco.push_back(atmos->interface_row(i,j,Q));
                        el_ctr++;

                        block->co.push_back(daatmFQ);
                        block->jco.push_back(atmos->interface_row(i,j,A));
                        el_ctr++;
                        break;
                    }
                }
            }
    }

    // auxiliary equation
    if (aux_ == 1)
    {
        int sr;
        double dQFG; // d / dQ (F_G)
        double dPFG; // d / dQ (F_G)
        double ICval, Mval;
        block->beg.push_back(el_ctr);
        for (int j = 0; j != mGlob_; ++j)
            for (int i = 0; i != nGlob_; ++i)
            {
                sr    = j*nGlob_ + i;            // global surface index

//                if ((*surfmask_)[sr] == 0)
                {
                    ICval = (*globalIntCoeff_)[sr];  // integral coefficient
                    Mval  = (*(*MsiG)(0))[sr];        // mask value

                    dQFG  = Mval * ICval * pQSnd_ * (-dEdq_);
                    block->co.push_back(dQFG);
                    block->jco.push_back(atmos->interface_row(i,j,Q));
                    el_ctr++;

                }
            }
        col = atmos->interface_row(0,0,P);

        double totalM = Utils::dot(intCoeff_, Msi);

        if (col >= 0)
        {
            dPFG  = totalM * pQSnd_ * eta_ * qdim_;
            block->co.push_back(dPFG);
            block->jco.push_back(col);
            el_ctr++;
        }
    }

    // final entry in beg ( == nnz)
    block->beg.push_back(el_ctr);

    // final check;
    assert( (int) block->co.size() == block->beg.back() );

    return block;
}

//=============================================================================
std::shared_ptr<Utils::CRSMat> SeaIce::getBlock(std::shared_ptr<Ocean> ocean)
{
    // initialize empty CRS matrix
    std::shared_ptr<Utils::CRSMat> block = std::make_shared<Utils::CRSMat>();

    int el_ctr = 0;

    int T = 5; // (1-based) in the Ocean, temperature is the fifth unknown
    int S = 6; // (1-based) in the Ocean, salinity is the sixth unknown

    // compute a few constant derivatives (see computeRHS)

    // d / dT (F_H)
    double dTFH = -1.0;

    // d / dS (F_H)
    double dSFH =  a0_ - ( rhoo_ * Lf_ / zeta_) * dEdT_ * a0_;

    // d / dS (F_T)
    double dSFT =  a0_;

//    int sr;
    for (int j = 0; j != mGlob_; ++j)
        for (int i = 0; i != nGlob_; ++i)
            for (int XX = 1; XX <= dof_; ++XX)
            {
//                sr = j*nGlob_ + i;            // global surface index

                block->beg.push_back(el_ctr);

//                if ((*surfmask_)[sr] == 0)
                {
                    switch (XX)
                    {
                    case SEAICE_HH_:
                        block->co.push_back(dTFH);
                        block->jco.push_back(ocean->interface_row(i,j,T));
                        el_ctr++;

                        block->co.push_back(dSFH);
                        block->jco.push_back(ocean->interface_row(i,j,S));
                        el_ctr++;
                        break;

                    case SEAICE_TT_:
                        block->co.push_back(dSFT);
                        block->jco.push_back(ocean->interface_row(i,j,S));
                        el_ctr++;
                        break;
                    }
                }
            }

    // Auxiliary equation

    // FIXME: build this block locally to avoid Msi gather
    Teuchos::RCP<Epetra_MultiVector> MsiG = Utils::AllGather(*interfaceM());

    if (aux_ == 1)
    {
        int sr;
        double dTFG; // d / dTo (F_G)
        double dSFG; // d / dSo (F_G)
        double ICval, Mval; 
        block->beg.push_back(el_ctr);
        for (int j = 0; j != mGlob_; ++j)
            for (int i = 0; i != nGlob_; ++i)
            {
                sr    = j*nGlob_ + i;             // global surface index
                ICval = (*globalIntCoeff_)[sr];   // integral coefficient
                Mval  = (*(*MsiG)(0))[sr];        // mask value

//                if ((*surfmask_)[sr] == 0)
                {
                    dTFG  = Mval * ICval * pQSnd_ * zeta_ * -1.0 / rhoo_ / Lf_;
                    block->co.push_back(dTFG);
                    block->jco.push_back( ocean->interface_row(i,j,T) );
                    el_ctr++;

                    dSFG  = Mval * ICval * pQSnd_ * zeta_ * a0_ / rhoo_ / Lf_;
                    block->co.push_back(dSFG);
                    block->jco.push_back(ocean->interface_row(i,j,S));
                    el_ctr++;
                }
            }
    }

    // final entry in beg ( == nnz)
    block->beg.push_back(el_ctr);

    // final check;
    assert( (int) block->co.size() == block->beg.back());

    return block;
}

//=============================================================================
void SeaIce::synchronize(std::shared_ptr<Ocean> ocean)
{
    // Obtain surface ocean temperature
    Teuchos::RCP<Epetra_Vector> sst = ocean->interfaceT();
    CHECK_MAP(sst, standardSurfaceMap_);
    sst_ = sst;

    // Obtain surface ocean salinity
    Teuchos::RCP<Epetra_Vector> sss = ocean->interfaceS();
    CHECK_MAP(sss, standardSurfaceMap_);
    sss_ = sss;

    // Get ocean parameters
    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 ;
    FNAME(getdeps)(&tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &pQSnd_);

}

//=============================================================================
void SeaIce::synchronize(std::shared_ptr<Atmosphere> atmos)
{
    // get atmosphere temperature
    Teuchos::RCP<Epetra_Vector> tatm  = atmos->interfaceT();
    CHECK_MAP(tatm, standardSurfaceMap_);
    tatm_ = tatm;

    // get atmosphere humidity
    Teuchos::RCP<Epetra_Vector> qatm  = atmos->interfaceQ();
    CHECK_MAP(qatm, standardSurfaceMap_);
    qatm_ = qatm;

    // get albedo
    Teuchos::RCP<Epetra_Vector> albe  = atmos->interfaceA();
    CHECK_MAP(albe, standardSurfaceMap_);
    albe_ = albe;

    // get precip
    Teuchos::RCP<Epetra_Vector> patm  = atmos->interfaceP();
    CHECK_MAP(patm, standardSurfaceMap_);
    patm_ = patm;

    Atmosphere::CommPars atmosPars;
    atmos->getCommPars(atmosPars);

    albe0_ = atmosPars.a0;
    albed_ = atmosPars.da;
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> SeaIce::interface(int XX)
{
    Teuchos::RCP<Epetra_Vector> out = Teuchos::rcp(new Epetra_Vector(*Maps_[XX]));
    CHECK_ZERO(out->Import(*state_, *Imps_[XX], Insert));

    assert(out->MyLength() >= 1);
    
    if (out->MyLength() == 1) // auxiliary unknown, convert to field
    {
        Teuchos::RCP<Epetra_Vector> converted =
            Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
        
        converted->PutScalar((*out)[0]);
        return converted;
    }
    else
        return out;
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> SeaIce::interfaceH()
{
    return interface(SEAICE_HH_);
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> SeaIce::interfaceQ()
{
    return interface(SEAICE_QQ_);
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> SeaIce::interfaceM()
{
    return interface(SEAICE_MM_);
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> SeaIce::interfaceT()
{
    return interface(SEAICE_TT_);
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> SeaIce::interfaceG()
{
    if (aux_ == 1)
        return interface(SEAICE_GG_);
    else
        return Teuchos::null;
}

//=============================================================================
// Build idealized local and global forcing vectors
void SeaIce::idealizedForcing()
{
    // extract views of sea surface temp sst
    //                  sea surface salinity sss
    //                  atmosphere temp tatm
    //                  atmosphere humidity qatm

    // a few idealized variations
    double tvar = 15.0;
    double svar = 1.0;
    double qvar = 1e-3;

    double *sst, *sss, *tatm, *qatm, *albe, *patm;
    localSST_->ExtractView(&sst);
    localSSS_->ExtractView(&sss);
    localAtmosT_->ExtractView(&tatm);
    localAtmosQ_->ExtractView(&qatm);
    localAtmosA_->ExtractView(&albe);
    localAtmosP_->ExtractView(&patm);

    int row;
    for (int j = 0; j != mLoc_; ++j)
        for (int i = 0; i != nLoc_; ++i)
        {
            row      = j * nLoc_ + i;

            sst[row]  = tvar * cos( PI_ * y_[j] / ymax_ ) - 10.0;
            sss[row]  = svar * cos( PI_ * y_[j] / ymax_ ) / cos( y_[j] );
            tatm[row] = tvar * cos( PI_ * y_[j] / ymax_ ) - 10.0;
            qatm[row] = qvar * cos( PI_ * y_[j] / ymax_ );
            albe[row] = 0.0;
            patm[row] = 1.0;
        }

    // Transfer data to non-overlapping vectors
    domain_->Assembly2StandardSurface(*localSST_,     *sst_);
    domain_->Assembly2StandardSurface(*localSSS_,     *sss_);
    domain_->Assembly2StandardSurface(*localAtmosT_, *tatm_);
    domain_->Assembly2StandardSurface(*localAtmosQ_, *qatm_);
    domain_->Assembly2StandardSurface(*localAtmosA_, *albe_);
    domain_->Assembly2StandardSurface(*localAtmosP_, *patm_);
}

//=============================================================================
// Create local grid points
void SeaIce::createGrid()
{
    dx_  =  (xmaxLoc_ - xminLoc_) / nLoc_;
    dy_  =  (ymaxLoc_ - yminLoc_) / mLoc_;

    for (int i = 0; i != nLoc_; ++i)
        x_.push_back(xminLoc_ + (i + 0.5) * dx_);

    for (int j = 0; j != mLoc_; ++j)
        y_.push_back(yminLoc_ + (j + 0.5) * dy_);
}

//=============================================================================
void SeaIce::createIntCoeff()
{
    // creating local integral coefficients
    localIntCoeff_->PutScalar(0.0);
    intCoeff_->PutScalar(0.0);
    int idx = 0;
    int sr;
    for (int j = 0; j != mLoc_; ++j)
        for (int i = 0; i != nLoc_; ++i)
        {
            sr = j*nLoc_ + i;
            if ((*localSurfmask_)[sr] == 0)
                (*localIntCoeff_)[idx] = cos(y_[j]) * dx_ * dy_;
            idx++;
        }

    // distribute to non-overlapping vector
    domain_->Assembly2StandardSurface(*localIntCoeff_, *intCoeff_);

    // obtain total area
    intCoeff_->Norm1(&totalArea_);

    // create global gathered vector
    globalIntCoeff_ =
        Teuchos::rcp(new Epetra_Vector( *(*Utils::AllGather(*intCoeff_))(0) ));
}

//=============================================================================
// Create matrix graph
void SeaIce::createMatrixGraph()
{
    // We do not have any spatial relation between the unknowns so
    // maxDeps is at most the number of unknowns.
    int maxDeps = dof_ + aux_;

    matrixGraph_ =
        Teuchos::rcp(new Epetra_CrsGraph(Copy, *standardMap_, maxDeps, false));

    int indices[maxDeps];

    // Get our local range in all directions
    // 0-based
    int I0 = domain_->FirstRealI();
    int J0 = domain_->FirstRealJ();
    int I1 = domain_->LastRealI();
    int J1 = domain_->LastRealJ();

    int pos; // position in indices array, not really useful here
    int gidU, gid0;

    for (int j = J0; j <= J1; ++j)
        for (int i = I0; i <= I1; ++i)
        {
            // first element at grid point
            gidU = find_row0(nGlob_, mGlob_, i, j, 1);

            // reference to offset with 1-based identifier
            gid0 = gidU - 1;

            // H-equation: connections to H and Q points
            pos  = 0;
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_HH_);
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_QQ_);

            CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                           gid0 + SEAICE_HH_, pos, indices));

            // Q-equation: connections to Q and T points
            pos  = 0;
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_QQ_);
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_TT_);

            CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                           gid0 + SEAICE_QQ_, pos, indices));

            // M-equation: connections to H and M points
            pos  = 0;
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_HH_);
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_MM_);

            CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                           gid0 + SEAICE_MM_, pos, indices));

            // T-equation: connections to H, Q and T points
            pos = 0;
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_HH_);
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_QQ_);
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_TT_);

            CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                           gid0 + SEAICE_TT_, pos, indices));
        }


    // Last ordinary row in the grid, i.e., the last non-auxiliary row.
    int last = find_row0(nGlob_, mGlob_, nGlob_-1, mGlob_-1, SEAICE_NUN_);

    // global auxiliary row
    int auxRow = (aux_ == 1) ? last + aux_ : -1;
    
    // add entries for auxiliary condition
    if ( ( aux_ == 1) && standardMap_->MyGID(auxRow) )
    {
        // For the integral correction dependencies in this model
        // exist at nGlob_ * mGlob_ * 3 + aux_ points. Three of our
        // unknowns are integrated: Q,M and T.
        int len = nGlob_ * mGlob_ * 3 + aux_;
        int auxInds[len];
        int gid = 0;
        pos = 0;
        for (int j = 0; j < mGlob_; ++j)
            for (int i = 0; i < nGlob_; ++i)
                for (int xx = SEAICE_QQ_; xx <= SEAICE_TT_; ++xx)
                {
                    gid = find_row0(nGlob_, mGlob_, i, j, xx);
                    auxInds[pos] = gid;
                    pos++;
                }

        auxInds[pos++] = auxRow;

        assert(len == pos);
        CHECK_ZERO(matrixGraph_->InsertGlobalIndices(auxRow, len, auxInds));
    }

    // Finalize matrixgraph
    CHECK_ZERO(matrixGraph_->FillComplete() );
    
}

//=============================================================================
void SeaIce::initializePrec()
{
    Ifpack Factory;
    string precType = "Amesos"; // direct solve on subdomains with some overlap
    int overlapLevel = params_->get("Ifpack overlap level", 2);

    INFO("SeaIce: preconditioner overlap level: " << overlapLevel);

    // Create preconditioner
    precPtr_ = Teuchos::rcp(Factory.Create(precType, jac_.get(), overlapLevel));
    precPtr_->Initialize();
    precPtr_->Compute();
    precInitialized_ = true;
}

//=============================================================================
void SeaIce::solve(Teuchos::RCP<Epetra_MultiVector> const &b)
{
    // when using the preconditioner as a solver make sure the overlap
    // is large enough (depending on number of cores obv).
    applyPrecon(*b, *sol_);

    // compute residual
    Teuchos::RCP<Epetra_Vector> tmp = getSolution('C');
    applyMatrix(*sol_, *tmp);
    tmp->Update(1.0, *b, -1.0);
    tmp->Scale(1.0/Utils::norm(b));
    double normRes = Utils::norm(tmp);

    std::cout << *tmp << std::endl;
    
    INFO("SeaIce: solve ||b-Ax|| / ||b|| = " << normRes);
}

//=============================================================================
void SeaIce::applyMatrix(Epetra_MultiVector const &in,
                         Epetra_MultiVector &out)
{
    TIMER_START("SeaIce: apply matrix...");

    jac_->Apply(in, out);

    TIMER_STOP("SeaIce: apply matrix...");
}

//=============================================================================
void SeaIce::applyPrecon(Epetra_MultiVector const &in,
                         Epetra_MultiVector &out)
{
    TIMER_START("SeaIce: apply preconditioner...");
    if (!precInitialized_)
    {
        initializePrec();
    }
    if (recomputePrec_)
    {
        INFO("SeaIce: recomputing prec");
        precPtr_->Initialize();
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

    // INFO("SeaIce: preconditioner residual: " << rnorm);

    TIMER_STOP("SeaIce: apply preconditioner...");
}

//=============================================================================
void SeaIce::initializeState()
{
    state_->PutScalar(0.0);

    int maxit = 5;
    int it = 0;
    computeRHS();
    computeJacobian();
    Teuchos::RCP<Epetra_Vector> b;
    Teuchos::RCP<Epetra_Vector> t = getSolution('C');
    for (; it < maxit; ++it)
    {
        b = getRHS('C');
        b->Scale(-1.0);

        Teuchos::RCP<Epetra_MultiVector> bmu = Utils::Gather(*b, 0);       

        DUMP_VECTOR("b", *b    );
        DUMPMATLAB("A" , *jac_ );
        
        solve(b);
        DUMP_VECTOR("x", *sol_ );
        
        std::cout << "||b|| " << Utils::norm(b) << std::endl;
        std::cout << "||x|| " << Utils::norm(sol_) << std::endl;

        double normInf = jac_->NormInf();
        std::cout << "||A||inf " << normInf << std::endl;
        double normOne = jac_->NormOne();
        std::cout << "||A||one " << normOne << std::endl;
        double normFro = jac_->NormFrobenius();
        std::cout << "||A||frob " << normFro << std::endl;

        applyMatrix(*sol_, *t);
        std::cout << "||Ax||2 " << Utils::norm(t) << std::endl;
        
        getchar();
        
        state_->Update(1.0, *sol_, 1.0);

        computeRHS();
        computeJacobian();

        INFO( "SeaIce::initializeState() norm F = "
              << Utils::norm( getRHS('V') ) );


    }
}

//=============================================================================
void SeaIce::preProcess()
{
}

//=============================================================================
void SeaIce::postProcess()
{
    // save state -> hdf5
    if (saveState_)
        saveStateToFile(outputFile_); // Save to hdf5

}

//=============================================================================
std::string const SeaIce::writeData(bool describe)
{
    std::ostringstream datastring;

    if (describe)
    {
        datastring << std::setw(_FIELDWIDTH_)
                   << "SI vol";

        return datastring.str();
    }
    else
    {    
        datastring.precision(_PRECISION_);   

        Teuchos::RCP<Epetra_Vector> M = interfaceM();
        Teuchos::RCP<Epetra_Vector> H = interfaceH();

        Teuchos::RCP<Epetra_Vector> restr =
            Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

        for (int i = 0; i != M->MyLength(); ++i)
            (*restr)[i] = (*M)[i]*(*H)[i];
        
        double SIV = Utils::dot(intCoeff_, restr);
        
        datastring << std::scientific << std::setw(_FIELDWIDTH_)
                   << SIV;

        return datastring.str();
    }
}

//=============================================================================
int SeaIce::find_row0(int n, int m, int i, int j, int XX)
{
    //////////////////////////////
    // 0-based
    //////////////////////////////
    if (XX >= SEAICE_GG_)
    {
        if (aux_ <= 0)
        {
            ERROR("SeaIce: invalid row", __FILE__, __LINE__);
            return -1;
        }
        else
            return n * m * dof_ + (XX - SEAICE_GG_);
    }
    else // ordinary 0-based find_row
        return dof_ * ( j*n + i) + XX - 1;
}

//=============================================================================
int SeaIce::find_row1(int n, int m, int i, int j, int XX)
{
    //////////////////////////////
    // 1-based
    //////////////////////////////
    if (XX >= SEAICE_GG_)
    {
        if (aux_ <= 0)
        {
            ERROR("SeaIce: invalid row", __FILE__, __LINE__);
            return -1;
        }
        else
            return n * m * dof_ + (XX - SEAICE_GG_) + 1;
    }
    else // ordinary 1-based find_row
        return dof_ * ( (j-1)*n + (i-1)) + XX;
}

