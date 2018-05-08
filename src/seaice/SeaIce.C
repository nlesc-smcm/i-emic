#include "SeaIce.H"
#include "SeaIceDefinitions.H"

//=============================================================================
// Constructor
SeaIce::SeaIce(Teuchos::RCP<Epetra_Comm> comm, ParameterList params)
    :
    params_          (params),
    comm_            (comm),
    n_               (params->get("Global Grid-Size n", 180)),
    m_               (params->get("Global Grid-Size m", 90)),
    periodic_        (params->get("Periodic", false)),

    taus_         (0.01),   // threshold ice thickness
    epsilon_      (1e2),    // Heavyside approximation steepness
    
    // background mean values
    t0o_   (params->get("background ocean temp t0o", 7)),
    t0a_   (params->get("background atmos temp t0a", 10)),
    t0i_   (params->get("background seaice temp t0i",-15)),
    tvar_  (params->get("ocean temp variation tvar", 15)),
    s0_    (params->get("ocean background salinity s0", 35)),
    svar_  (params->get("ocean salinity variation svar", 1)),
    q0_    (params->get("atmos background humidity q0", 1e-3)),
    qvar_  (params->get("atmos humidity variation qvar", 5e-4)),
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

    // combined parameter
    zeta_  (ch_ * utau_ * rhoo_ * cpo_),

    // sublimation constants, parameters for saturation humidity over
    // ice
    c1_    (params->get("c1", 3.8e-3)),
    c2_    (params->get("c2", 21.87)),
    c3_    (params->get("c3", 265.5)),

    ce_    (params->get("Dalton number", 1.3e-03)),
    uw_    (params->get("mean atmospheric surface wind speed, ms^{-1}", 8.5)),

// typical vertical velocity
    eta_   ( ( rhoa_ / rhoo_ ) * ce_ * uw_),

// Shortwave radiation constants and functions
    alpha_   (params->get("albedo", 0.3)),
    sun0_    (params->get("solar constant", 1360)),
    c0_      (params->get("atmospheric absorption coefficient", 0.43)),
    Ch_      (params->get("Ch", 1.22e-3)),
    cpa_     (params->get("heat capacity", 1000)),

// exchange coefficient 
    muoa_    (rhoa_ * Ch_ * cpa_ * uw_)
    
{
    // Background sublimation and derivatives: calculate background
    // saturation specific humidity according to [Bolton,1980], T in
    // \deg C
    auto qsi = [&] (double t0i)
        {
            return c1_ * exp(c2_ * t0i / (t0i + c3_) );
        };

    auto dqsi = [&] (double t0i)
        {
            return (c1_ * c2_ * c3_) / pow(t0i + c3_, 2) * 
            exp( (c2_ * t0i) / (t0i + c3_) );
        };

    // Background sublimation and derivatives
    E0_    =  eta_ * ( qsi(t0i_) - q0_ );
    dEdT_  =  eta_ *  dqsi(t0i_);
    dEdq_  =  eta_ * -1;

    // background heat flux variation
    Qvar_  = zeta_;

    // background heat flux
    Q0_    = zeta_ * (freezingT(0) - t0o_) - rhoo_ * Lf_ * E0_;    
    
    xmin_ = params->get("Global Bound xmin", 286.0) * PI_ / 180.0;
    xmax_ = params->get("Global Bound xmax", 350.0) * PI_ / 180.0;
    ymin_ = params->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
    ymax_ = params->get("Global Bound ymax", 80.0)  * PI_ / 180.0;

    dof_  = SEAICE_NUN_;

    // Create domain object
    domain_ = Teuchos::rcp(new TRIOS::Domain(n_, m_, 1, dof_,
                                             xmin_, xmax_, ymin_, ymax_,
                                             periodic_, 1.0, comm_));

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

    // Local (overlapping) vectors
    localState_  = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localRHS_    = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localDiagB_  = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localSol_    = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));

    // External data (overlapping)
    localSST_    = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localSSS_    = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosT_ = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosQ_ = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));

    // Create local computational grid in x_, y_
    createGrid();
    
}

//=============================================================================
void SeaIce::computeRHS()
{
    // zero rhs vector
    localRHS_->PutScalar(0.0);

    // obtain view of rhs, state and external data
    double *rhs, *state, *sst, *sss, *tatm, *qatm;
    localRHS_->ExtractView(&rhs);
    localState_->ExtractView(&state);
    localSST_->ExtractView(&sst);
    localSSS_->ExtractView(&sss);
    localAtmosT_->ExtractView(&tatm);
    localAtmosQ_->ExtractView(&qatm);

    // row indices for rhs and data
    int rr, dr;
    double Tsi, Hval, Qval, Mval, val;
    for (int j = 0; j != mLoc_; ++j)
        for (int i = 0; i != nLoc_; ++i)
        {
            dr = j*nLoc_ + i;
            
            Hval = state[find_row(i, j, SEAICE_HH_)];
            Qval = state[find_row(i, j, SEAICE_QQ_)];
            Mval = state[find_row(i, j, SEAICE_MM_)];

            Tsi = iceSurfT(Qval, Hval, sss[dr]);
            
            for (int XX = 1; XX <= dof_; ++XX)
            {
                switch (XX)
                {

                case SEAICE_HH_:      // H row (thickness)

                    val = freezingT(sss[dr]) - sst[dr] - t0o_ - 
                        Q0_ / zeta_ - Qvar_ / zeta_ * Qval - 
                        ( rhoo_ * Lf_ / zeta_) * 
                        ( E0_ + dEdT_ * Tsi + dEdq_ * qatm[dr] );

                case SEAICE_QQ_:       // Q row (heat flux)
                    
                    val = 1 / muoa_ * (Q0_ + Qvar_ * Qval) - 
                        (sun0_ / 4 / muoa_) * shortwaveS(y_[j]) * (1-alpha_) * c0_ + 
                        (t0i_ + Tsi - tatm[dr] - t0a_) + 
                        (rhoo_ * Ls_ / muoa_) * 
                        (E0_ + dEdT_ * Tsi + dEdq_ * qatm[dr]);                  

                case SEAICE_MM_:       // M row (mask)
                    
                    val = Mval - 
                        (1/2) * (1 + tanh(epsilon_ * Hval) );
                    
                }
                
                rr = find_row(i, j, XX);
                rhs[rr] = val;
            }
        }
    domain_->Assembly2StandardSurface(*localRHS_, *rhs_);
}

//=============================================================================
// Build idealized local and global forcing vectors
void SeaIce::idealizedForcing()
{
    // extract views of sea surface temp sst
    //                  sea surface salinity sss
    //                  atmosphere temp tatm
    //                  atmosphere humidity qatm

    double *sst, *sss, *tatm, *qatm;
    localSST_->ExtractView(&sst);
    localSSS_->ExtractView(&sss);
    localAtmosT_->ExtractView(&tatm);
    localAtmosQ_->ExtractView(&qatm);

    int row;
    for (int j = 0; j != mLoc_; ++j)
        for (int i = 0; i != nLoc_; ++i)
        {
            row      = j * nLoc_ + i;

            sst[row]  = tvar_ * cos( PI_ * y_[j] / ymax_ );
            sss[row]  = svar_ * cos( PI_ * y_[j] / ymax_ ) / cos( y_[j] );
            tatm[row] = tvar_ * cos( PI_ * y_[j] / ymax_ );
            qatm[row] = qvar_ * cos( PI_ * y_[j] / ymax_ );
        }

    // Transfer data to non-overlapping vectors
    domain_->Assembly2StandardSurface(*localSST_,     *sst_);
    domain_->Assembly2StandardSurface(*localSSS_,     *sss_);
    domain_->Assembly2StandardSurface(*localAtmosT_, *tatm_);
    domain_->Assembly2StandardSurface(*localAtmosQ_, *qatm_);
}

//=====Loc========================================================================
// Create local grid points
void SeaIce::createGrid()
{
    // local dimensions
    xminLoc_ = domain_->XminLoc();
    xmaxLoc_ = domain_->XmaxLoc();
    yminLoc_ = domain_->YminLoc();
    ymaxLoc_ = domain_->YmaxLoc();
    
    // local grid dimensions
    nLoc_  =  domain_->LocalN();
    mLoc_  =  domain_->LocalM();

    dx_  =  (xmaxLoc_ - xminLoc_) / nLoc_;
    dy_  =  (ymaxLoc_ - yminLoc_) / mLoc_;

    for (int i = 0; i != nLoc_; ++i)
        x_.push_back(xminLoc_ + (i + 0.5) * dx_);

    for (int j = 0; j != mLoc_; ++j)
        y_.push_back(yminLoc_ + (j + 0.5) * dy_);
}
