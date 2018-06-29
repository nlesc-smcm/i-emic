#include "SeaIce.H"
#include "Ocean.H"
#include "Atmosphere.H"

//=============================================================================
// Constructor
SeaIce::SeaIce(Teuchos::RCP<Epetra_Comm> comm, ParameterList params)
    :
    params_          (params),
    comm_            (comm),
    nGlob_           (params->get("Global Grid-Size n", 16)),
    mGlob_           (params->get("Global Grid-Size m", 16)),
    periodic_        (params->get("Periodic", false)),

    precInitialized_ (false),
    recomputePrec_   (false),
    
    taus_         (0.01),   // threshold ice thickness
    epsilon_      (1e-2),    // Heavyside approximation steepness

    // background mean values
    t0o_   (params->get("background ocean temp t0o", 7)),
    t0a_   (params->get("background atmos temp t0a", 10)),
    t0i_   (params->get("background seaice temp t0i",-15)),
    s0_    (params->get("ocean background salinity s0", 35)),
    q0_    (params->get("atmos background humidity q0", 1e-3)),
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
    albe0_   (params->get("reference albedo", 0.3)),
    albed_   (params->get("albedo excursion", 0.5)),
    sun0_    (params->get("solar constant", 1360)),
    c0_      (params->get("atmospheric absorption coefficient", 0.43)),
    Ch_      (params->get("Ch", 1.22e-3)),
    cpa_     (params->get("heat capacity", 1000)),

// exchange coefficient
    muoa_    (rhoa_ * Ch_ * cpa_ * uw_)   

{
    // Continuation parameters
    allParameters_ = { "Combined Forcing",
                       "Solar Forcing",
                       "Mask Forcing" };

    parName_  = params->get( "Continuation parameter",
                             allParameters_[0] );

    comb_     = params->get(allParameters_[0], 0.0);
    sunp_     = params->get(allParameters_[1], 1.0);
    maskf_    = params->get(allParameters_[2], 1.0);               

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

    dof_      = SEAICE_NUN_;
    dimGlob_  = mGlob_ * nGlob_ * dof_;

    // Create domain object
    domain_ = Teuchos::rcp(new TRIOS::Domain(nGlob_, mGlob_, 1, dof_,
                                             xmin_, xmax_, ymin_, ymax_,
                                             periodic_, 1.0, comm_));

    // Compute 2D decomposition
    domain_->Decomp2D();

    // local dimensions
    xminLoc_ = domain_->XminLoc();
    xmaxLoc_ = domain_->XmaxLoc();
    yminLoc_ = domain_->YminLoc();
    ymaxLoc_ = domain_->YmaxLoc();

    // local grid dimensions
    nLoc_   =  domain_->LocalN();
    mLoc_   =  domain_->LocalM();
    dimLoc_ =  mLoc_ * nLoc_ * dof_;

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
    albe_       = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

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
    localAtmosA_ = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));

    // Create local computational grid in x_, y_
    createGrid();

    // Construct local dependency grid:
    Al_ = std::make_shared<DependencyGrid>(nLoc_, mLoc_, 1, 1, dof_);

    createMatrixGraph();

    // Initialize Jacobian
    jac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *matrixGraph_));

    // Initialize state
    initializeState();

    // Create importers for communication with other models
    int XX; 
    for (int i = 0; i != dof_; ++i)
    {
        XX = SEAICE_HH_ + i;
        Maps_[XX] = Utils::CreateSubMap(*standardMap_, dof_, XX);
        Imps_[XX] = Teuchos::rcp(new Epetra_Import(*Maps_[XX], *standardMap_));
    }
}

//=============================================================================
void SeaIce::computeRHS()
{
    // zero rhs vector
    localRHS_->PutScalar(0.0);

    // obtain view of rhs, state and external data
    double *rhs, *state, *sst, *sss, *tatm, *qatm, *albe;
    localRHS_->ExtractView(&rhs);

    domain_->Standard2Assembly(*state_, *localState_);

    domain_->Standard2AssemblySurface(*sst_,   *localSST_);
    domain_->Standard2AssemblySurface(*sss_,   *localSSS_);
    domain_->Standard2AssemblySurface(*tatm_,  *localAtmosT_);
    domain_->Standard2AssemblySurface(*qatm_,  *localAtmosQ_);
    domain_->Standard2AssemblySurface(*albe_,  *localAtmosA_);

    localState_->ExtractView(&state);
    
    localSST_->ExtractView(&sst);
    localSSS_->ExtractView(&sss);
    
    localAtmosT_->ExtractView(&tatm);
    localAtmosQ_->ExtractView(&qatm);
    localAtmosA_->ExtractView(&albe);

    // row indices for rhs and data
    int rr, sr;
    double Tsi, Hval, Qval, Mval, Tval, val;
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

                    val = comb_ * ( freezingT(sss[sr]) - sst[sr]
                                    - t0o_ - Q0_ / zeta_ )
                        - Qvar_ / zeta_ * Qval -
                        ( rhoo_ * Lf_ / zeta_) *
                        ( comb_ * E0_ + dEdT_ * Tsi + comb_ * dEdq_ * qatm[sr] );
                    
                    break;

                case SEAICE_QQ_:   // Q row (heat flux)

                    val = comb_ / muoa_ * Q0_ + Qvar_ / muoa_ * Qval -
                        (comb_ * sunp_ * sun0_ / 4. / muoa_) * shortwaveS(y_[j]) *
                        (1. - albe0_ - albed_*albe[sr]) * c0_ +
                        (Tval + comb_*(t0i_ - tatm[sr] - t0a_ )) +
                        (rhoo_ * Ls_ / muoa_) *
                        (comb_*E0_ + dEdT_ * Tval + comb_*dEdq_ * qatm[sr]);

                    break;

                case SEAICE_MM_:   // M row (mask)

                    val = Mval -
                        (comb_*maskf_/2.) * (1. + tanh( Hval / epsilon_ ) );

                    break;

                case SEAICE_TT_:   // T row (surface temperature)

                    val = comb_* ( freezingT(sss[sr]) - t0i_ +
                                   (Q0_*H0_ + H0_*Qvar_*Qval + Q0_*Hval) / Ic_ )
                        - Tval ;

                    break;
                }

                rr = find_row0(nLoc_, mLoc_, i, j, XX);
                rhs[rr] = val;
            }
        }

    domain_->Assembly2Standard(*localRHS_, *rhs_);
    INFO(" seaic F = " << Utils::norm(rhs_));
}

//=============================================================================
void SeaIce::computeLocalJacobian()
{
    int H = SEAICE_HH_;
    int Q = SEAICE_QQ_;
    int M = SEAICE_MM_;
    int T = SEAICE_TT_;

    // reset entries
    Al_->zero();

    // obtain local state
    domain_->Standard2Assembly(*state_, *localState_);
    double *state;
    localState_->ExtractView(&state);

    // our range is the entire local domain (1-based)
    int range[8] = {1, nLoc_, 1, mLoc_, 1, 1, 1, 1};

    // initialize dependencies
    double HH_HH, HH_QQ, QQ_QQ, QQ_TT, MM_MM, TT_HH, TT_QQ, TT_TT;
    Atom MM_HH(nLoc_, mLoc_, 1, 1);

    // dHdt equation ----------------------------------
    HH_HH = -rhoo_ * Lf_ / zeta_ * dEdT_ * Q0_ / Ic_;
    HH_QQ = -Qvar_ / zeta_ -
        rhoo_ * Lf_ / zeta_ * dEdT_ * H0_ * Qvar_ / Ic_;

    Al_->set(range, H, H, HH_HH);
    Al_->set(range, H, Q, HH_QQ);

    // Qtsa equation ---------------------------------
    QQ_QQ = Qvar_ / muoa_;
    QQ_TT = 1.0 + rhoo_ * Ls_ / muoa_ * dEdT_;

    Al_->set(range, Q, Q, QQ_QQ);
    Al_->set(range, Q, T, QQ_TT);

    // Msi equation -----------------------------------
    // fill atom for nonlinear contribution H
    int ind;     // state index
    double val;
    for (int j = 1; j <= mLoc_; ++j)
        for (int i = 1; i <= nLoc_; ++i)
        {
            ind  = find_row1(nLoc_, mLoc_, i, j, H); // H row
            val  = -(1 / (epsilon_ * 2.0) ) *
                ( 1.0 - pow(tanh((state[ind]) / epsilon_), 2) );
            MM_HH.set( i, j, 1, 1, val);
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

                for (int B = 1; B <= dof_; ++B)
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

// ---------------------------------------------------------------------------
// Set specific continuation parameter
void SeaIce::setPar(std::string const &parName, double value)
{
    parName_ = parName; // Overwrite our parameter name

    if (parName.compare(allParameters_[0]) == 0)
        comb_ = value;
    else if (parName.compare(allParameters_[1]) == 0)
        sunp_ = value;
    else if (parName.compare(allParameters_[2]) == 0)
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
    int A = ATMOS_AA_; // (1-based) in the Atmosphere, humidity is the second unknown

    // compute a few constant derivatives (see computeRHS)

    // d / dq_atm (F_H)
    double dqatmFH = -(rhoo_ * Lf_ / zeta_) * dEdq_;

    // d / dt_atm (F_Q)
    double dtatmFQ = -1;

    // d / dq_atm (F_Q)
    double dqatmFQ =  (rhoo_ * Ls_ / muoa_) * dEdq_;

    // d / da_atm (F_Q)
    double daatmFQ;
    double tmp = 0.0;

    for (int j = 0; j != mGlob_; ++j)
    {
        int gid = j * nGlob_;                      
        int lid = standardSurfaceMap_->LID(gid);   
        
        if (lid >= 0)
            tmp = (sun0_ / 4. / muoa_) *
                shortwaveS(y_[lid / nLoc_]) * albed_ * c0_;

        comm_->SumAll(&tmp, &daatmFQ, 1);
        
        for (int i = 0; i != nGlob_; ++i)
            for (int XX = 1; XX <= dof_; ++XX)
            {
                block->beg.push_back(el_ctr);

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
    double dTFH = -1;

    // d / dS (F_H)
    double dSFH = a0_ - ( rhoo_ * Lf_ / zeta_) * dEdT_ * a0_;

    // d / dS (F_T)
    double dSFT = a0_;

    for (int j = 0; j != mGlob_; ++j)
        for (int i = 0; i != nGlob_; ++i)
            for (int XX = 1; XX <= dof_; ++XX)
            {
                block->beg.push_back(el_ctr);

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
    return out;
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
    double qvar = 5e-4;
    
    double *sst, *sss, *tatm, *qatm, *albe;
    localSST_->ExtractView(&sst);
    localSSS_->ExtractView(&sss);
    localAtmosT_->ExtractView(&tatm);
    localAtmosQ_->ExtractView(&qatm);
    localAtmosA_->ExtractView(&albe);

    int row;
    for (int j = 0; j != mLoc_; ++j)
        for (int i = 0; i != nLoc_; ++i)
        {
            row      = j * nLoc_ + i;

            sst[row]  = tvar * cos( PI_ * y_[j] / ymax_ );
            sss[row]  = svar * cos( PI_ * y_[j] / ymax_ ) / cos( y_[j] );
            tatm[row] = tvar * cos( PI_ * y_[j] / ymax_ );
            qatm[row] = qvar * cos( PI_ * y_[j] / ymax_ );
            albe[row] = 0.0;
        }

    // Transfer data to non-overlapping vectors
    domain_->Assembly2StandardSurface(*localSST_,     *sst_);
    domain_->Assembly2StandardSurface(*localSSS_,     *sss_);
    domain_->Assembly2StandardSurface(*localAtmosT_, *tatm_);
    domain_->Assembly2StandardSurface(*localAtmosQ_, *qatm_);
    domain_->Assembly2StandardSurface(*localAtmosA_, *albe_);
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
// Create matrix graph
void SeaIce::createMatrixGraph()
{
    // We do not have any spatial relation between the unknowns so
    // maxDeps is at most the number of unknowns.
    int maxDeps = dof_;

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
            gidU = find_row0(nGlob_, mGlob_, i, j, SEAICE_HH_);
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

    // Finalize matrixgraph
    CHECK_ZERO(matrixGraph_->FillComplete() );
}

//=============================================================================
void SeaIce::initializePrec()
{
    Ifpack Factory;
    string precType = "Amesos"; // direct solve on subdomains with some overlap
    int overlapLevel = params_->get("Ifpack overlap level", 0);

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
    sol_->PutScalar(0.0);

    int maxit = 5;
    int it = 0;
    computeRHS();
    for (; it < maxit; ++it)
    {
        computeJacobian();
        rhs_->Scale(-1.0);
        solve(rhs_);
        state_->Update(1.0, *sol_, 1.0);
        computeRHS();
        INFO("SeaIce::initializeState() norm F = " << Utils::norm(rhs_));
    }
}


//=============================================================================
void SeaIce::preProcess()
{}

//=============================================================================
void SeaIce::postProcess()
{}

