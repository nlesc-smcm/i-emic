#include "SeaIce.H"
#include "SeaIceDefinitions.H"

//=============================================================================
// Constructor
SeaIce::SeaIce(Teuchos::RCP<Epetra_Comm> comm, ParameterList params)
    :
    params_          (params),
    comm_            (comm),
    nGlob_           (params->get("Global Grid-Size n", 16)),
    mGlob_           (params->get("Global Grid-Size m", 16)),
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
    dimLoc_ = mLoc_ * nLoc_ * dof_;

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

    // Construct local dependency grid:
    Al_ = std::make_shared<DependencyGrid>(nLoc_, mLoc_, 1, 1, dof_);

    createMatrixGraph();
    
    // Initialize Jacobian
    jac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *matrixGraph_));

}

//=============================================================================
void SeaIce::computeRHS()
{
    // zero rhs vector
    localRHS_->PutScalar(0.0);

    // obtain view of rhs, state and external data
    double *rhs, *state, *sst, *sss, *tatm, *qatm;
    localRHS_->ExtractView(&rhs);

    domain_->Standard2Assembly(*state_, *localState_);
    
    domain_->Standard2AssemblySurface(*sst_,   *localSST_);
    domain_->Standard2AssemblySurface(*sss_,   *localSSS_);
    domain_->Standard2AssemblySurface(*tatm_,  *localAtmosT_);
    domain_->Standard2AssemblySurface(*qatm_,  *localAtmosQ_);
    
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
            
            Hval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_HH_)];
            Qval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_QQ_)];
            Mval = state[find_row0(nLoc_, mLoc_, i, j, SEAICE_MM_)];

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

                    break;

                case SEAICE_QQ_:       // Q row (heat flux)
                    
                    val = 1. / muoa_ * (Q0_ + Qvar_ * Qval) - 
                        (sun0_ / 4. / muoa_) * shortwaveS(y_[j]) * (1.-alpha_) * c0_ + 
                        (t0i_ + Tsi - tatm[dr] - t0a_) + 
                        (rhoo_ * Ls_ / muoa_) * 
                        (E0_ + dEdT_ * Tsi + dEdq_ * qatm[dr]);

                    break; 

                case SEAICE_MM_:       // M row (mask)
                    
                    val = Mval - 
                        (1./2.) * (1. + tanh(epsilon_ * Hval) );

                    break;                    
                }
                
                rr = find_row0(nLoc_, mLoc_, i, j, XX);
                rhs[rr] = val;
            }
        }

    domain_->Assembly2Standard(*localRHS_, *rhs_);
}

//=============================================================================
void SeaIce::computeLocalJacobian()
{
    int HH = SEAICE_HH_;
    int QQ = SEAICE_QQ_;
    int MM = SEAICE_MM_;

    // reset entries
    Al_->zero();

    // obtain local state
    domain_->Standard2Assembly(*state_, *localState_);
    double *state;
    localState_->ExtractView(&state);

    // our range is the entire local domain (1-based)
    int range[8] = {1,nLoc_,1,mLoc_,1,1,1,1};

    // initialize dependencies
    double HH_HH, HH_QQ, QQ_HH, QQ_QQ, MM_MM;
    Atom MM_HH(nLoc_, mLoc_, 1, 1);

    // dHdt equation ----------------------
    HH_HH = -rhoo_ * Lf_ / zeta_ * dEdT_ * Q0_ / Ic_;
    Al_->set(range, HH, HH, HH_HH);

    HH_QQ = -Qvar_ / zeta_ - 
        rhoo_ * Lf_ / zeta_ * dEdT_ * H0_ * Qvar_ / Ic_;
    Al_->set(range, HH, QQ, HH_QQ);

    // Qtsa equation ----------------------
    QQ_HH = Q0_ / Ic_ + 
        rhoo_ * Ls_ / muoa_ * dEdT_ * Q0_ / Ic_;
    Al_->set(range, QQ, HH, QQ_HH);

    QQ_QQ = Qvar_ / muoa_ + 
        H0_ * Qvar_ / Ic_ + 
        rhoo_ * Ls_ / muoa_ * dEdT_ * H0_ * Qvar_ / Ic_;
    Al_->set(range, QQ, QQ, QQ_QQ);

    // Msi equation ----------------------
    // fill atom for nonlinear contribution H
    int ind;     // state index
    double val;        
    for (int j = 1; j <= mLoc_; ++j)
        for (int i = 1; i <= nLoc_; ++i)
        {
            ind  = find_row1(nLoc_, mLoc_, i, j, HH); // H row
            val  = -(epsilon_ / 2.0) * 
                ( 1.0 - pow(tanh(epsilon_ * (state[ind])), 2) );
            MM_HH.set( i, j, 1, 1, val);
        }
    Al_->set(range, MM, HH, MM_HH);

    MM_MM = 1.0;
    Al_->set(range, MM, MM, MM_MM);
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
            // H-equation
            gidU = find_row0(nGlob_, mGlob_, i, j, SEAICE_HH_);
            gid0 = gidU - 1;

            pos  = 0;
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_HH_);
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_QQ_);

            CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                           gid0 + SEAICE_HH_, pos, indices));

            // Q-equation
            pos  = 0;
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_HH_);
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_QQ_);

            CHECK_ZERO(matrixGraph_->InsertGlobalIndices(
                           gid0 + SEAICE_QQ_, pos, indices));

            // M-equation
            pos  = 0;
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_HH_);
            indices[pos++] = find_row0(nGlob_, mGlob_, i, j, SEAICE_MM_);

            int ierr = matrixGraph_->InsertGlobalIndices(
                gid0 + SEAICE_MM_, pos, indices);

            if (ierr != 0)
            {
                std::cout << "Error " << ierr
                          << " in matrixGraph_->InsertGlobalIndices " << std::endl;
                std::cout << "i: " << i << " j: " << j << std::endl;
                std::cout << "matrixGraph_->IndicesAreLocal() = "
                          << matrixGraph_->IndicesAreLocal() << std::endl;
                std::cout << "indices: ";
                for (int ii = 0; ii != pos; ++ii)
                {                
                    std::cout << indices[ii] << " ";
                }
                std::cout << std::endl;

                std::cout << "GRID = " << gid0 + SEAICE_MM_ << std::endl;
                std::cout << "LRID = " << matrixGraph_->LRID(gid0 + SEAICE_MM_);
                std::cout << std::endl;

                ERROR(" Error in InsertGlobalIndices ", __FILE__, __LINE__);
            }
            
        }
    
    // Finalize matrixgraph
    CHECK_ZERO(matrixGraph_->FillComplete() );
}
