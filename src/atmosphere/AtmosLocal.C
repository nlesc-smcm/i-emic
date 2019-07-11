#include "AtmosLocal.H"
#include "AtmosphereDefinitions.H"
#include "DependencyGrid.H"

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm> // std::fill in assemble

#include "GlobalDefinitions.H"
#include "my_f2c.H"

extern "C" _SUBROUTINE_(getdeps)(double*, double*, double*,
                                 double*, double*, double*, double *);

//==================================================================
// Constructor for use with parallel atmosphere
AtmosLocal::AtmosLocal(int n, int m, int l, bool periodic,
                       double xmin, double xmax, double ymin, double ymax,
                       Teuchos::RCP<Teuchos::ParameterList> params)
    :
    params_   (params),

// grid --------------------------------------------------------------------
    n_               (n),
    m_               (m),
    l_               (l),

    aux_             (params->get("Auxiliary unknowns", 1)),

    xmin_            (xmin),
    xmax_            (xmax),
    ymin_            (ymin),
    ymax_            (ymax),

    periodic_        (periodic)
{
    INFO("AtmosLocal: constructor for parallel use...");

    // Define domain
    xmin_glob_ = params->get("Global Bound xmin", 286.0) * PI_ / 180.0;
    xmax_glob_ = params->get("Global Bound xmax", 350.0) * PI_ / 180.0;
    ymin_glob_ = params->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
    ymax_glob_ = params->get("Global Bound ymax", 74.0)  * PI_ / 180.0;

    // Set non-grid (domain decomposition) related parameters
    setParameters(params);

    // Notify that we are now in parallel mode
    parallel_ = true;

    // Other setup things
    setup();

    INFO("AtmosLocal: constructor for parallel use done");
}

//==================================================================
// Constructor for standalone serial use
AtmosLocal::AtmosLocal(Teuchos::RCP<Teuchos::ParameterList> params)
    :
    params_          (params),

// grid ------------------------------------------------------------------
    n_               (params->get("Global Grid-Size n", 16)),
    m_               (params->get("Global Grid-Size m", 16)),
    l_               (params->get("Global Grid-Size l", 1)),
    aux_             (0),
    periodic_        (params->get("Periodic", false))
{
    INFO("AtmosLocal: constructor...");

    // Set non-grid (domain decomposition) related parameters
    setParameters(params);

    // Define domain
    xmin_glob_ = params->get("Global Bound xmin", 286.0) * PI_ / 180.0;
    xmax_glob_ = params->get("Global Bound xmax", 350.0) * PI_ / 180.0;
    ymin_glob_ = params->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
    ymax_glob_ = params->get("Global Bound ymax", 74.0)  * PI_ / 180.0;

    // Global domain is the local domain in the serial case
    xmin_ = xmin_glob_;
    xmax_ = xmax_glob_;
    ymin_ = ymin_glob_;
    ymax_ = ymax_glob_;

    parallel_ = false; // this is not the constructor for parallel use.

    // More factorized setup stuff
    setup();

    // Set row for integral condition (only in serial case)
    // use the final q-row (1-based!!)
    rowIntCon_ = find_row(n_, m_, l_, ATMOS_QQ_);

    // Setup coefficients for integrals
    setupIntCoeff();

    INFO("AtmosLocal: constructor... done");
}

//==================================================================
void AtmosLocal::setParameters(Teuchos::RCP<Teuchos::ParameterList> params)
{
    // physics ------------------------------------------------------------------
    rhoa_            = params->get("atmospheric density",1.25);
    rhoo_            = params->get("oceanic density",1024);
    hdima_           = params->get("atmospheric scale height",8400.);
    hdimq_           = params->get("humidity scale height",1800.);
    hdim_            = params->get("vertical length scale",4000.);
    cpa_             = params->get("heat capacity",1000.);
    D0_              = params->get("temperature eddy diffusivity",3.1e+06);
    kappa_           = params->get("humidity eddy diffusivity",1e+06);
    arad_            = params->get("radiative flux param A",212.0);
    brad_            = params->get("radiative flux param B",1.5);
    sun0_            = params->get("solar constant",1360.);
    c0_              = params->get("atmospheric absorption coefficient",0.43);
    ce_              = params->get("Dalton number",1.3e-03);
    ch_              = params->get("exchange coefficient ch",0.94 * ce_);
    uw_              = params->get("mean atmospheric surface wind speed",8.5);
    t0a_             = params->get("background temperature atmosphere",15.0); //(C)
    t0o_             = params->get("background temperature ocean",15.0);      //(C)
    t0i_             = params->get("background temperature seaice",-5.0);      //(C)
    tdim_            = params->get("temperature scale", 1.0); // ( not used)
    q0_              = params->get("atmos reference humidity",2e-3); // (kg/kg)
    qdim_            = params->get("atmos humidity scale", 1e-3);  // (kg/kg)
    lv_              = params->get("latent heat of vaporization", 2.5e06); // (J/kg)

    udim_            = params->get("horizontal velocity of the ocean", 0.1e+00);
    r0dim_           = params->get("radius of the earth", 6.37e+06);

    // Albedo parameters
    a0_              = params->get("reference albedo", 0.3);
    da_              = params->get("albedo excursion", 0.5);

    tauf_            = params->get("restoring timescale tauf (in days)", 1.0);
    tauc_            = params->get("restoring timescale tauc (in days)", 1.0);

    // restoring timescales in model time
    tauf_            = (tauf_ * 3600. * 24. * udim_) / r0dim_;
    tauc_            = (tauc_ * 3600. * 24. * udim_) / r0dim_;

    Tm_              = params->get("melt temperature threshold (deg C)", 0.0);
    Tr_              = params->get("rain/snow temperature threshold (deg C)", 1.0);
    Pa_              = params->get("accumulation precipitation threshold (m/y)", 0.2);
    epm_             = params->get("melt threshold width (deg C)", 5.0);
    epr_             = params->get("rain/snow threshold width (deg C)", 1.0);
    epa_             = params->get("accumulation threshold width (m/y)", 0.1);

    // continuation ----------------------------------------------------------------
    allParameters_   = { "Combined Forcing",
                         "Solar Forcing",
                         "Longwave Forcing",
                         "Humidity Forcing",
                         "Latent Heat Forcing",
                         "Albedo Forcing",
                         "T Eddy Diffusivity"};

    // starting values
    int ctr = 0;
    comb_            = params->get(allParameters_[ctr++], 0.0);
    sunp_            = params->get(allParameters_[ctr++], 1.0);
    lonf_            = params->get(allParameters_[ctr++], 1.0);
    humf_            = params->get(allParameters_[ctr++], 1.0);
    latf_            = params->get(allParameters_[ctr++], 1.0);
    albf_            = params->get(allParameters_[ctr++], 1.0);
    tdif_            = params->get(allParameters_[ctr++], 1.0);
}

//==================================================================
void AtmosLocal::setup()
{
    // Filling the coefficients
    muoa_ =  rhoa_ * ch_ * cpa_ * uw_;
    amua_ = (arad_ + brad_ * t0a_) / muoa_;
    bmua_ =  brad_ / muoa_;
    Ai_   =  rhoa_ * hdima_ * cpa_ * udim_ / (r0dim_ * muoa_);
    Ad_   =  rhoa_ * hdima_ * cpa_ * D0_ / (muoa_ * r0dim_ * r0dim_);
    As_   =  sun0_ * (1 - c0_) / (4 * muoa_);

    // Filling more coefficients (humidity specific)
    eta_ =  (rhoa_ / rhoo_) * ce_ * uw_;

    // nuq will be our continuation parameter for the humidity related
    // physics and is adjusted in setpar.
    nuq_ = comb_* humf_ * (eta_ / hdimq_ ) * (rhoo_ / rhoa_) * (r0dim_ / udim_);

    Phv_ =  kappa_ / (udim_ * r0dim_);

    // Parameters for saturation humidity over ocean and ice
    double c1 = 3.8e-3; // (kg / kg)
    double c2 = 21.87;  //
    double c3 = 265.5;  // (K)
    double c4 = 17.67;  //
    double c5 = 243.5;  // (K)

    // background ice temperature is chosen such that background
    // evaporation and sublimation cancel.
    // t0i_ = c3*c4*t0o_ / (c2*c5+(c2-c4)*t0o_);

    // Calculate background saturation specific humidity according to
    // [Bolton,1980], T in \deg C
    qso_   = c1 * exp( c4 * t0o_ / (t0o_ + c5) );
    qsi_   = c1 * exp( c2 * t0i_ / (t0i_ + c3) );

    // Dimensional background evaporation and precipitation
    Eo0_ = eta_ * ( qso_ - q0_ );
    Ei0_ = eta_ * ( qsi_ - q0_ );

    // Background correction factor for sublimation
    Cs_ = (Ei0_ - Eo0_) / eta_ / qdim_;

    // Backgr. precipitation is taken equal to backgr. evaporation over ocean
    Po0_ = Eo0_;

    // Convert threshold precipitation to threshold P deviation
    // Pa_  = Pa_ / 3600. / 24. / 365.;    // convert to m/s (deprecated)
    // Pa_  = (Pa_ - Po0_) / eta_ / qdim_; // convert to deviation (deprecated)
    // epa_ = epa_ / 3600. / 24. / 365.;   // convert threshold width to m/s
    // epa_ = epa_ / eta_ / qdim_;         // convert to deviation (deprecated)

    // Convert threshold temperatures to T deviations
    Tr_ = Tr_ - t0o_;
    Tm_ = Tm_ - t0o_;

    // Calculate saturation humidity derivatives at ref. temps
    dqso_  = (c1 * c4 * c5) / pow(t0o_ + c5, 2);
    dqso_ *= exp( (c4 * t0o_) / (t0o_ + c5) );

    dqso_  = 5e-4; // hack

    dqsi_  = (c1 * c2 * c3) / pow(t0i_ + c3, 2);
    dqsi_ *= exp( (c2 * t0i_) / (t0i_ + c3) );

    // large temperature response in E
    double EdevT = eta_ * qdim_ * ( tdim_ / qdim_ * dqso_ * -17 );

    // large humidity response in E
    double Edevq = eta_ * qdim_ * ( -5.0 );

    // latent heat due to precipitation coeff
    lvscale_ = rhoo_ * lv_ / muoa_ ;

    // Get ocean parameters
    double tmp1, tmp2, tmp3, tmp4, tmp5;
    FNAME(getdeps)(&Ooa_, &Os_, &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);

    INFO("AtmosLocal computed parameters: ");
    INFO("       mu   = " << muoa_);
    INFO("   B / mu   = " << bmua_);
    INFO("       Ad   = " << Ad_);
    INFO("      Phv   = " << Phv_);
    INFO("      Ooa   = " << Ooa_);
    INFO("      nuq   = " << nuq_);
    INFO("      eta   = " << eta_);
    INFO("      qso   = " << qso_);
    INFO("      qsi   = " << qsi_);
    INFO("     dqso   = " << dqso_);
    INFO("     dqsi   = " << dqsi_);
    INFO("     qdim   = " << qdim_);
    INFO("   A*DpDq   = " << -eta_ * nuq_);
    INFO("    DqDt0   = " << nuq_ * tdim_ / qdim_ * dqso_);
    INFO("    DqDti   = " << nuq_ * tdim_ / qdim_ * dqsi_);
    INFO("  lvscale   = " << lvscale_);
    INFO("      Eo0   = " << Eo0_ << " m/s");
    INFO("      Ei0   = " << Ei0_ << " m/s");
    INFO("       Cs   = " << Cs_);
    INFO("    EdevT   = " << EdevT);
    INFO("    Edevq   = " << Edevq);
    INFO("      Po0   = " << Po0_ << " m/s = "
         << Po0_ * 3600 * 24 * 365 << " m/y (background precipitation)" );
    INFO("     QLH0   = " << Po0_ * rhoo_ * lv_ << " Wm^{-2} (background latent heat flux)");

    INFO(std::endl << "AtmosLocal all xml parameters: ");
    INFO(*params_);
    INFO(std::endl);

    np_  = ATMOS_NP_;   // all neighbouring points including the center
    nun_ = ATMOS_NUN_;  // ATMOS_TT_, ATMOS_QQ_, ATMOS_AA_ (ATMOS_PP_
                        // exists at only aux_ points).

    // Problem size
    dim_ = m_ * n_ * l_ * nun_ + aux_;

    // Initialize state, rhs and solution of linear solve with zeros
    rhs_   = std::make_shared<std::vector<double> >(dim_, 0.0);
    sol_   = std::make_shared<std::vector<double> >(dim_, 0.0);
    state_ = std::make_shared<std::vector<double> >(dim_, 0.0);
    diagB_ = std::make_shared<std::vector<double> >(dim_, 0.0);

    // Initialize fields for evaporation and precipitation at surface
    E_   = std::make_shared<std::vector<double> >(m_ * n_, 0.0);

    // Precipitation
    P_   = std::make_shared<std::vector<double> >(m_ * n_, 0.0);

    // Precipitation distribution function
    Pdist_ = std::vector<double>(m_ * n_, 1.0);

    // Initialize land surface temperature
    lst_ = std::make_shared<std::vector<double> >(m_ * n_, 0.0);

    // Initialize ocean surface temperature
    sst_ = std::make_shared<std::vector<double> >(m_ * n_, 0.0);

    // Initialize sea ice surface temperature
    sit_ = std::make_shared<std::vector<double> >(m_ * n_, 0.0);

    // Sea ice mask
    Msi_ = std::make_shared<std::vector<double> >(m_ * n_, 0.0);

    // Initialize surface mask
    surfmask_ = std::make_shared<std::vector<int> >(m_ * n_, 0);

    // Initialize forcing with zeros
    frc_ = std::vector<double>(dim_, 0.0);

    if (!parallel_)
    {
        // Number of super and sub-diagonal bands in banded matrix
        ksub_ = nun_ * std::max(n_, m_);
        ksup_ = nun_ * std::max(n_, m_);

        // Leading dimension of banded matrix
        ldimA_  = 2 * ksub_ + 1 + ksup_;

        // Initialize banded storage
        bandedA_ = std::vector<double>(ldimA_  * dim_, 0.0);
        buildLU_ = true;
        // Create pivot array for use in lapack
        ipiv_ = std::vector<int> (dim_+1, 0);
    }

    // Construct dependency grid:
    Al_ = std::make_shared<DependencyGrid>(n_, m_, l_, np_, nun_ + aux_);

    // Set the grid increments
    dx_ = (xmax_ - xmin_) / n_;
    dy_ = (ymax_ - ymin_) / m_;

    // Fill x
    xu_.reserve(n_+1);
    xc_.reserve(n_+1);
    for (int i = 0; i != n_+1; ++i)
    {
        xu_.push_back(xmin_ + i * dx_);
        xc_.push_back(xmin_ + (i - 0.5) * dx_);
    }

    // Fill y and latitude-based arrays
    yv_.reserve(m_+1);
    yc_.reserve(m_+1);

    datc_.reserve(m_+1);
    datv_.reserve(m_+1);
    suna_.reserve(m_+1);
    suno_.reserve(m_+1);

    for (int j = 0; j != m_+1; ++j)
    {
        yv_.push_back( ymin_ + j * dy_ );
        yc_.push_back( ymin_ + ( j - 0.5) * dy_ );

        datc_.push_back(0.9 + 1.5 * exp(-12 * yc_[j] * yc_[j] / PI_));
        datv_.push_back(0.9 + 1.5 * exp(-12 * yv_[j] * yv_[j] / PI_));
        suna_.push_back(As_*(1 - .482 * (3 * pow(sin(yc_[j]), 2) - 1.) / 2.));
        suno_.push_back(Os_*(1 - .482 * (3 * pow(sin(yc_[j]), 2) - 1.) / 2.));
    }

    if (periodic_)
    {
        INFO("Periodicity enabled, which should only be the case in a serial run.");
    }
}

//-----------------------------------------------------------------------------
// Destructor
AtmosLocal::~AtmosLocal()
{
    INFO("AtmosLocal destructor");
}

//------------------------------------------------------------------
void AtmosLocal::setupIntCoeff()
{
    // Create serial integration coefficients for integral condition on q
    intcondCoeff_ = std::make_shared<std::vector<double> >(dim_, 0.0);

    // obtain indices and values for integration
    std::vector<double> vals, inds;
    integralCoeff(vals, inds);

    // Fill coefficients
    for (size_t idx = 0; idx != inds.size(); ++idx)
    {
        (*intcondCoeff_)[inds[idx]-1] = vals[idx];
    }

    // Create serial integration coefficients for precipitation integral
    // Use 1 dof and ignore land points
    pIntCoeff_ = std::make_shared<std::vector<double> >(m_ * n_, 0.0);
    integralCoeff(vals, inds, 1);

    // test indices
    assert(inds.back()-1 < pIntCoeff_->size());

    // Fill coefficients
    for (size_t idx = 0; idx != inds.size(); ++idx)
    {
        (*pIntCoeff_)[inds[idx]-1] = vals[idx];
    }

    // Total surface area (ignoring land)
    totalArea_ = Utils::sum(*pIntCoeff_);

    INFO("AtmosLocal: total E,P area = " << totalArea_);
}

//------------------------------------------------------------------
void AtmosLocal::idealized(double precip)
{
    // put idealized values in atmosphere
    double value;
    int rowSST, rowTT, rowQQ, rowPP, rowAA;
    for (int i = 1; i <= n_; ++i)
        for (int j = 1; j <= m_; ++j)
        {
            value = cos(PI_*(yc_[j]-ymin_glob_)/(ymax_glob_-ymin_glob_));

            rowSST = n_*(j-1) + (i-1);
            rowTT  = find_row(i,j,l_,ATMOS_TT_)-1;
            rowQQ  = find_row(i,j,l_,ATMOS_QQ_)-1;

            rowAA  = find_row(i,j,l_,ATMOS_AA_)-1;

            // These values are chosen such that the integrated
            // evaporation is zero.
            (*sst_) [rowSST] = value;
            (*state_)[rowTT] = value;
            (*state_)[rowQQ] = value * tdim_ * dqso_ / qdim_;
            (*state_)[rowAA] = a0_;
        }

    // Compute evaporation based on idealized sst and q
    computeEvaporation();

    // Set idealized precipitation
    if (aux_ == 1)
    {
        rowPP = find_row(n_ ,m_, l_, ATMOS_PP_) - 1;
        (*state_)[rowPP] = precip;
    }
}

//-----------------------------------------------------------------------------
void AtmosLocal::zeroState()
{
    // Set state to zero
    int dim = n_ * m_ * l_;
    *state_ = std::vector<double>(dim, 0.0);
}

//-----------------------------------------------------------------------------
void AtmosLocal::zeroOcean()
{
    // Set sst to zero
    *sst_ = std::vector<double>(m_ * n_, 0.0);
}

//-----------------------------------------------------------------------------
void AtmosLocal::setOceanTemperature(std::vector<double> const &sst)
{
    assert((int) sst.size() == n_ * m_);
    *sst_ = sst;
}

//-----------------------------------------------------------------------------
void AtmosLocal::setSeaIceTemperature(std::vector<double> const &sit)
{
    assert((int) sit.size() == n_ * m_);
    *sit_ = sit;
}

//-----------------------------------------------------------------------------
void AtmosLocal::setSeaIceMask(std::vector<double> const &Msi)
{
    assert((int) Msi.size() == n_ * m_);
    *Msi_ = Msi;
}

//-----------------------------------------------------------------------------
void AtmosLocal::fillPdist(double *Pdist)
{
    int sr, pos = 0;
    bool on_land;
    double y;
    for (int j = 1; j <= m_; ++j)
        for (int i = 1; i <= n_; ++i)
        {
            sr = n_*(j-1) + (i-1); // surface row
            on_land = (*surfmask_)[sr];
            if (!on_land)
            {
                y = yc_[j];
                Pdist[pos] = 2*exp(-pow(6*y,2))+pow(sin(2.0*y),2);
            }
            else
            {
                Pdist[pos] = 0.0;
            }
            pos++;
        }
}

//-----------------------------------------------------------------------------
void AtmosLocal::setPdist(double *Pdist)
{
    int sr, pos = 0;
    bool on_land;
    for (int j = 1; j <= m_; ++j)
        for (int i = 1; i <= n_; ++i)
        {
            sr = n_*(j-1) + (i-1); // surface row
            on_land = (*surfmask_)[sr];
            if (on_land)
                Pdist_[pos] = 0.0;
            else
                Pdist_[pos] = Pdist[pos];
            pos++;
        }
}

//-----------------------------------------------------------------------------
void AtmosLocal::getCommPars(AtmosLocal::CommPars &parStruct) const
{
    parStruct.tdim = tdim_;
    parStruct.qdim = qdim_;
    parStruct.nuq  = nuq_;
    parStruct.eta  = eta_;
    parStruct.dqso = dqso_;
    parStruct.dqsi = dqsi_;
    parStruct.dqdt = nuq_ * tdim_ / qdim_ * dqso_ ;
    parStruct.Eo0  = Eo0_;
    parStruct.Ei0  = Ei0_;
    parStruct.Cs   = Cs_;
    parStruct.t0o  = t0o_;
    parStruct.t0i  = t0i_;
    parStruct.a0   = a0_;
    parStruct.da   = da_;
    parStruct.tauf = tauf_;
    parStruct.tauc = tauc_;
    parStruct.comb = comb_;
    parStruct.albf = albf_;
}

//-----------------------------------------------------------------------------
void AtmosLocal::integralCoeff(std::vector<double> &val,
                               std::vector<double> &ind, int nun, bool ignoreLand)
{
    // Clear arrays
    val.clear();
    ind.clear();
    // Assuming that nun > 1 implies we want the integral coefficients
    // at the QQ points
    int XX = (nun == 1) ? ATMOS_TT_ : ATMOS_QQ_;

    // Obtain values and indices to compute integral
    // 1-based!
    for (int k = 1; k <= l_; ++k)
        for (int j = 1; j <= m_; ++j)
            for (int i = 1; i <= n_; ++i)
            {
                if (ignoreLand && (*surfmask_)[(j-1)*n_+(i-1)])
                    continue;

                val.push_back(cos(yc_[j]) * dx_ * dy_);
                ind.push_back(FIND_ROW_ATMOS1(nun, n_, m_, l_, i, j, k, XX));
            }
}

//-----------------------------------------------------------------------------
void AtmosLocal::computeJacobian()
{
    TIMER_START("AtmosLocal: compute Jacobian...");

    // Reset dependency grid
    Al_->zero();

    //------------------------------------------------------------------
    // Temperature equation dependencies
    Atom tc (n_, m_, l_, np_);
    Atom tc2(n_, m_, l_, np_);
    Atom txx(n_, m_, l_, np_);
    Atom tyy(n_, m_, l_, np_);

    discretize(1, tc);  // sensible heat flux component
    discretize(2, tc2); // outgoing longwave radiation component
    discretize(3, txx); // longitudinal diffusion
    discretize(4, tyy); // meridional diffusion

    // Combine atoms:
    // Al(:,:,:,:,TT,TT) = Ad * (txx + tyy) - tc - bmua*tc2
    txx.update(tdif_*Ad_, tdif_*Ad_, tyy, -1.0, tc, -bmua_, tc2);

    // Set temperature atom in dependency grid
    Al_->set({1,n_,1,m_,1,l_,1,np_}, ATMOS_TT_, ATMOS_TT_, txx);

    if (aux_ == 1) // latent heat due to precipitation
    {
        Atom TT_PP(n_, m_, l_, np_);

        // Central P dependency can be implemented in this way. Note
        // that we can add P dependencies in non-P rows, but it is not
        // possible to use the stencil approach to add dependencies to
        // a P-row, as there is only a single row corresponding to,
        // e.g., ATMOS_PP_. The nested loop in assemble would give
        // some trouble. The final P-row is implemented from the
        // parallel side. (Similar to the integral condition.)

        int sr;
        double value;
        for (int j = 1; j <= m_; ++j)
            for (int i = 1; i <= n_; ++i)
            {
                sr = n_*(j-1) + (i-1); // surface row
                value = comb_ * latf_ * lvscale_ *
                    eta_ * qdim_ * Pdist_[sr];
                TT_PP.set(i,j,l_, 5, value);
            }
        Al_->set({1,n_,1,m_,1,l_,1,np_}, ATMOS_TT_, ATMOS_PP_, TT_PP);
    }

    // Temperature equation: albedo dependence TT_AA
    Atom TT_AA(n_, m_, l_, np_);
    bool on_land;
    int sr, pr;
    double dTdA;
    for (int j = 1; j <= m_; ++j)
        for (int i = 1; i <= n_; ++i)
        {
            sr = n_*(j-1) + (i-1); // surface row
            on_land = (*surfmask_)[sr];
            if (on_land)
                dTdA =   dTldA(j) + dTadA(j);
            else // above ocean / sea ice
                dTdA =   dTadA(j);

            TT_AA.set(i, j, l_, 5, dTdA);
        }

    Al_->set({1,n_,1,m_,1,l_,1,np_}, ATMOS_TT_, ATMOS_AA_, TT_AA);

    //------------------------------------------------------------------
    // Humidity equation
    Atom qc (n_, m_, l_, np_);
    Atom qxx(n_, m_, l_, np_);
    Atom qyy(n_, m_, l_, np_);
    Atom QQ_PP(n_, m_, l_, np_);

    discretize(1, qc );  // Evaporation component (over land no evaporation)
    discretize(5, qxx);  // longitudinal diffusion
    discretize(6, qyy);  // meridional diffusion

    // Combine atoms
    qxx.update(Phv_, Phv_, qyy, -nuq_, qc);

    // Set humidity atom in dependency grid
    // Al(:,:,:,:,QQ,QQ) = Phv * (qxx + qyy) - nuq * qc
    Al_->set( {1,n_,1,m_,1,l_,1,np_}, ATMOS_QQ_, ATMOS_QQ_, qxx);
    if (aux_ == 1)
    {
        // Al(:,:,:,:,QQ, PP) = - nuq * p(theta) * P, :
        int sr;
        double value;
        for (int j = 1; j <= m_; ++j)
            for (int i = 1; i <= n_; ++i)
            {
                sr = n_*(j-1) + (i-1); // surface row
                value = -nuq_ * Pdist_[sr];
                QQ_PP.set(i,j,l_, 5, value);
            }

        Al_->set( {1,n_,1,m_,1,l_,1,np_}, ATMOS_QQ_, ATMOS_PP_, QQ_PP);
    }

    //------------------------------------------------------------------
    // Albedo equation, AA_AA dependency
    Atom AA_AA(n_, m_, l_, np_);
    Atom AA_TT(n_, m_, l_, np_);
    Atom AA_PP(n_, m_, l_, np_);

    double dAdA, dAdT, dAdP;
    double A, Ta, P = 0;
    int ar,tr;
    for (int j = 1; j <= m_; ++j)
        for (int i = 1; i <= n_; ++i)
        {
            sr = n_*(j-1) + (i-1); // surface row
            ar = find_row(i, j, l_, ATMOS_AA_) - 1;
            tr = find_row(i, j, l_, ATMOS_TT_) - 1;
            if (aux_ == 1)
                pr = find_row(i, j, l_, ATMOS_PP_) - 1; // precipitation row
            else
                pr = -1;

            on_land = (*surfmask_)[sr];
            A  = (*state_)[ar];
            Ta = (*state_)[tr];
            if (pr >= 0)
                P  = (*state_)[pr];

            if (on_land)
            {
                dAdA = (comb_*albf_*daFdA(A, Ta, P, i, j) - 1) / tauf_;
                dAdP = (comb_*albf_*daFdP(A, Ta, P, i, j)) / tauf_;
                dAdT = (comb_*albf_*daFdT(A, Ta, P, i, j)) / tauf_;
            }
            else
            {
                dAdA = -1 / tauc_;
                dAdP = 0.0;
                dAdT = 0.0;
            }

            AA_AA.set(i, j, l_, 5, dAdA);
            AA_PP.set(i, j, l_, 5, dAdP);
            AA_TT.set(i, j, l_, 5, dAdT);
        }

    Al_->set({1,n_,1,m_,1,l_,1,np_}, ATMOS_AA_, ATMOS_AA_, AA_AA);
    Al_->set({1,n_,1,m_,1,l_,1,np_}, ATMOS_AA_, ATMOS_PP_, AA_PP);
    Al_->set({1,n_,1,m_,1,l_,1,np_}, ATMOS_AA_, ATMOS_TT_, AA_TT);

    // Apply boundary conditions to stencil
    boundaries();

    // Assemble stencil into CRS struct
    assemble();

    if (!parallel_) buildLU_ = true;

    TIMER_STOP("AtmosLocal: compute Jacobian...");
}

//-----------------------------------------------------------------------------
std::shared_ptr<Utils::CRSMat> AtmosLocal::getJacobian()
{
    std::shared_ptr<Utils::CRSMat> jac = std::make_shared<Utils::CRSMat>();
    jac->co  = co_;
    jac->jco = jco_;
    jac->beg = beg_;
    return jac;
}

//-----------------------------------------------------------------------------
void AtmosLocal::computeMassMat()
{
    TIMER_START("AtmosLocal: compute mass matrix...");

    // Only TT and QQ equations. Integral equations give a zero
    // diagonal element. This is taken care of on the parallel side.
    int rowTT, rowQQ, rowAA;
    for (int i = 1; i <= n_; ++i)
        for (int j = 1; j <= m_; ++j)
            for (int k = 1; k <= l_; ++k)
            {
                rowTT = find_row(i, j, k, ATMOS_TT_)-1;
                rowQQ = find_row(i, j, k, ATMOS_QQ_)-1;
                rowAA = find_row(i, j, k, ATMOS_AA_)-1;
                (*diagB_)[rowTT] = Ai_;
                (*diagB_)[rowQQ] = 1.0;
                (*diagB_)[rowAA] = 1.0;
           }

    TIMER_STOP("AtmosLocal: compute mass matrix...");
}

//-----------------------------------------------------------------------------
void AtmosLocal::computeRHS()
{
    TIMER_START("AtmosLocal: compute RHS...");

    std::fill(rhs_->begin(), rhs_->end(), 0.0);

    // In parallel, computing E and P are directed by AtmospherePar.
    // Otherwise we do this here
    if (!parallel_)
    {
        computeEvaporation();   //
        computePrecipitation(); //
    }

    // Compute new Jacobian
    computeJacobian();

    // Compute the forcing
    forcing();

    // test heat fluxes
    int idx = 5;
    int jdx = 1;
    while ((*surfmask_)[(jdx-1)*n_+(idx-1)])
    {
        idx++; jdx++;
    }

    // Compute the right hand side rhs_
    double value;
    int row;
    for (int i = 1; i <= n_; ++i)
        for (int j = 1; j <= m_; ++j)
            for (int XX = 1; XX <= nun_; ++XX)
            {
                row   = find_row(i, j, l_, XX);
                value = 0.0;

                // Nonlinear albedo equation is at this point computed
                // in forcing. FIXME, this is a HACK. I need to think
                // about this a little longer.
                if (XX != ATMOS_AA_)
                    value += matvec(row);

                value += frc_[row-1];
                (*rhs_)[row-1] = value;
            }

    // Check integral condition (from matvec vs dot)
    if (!parallel_)
    {
        double integral = Utils::dot(*intcondCoeff_, *state_);

        if (std::abs((*rhs_)[rowIntCon_-1] - integral) > 1e-7)
        {
            std::cout << "rhs[rowIntCon] = " << (*rhs_)[rowIntCon_-1];
            std::cout << " integral = " << integral;
            std::cout << " rowIntCon = " << rowIntCon_ - 1 << std::endl;

            for (size_t i = 0; i != intcondCoeff_->size(); ++i)
                std::cout << (*intcondCoeff_)[i] << " " << (*state_)[i] << std::endl;

            ERROR("Error in integral condition", __FILE__, __LINE__);
        }
    }

    TIMER_STOP("AtmosLocal: compute RHS...");
}

//-----------------------------------------------------------------------------
double AtmosLocal::matvec(int row)
{
    TIMER_START("AtmosLocal: matvec...");

    // Returns inner product of a row in the matrix with the state.
    // > ugly stuff with 1 to 0 based...
    int first = beg_[row-1];
    int last  = beg_[row] - 1;
    double result = 0.0;
    for (int j = first; j <= last; ++j)
        result += co_[j-1] * (*state_)[jco_[j-1]-1];

    TIMER_STOP("AtmosLocal: matvec...");
    return result;
}

//-----------------------------------------------------------------------------
// Note that we are slightly messing up the philosophy here by adding local state
// dependencies to the forcing.
void AtmosLocal::forcing()
{
    double value, Ts, Eo, Ei;
    double Ta, P = 0.0, A, QSW;
    int tr, hr, sr, pr, ar; // indices
    bool on_land;

    if (std::abs(Ooa_) < 1e-8)
        WARNING(" Ooa_ too small", __FILE__, __LINE__);

    for (int j = 1; j <= m_; ++j)
        for (int i = 1; i <= n_; ++i)
        {
            tr = find_row(i, j, l_, ATMOS_TT_) - 1; // temperature row
            hr = find_row(i, j, l_, ATMOS_QQ_) - 1; // humidity row
            ar = find_row(i, j, l_, ATMOS_AA_) - 1; // albedo row

            if (aux_ == 1)
                pr = find_row(i, j, l_, ATMOS_PP_) - 1; // precipitation row
            else
                pr = -1;

            sr = n_*(j-1) + (i-1); // plain surface row

            // get albedo at this grid point (state component)
            A  = (*state_)[ar];

            // get atmospheric temp at this grid point (state component)
            Ta = (*state_)[tr];

            // ------------ Temperature forcing
            // Apply surface mask and calculate land temperatures.
            // This is a copy of legacy stuff, could be put more
            // clearly.

            // linear component of incoming shortwave radiation
            QSW = suna_[j] * (1 - a0_);

            // above land
            on_land = (*surfmask_)[sr];
            if (on_land)
            {
                // Simplified expression by equating sensible and
                // shortwave heat flux from the atmosphere into the
                // land.
                value = comb_ * sunp_ * suno_[j] * (1 - a0_) / Ooa_ ;

                // set forcing
                value += comb_ * (sunp_*QSW - lonf_*amua_);

                // set land temperature
                (*lst_)[sr] = Tl(A, Ta, j);
            }
            else // above ocean / sea ice
            {
                // Sensible heat flux surface component. Sea surface
                // temperature is corrected with sea ice surface
                // temperature (sit - sst) when Msi = 1. In that case
                // the background values need to be corrected as well
                // (t0i_ - t0o_).
                Ts = (*sst_)[sr] +
                    (*Msi_)[sr] * ((*sit_)[sr] - ((*sst_)[sr]) + t0i_ - t0o_);

                value = Ts + comb_ * (sunp_*QSW - lonf_*amua_);

                // latent heat due to precipitation (background contribution)
                value += comb_ * latf_ * lvscale_ * Pdist_[sr] * Po0_;
            }

            frc_[tr] = value;

            // ------------ Humidity forcing
            // Again, check whether we are above land
            if (on_land)
                value = 0;
            else
            {
                // Evaporation/sublimation forcing
                Eo = ( tdim_ / qdim_ ) * dqso_ * (*sst_)[sr];
                Ei = ( tdim_ / qdim_ ) * dqsi_ * (*sit_)[sr];
                value = nuq_ *  (Eo + (*Msi_)[sr]*(Ei - Eo + Cs_) );
            }

            frc_[hr] = value;

            // ------------ Albedo forcing
            //
            // Here we use the full nonlinear discretization as the
            // tanh switching behaviour cannot be linearized.

            if (on_land)
            {
                // The albedo parametrization accepts the global
                // precipitation anomaly (state component), but uses
                // the full dimensional and spatially distributed
                // value internally. Derivatives w.r.t. state can be
                // computed using finite differences, see e.g. daFdP()
                if (pr >= 0)
                    P  = (*state_)[pr];

                value = ( comb_ * albf_ * aF(A,Ta,P,i,j) - A ) / tauf_;
            }
            else
            {
                value = ( comb_ * albf_ * (*Msi_)[sr] - A )  / tauc_;
            }

            frc_[ar] = value;
        }

    // adjust to allow for integral condition in the serial case
    if (!parallel_)
        frc_[rowIntCon_-1] = 0.0;
}

//-----------------------------------------------------------------------------
void AtmosLocal::getFluxes(double *lwflux, double *swflux,
                           double *shflux, double *lhflux)
{
    // Temperature fluxes into the atmosphere:
    // -QLW + QSW + QSH + QLH
    int pos = 0;
    int sr,tr,ar,pr;
    double Ta,A,P;
    bool on_land;
    for (int j = 1; j <= m_; ++j)
        for (int i = 1; i <= n_; ++i)
        {
            sr = n_*(j-1) + (i-1); // plain surface row
            tr = find_row(i, j, l_, ATMOS_TT_) - 1; // temperature row
            ar = find_row(i, j, l_, ATMOS_AA_) - 1; // albedo row

            if (aux_ == 1)
                pr = find_row(i, j, l_, ATMOS_PP_) - 1; // precipitation row
            else
            {
                pr = -1;
                WARNING("Latent heat due to precipitation not defined", __FILE__, __LINE__);
            }

            Ta = (*state_)[tr]; // atmos temp
            A  = (*state_)[ar]; // albedo
            P  = (*state_)[pr]; // global precipitation

            // long wave radiative flux
            lwflux[pos] = -muoa_ * (comb_*lonf_*amua_ + bmua_*Ta);

            // short wave radiative flux
            swflux[pos] = comb_*sunp_*muoa_*suna_[j]*((1 - a0_) - da_*A);

            // sensible heat flux
            on_land = (*surfmask_)[sr];
            if (on_land)
                shflux[pos] = muoa_ * ((*lst_)[sr] - Ta);
            else
                shflux[pos] = muoa_ * ((*sst_)[sr] - Ta + (*Msi_)[sr] *
                                       ((*sit_)[sr] - ((*sst_)[sr]) + t0i_ - t0o_));

            // latent heat due to precipitation
            if ((pr >= 0) && (!on_land))
            {
                lhflux[pos] = comb_ * latf_ * rhoo_ * lv_ *
                    Pdist_[sr] * ( Po0_ + eta_ * qdim_ * P);
            }

            pos++;
        }
}

//-----------------------------------------------------------------------------
// Here we calculate the fully dimensional evaporation
void AtmosLocal::computeEvaporation()
{
    int hr, sr, ctr = 0;
    double Eo, Ei, M, q;

    for (int j = 1; j <= m_; ++j)
        for (int i = 1; i <= n_; ++i)
        {
            sr = n_*(j-1) + (i-1);
            hr = find_row(i, j, l_, ATMOS_QQ_) - 1;

            // reset vector element
            (*E_)[sr] = 0.0;

            if ((*surfmask_)[(j-1)*n_+(i-1)])
                continue; // do nothing

            // Mask value
            M = (*Msi_)[sr];

            // Humidity value
            q = (*state_)[hr];

            // Compute evaporation/sublimation based on surface
            // temperature (sst or sit).
            Eo = (tdim_ / qdim_) * dqso_ * (*sst_)[sr];
            Ei = (tdim_ / qdim_) * dqsi_ * (*sit_)[sr];
            (*E_)[sr] = Eo - q + M * (Ei - Eo + Cs_);

            // Create dimensional value
            (*E_)[sr] = Eo0_ + eta_ * qdim_ * (*E_)[sr];

            ctr++;
        }
}

//-----------------------------------------------------------------------------
// Only for serial use. In a parallel setting, precipitation is a
// global integral governed by AtmospherePar.
void AtmosLocal::computePrecipitation()
{
    if (parallel_)
    {
        WARNING("Function should not be called in parallel.", __FILE__, __LINE__);
        return; // do nothing
    }

    // Assuming E_ is dimensional, this returns a constant dimensional
    // P_. P_ may be spatially distributed with a function f (Pdist)
    // that satisfies int f dA = int 1 dA
    double integral = Utils::dot(*pIntCoeff_, *E_) / totalArea_;

    // fill Pdist array
    fillPdist(&Pdist_[0]);
    double intPdist = Utils::dot(*pIntCoeff_, Pdist_) / totalArea_;
    for (auto &e: Pdist_)
        if (std::abs(e) > 0.0)
            e += 1 - intPdist;
    
    int sr; // surface row
    for (int j = 0; j != m_; ++j)
        for (int i = 0; i != n_; ++i)
        {
            sr = j*n_+i;

            // reset vector element
            (*P_)[sr] = 0.0;

            if ((*surfmask_)[sr])
                continue; // leave it

            (*P_)[sr] = Pdist_[sr] * integral;
        }
}

//-----------------------------------------------------------------------------
double AtmosLocal::aF(double A, double Ta, double P, int i, int j)
{
    // Create dimensional P (m/y) including spatial
    // distribution
    double dimP = 3600. * 24. * 365. * Pdist_[n_*(j-1)+(i-1)] *
        (Po0_ + eta_ * qdim_ * P);

    return
        H(Tm_ - Tl(A,Ta,j), epm_) *
        H(Tr_ - Tl(A,Ta,j), epr_) *
        H(dimP - Pa_, epa_);

    // return 5*Ta + 4*A + 3*P;

    // return H(Tm_ - Tl(A,Ta,j), epm_) *
    //     H(Tr_ - Tl(A,Ta,j), epr_) *
    //     H(P - Pa_, epa_);

    // return H(Tm_ - Tl(A,Ta,j), epm_) *
    //     H(Tr_ - Tl(A,Ta,j), epr_) * H(P - Pa_, epa_);
}

//-----------------------------------------------------------------------------
void AtmosLocal::discretize(int type, Atom &atom)
{
    switch (type)
    {
        double val2, val4, val5, val6, val8;
        double cosdx2i, dy2i;
    case 1: // tc (without land points)
        atom.set({1,n_,1,m_,1,l_}, 5, 1.0);

        // Apply land mask
        for (int k = 1; k <= l_; ++k)
            for (int j = 1; j <= m_; ++j)
                for (int i = 1; i <= n_; ++i)
                    if ((*surfmask_)[(j-1)*n_+(i-1)])
                        atom.set(i, j, k, 5, 0.0);
        break;

    case 2: // tc2 (with land points)
        atom.set({1,n_,1,m_,1,l_}, 5, 1.0);
        break;

    case 3: // txx
        for (int i = 1; i <= n_; ++i)
            for (int j = 1; j <= m_; ++j)
            {
                cosdx2i = 1.0 / pow(cos(yc_[j]) * dx_, 2);
                val2 = datc_[j] * cosdx2i;
                val8 = val2;
                val5 = -2 * val2;

                for (int k = 1; k <= l_; ++k)
                {
                    atom.set(i, j, k, 2, val2);
                    atom.set(i, j, k, 8, val8);
                    atom.set(i, j, k, 5, val5);
                }
            }
        break;

    case 4: // tyy
        dy2i = 1.0 / pow(dy_, 2);
        for (int i = 1; i <= n_; ++i)
            for (int j = 1; j <= m_; ++j)
            {
                val4 = dy2i * datv_[j-1] * cos(yv_[j-1]) / cos(yc_[j]);
                val6 = dy2i * datv_[j]   * cos(yv_[j])   / cos(yc_[j]);
                val5 = -(val4 + val6);

                for (int k = 1; k <= l_; ++k)
                {
                    atom.set(i, j, k, 4, val4);
                    atom.set(i, j, k, 6, val6);
                    atom.set(i, j, k, 5, val5);
                }
            }
        break;

    case 5: // qxx
        for (int i = 1; i <= n_; ++i)
            for (int j = 1; j <= m_; ++j)
            {
                cosdx2i = 1.0 / pow(cos(yc_[j]) * dx_, 2);
                val2 = cosdx2i;
                val8 = val2;
                val5 = -2 * val2;

                for (int k = 1; k <= l_; ++k)
                {
                    atom.set(i, j, k, 2, val2);
                    atom.set(i, j, k, 8, val8);
                    atom.set(i, j, k, 5, val5);
                }
            }
        break;

    case 6: // tyy
        dy2i = 1.0 / pow(dy_, 2);
        for (int i = 1; i <= n_; ++i)
            for (int j = 1; j <= m_; ++j)
            {
                val4 = dy2i * cos(yv_[j-1]) / cos(yc_[j]);
                val6 = dy2i * cos(yv_[j])   / cos(yc_[j]);
                val5 = -(val4 + val6);

                for (int k = 1; k <= l_; ++k)
                {
                    atom.set(i, j, k, 4, val4);
                    atom.set(i, j, k, 6, val6);
                    atom.set(i, j, k, 5, val5);
                }
            }
        break;
    }
}


//-----------------------------------------------------------------------------
void AtmosLocal::assemble()
{
    // Create CRS matrix storage and/or padded banded storage

    // clear old CRS matrix
    beg_.clear();
    co_.clear();
    jco_.clear();

    // Clear banded storage (just to be sure)
    if (!parallel_) std::fill(bandedA_.begin(), bandedA_.end(), 0.0);

    int i2,j2,k2; // will contain neighbouring grid pointes given by shift()
    int row;
    int rowb; // for banded storage
    int col;
    int colb; // for banded storage
    int idx;
    int kdiag = ksub_ + ksup_ + 1; // for banded storage
    int elm_ctr = 1;
    double value;
    for (int k = 1; k <= l_; ++k)
        for (int j = 1; j <= m_; ++j)
            for (int i = 1; i <= n_; ++i)
                for (int A = 1; A <= nun_; ++A)
                {
                    // Filling new row:
                    //  find the row corresponding to A at (i,j,k):
                    row = find_row(i, j, k, A);
                    //  put element counter in beg:
                    beg_.push_back(elm_ctr);
                    for (int loc = 1; loc <= np_; ++loc)
                    {
                        // find index of neighbouring point loc
                        shift(i,j,k,i2,j2,k2,loc);
                        for (int B = 1; B <= (nun_ + aux_); ++B)
                        {
                            value = Al_->get(i,j,k,loc,A,B);
                            if (std::abs(value) > 0)
                            {
                                // CRS --------------------------------------
                                co_.push_back(value);
                                col = find_row(i2,j2,k2,B);
                                jco_.push_back(col);

                                // increment the element counter
                                ++elm_ctr;

                                if (!parallel_)
                                {
                                    // BND --------------------------------------
                                    // get row index for banded storage
                                    rowb = row - col + kdiag;

                                    // put matrix values in column major fashion
                                    // in the array
                                    //  > go from 1 to 0-based
                                    colb = col;
                                    idx  = rowb + (colb - 1) * ldimA_ - 1;
                                    bandedA_[idx] = value;
                                }
                            }
                        }
                    }
                }

    // final element of beg
    beg_.push_back(elm_ctr);

    // in the serial case we replace the final row with coefficients for the integral condition
    if (!parallel_) intcond();

}

//----------------------------------------------------------------------------
void AtmosLocal::intcond()
{
    // We replace a row in the crs struct
    int first = beg_[rowIntCon_ - 1];
    int last  = beg_[rowIntCon_];

    // create tmp beg array relative to integral condition row
    std::vector<int> tmp;
    for (int i = rowIntCon_; i < (int) beg_.size(); ++i)
        tmp.push_back(beg_[i]-beg_[rowIntCon_]);

    // erase intcon row in co_ and jco and end pointer in beg
    co_.erase  ( co_.begin() + first - 1,  co_.begin() + last - 1);
    jco_.erase (jco_.begin() + first - 1, jco_.begin() + last - 1);
    beg_.erase (beg_.begin() + rowIntCon_);

    // append crs arrays with integral coefficients and column indices
    std::vector<double> vals,inds;
    integralCoeff(vals, inds);

    co_.insert  (co_.begin() + first - 1, vals.begin(), vals.end());
    jco_.insert (jco_.begin() + first - 1, inds.begin(), inds.end());
    beg_.insert (beg_.begin() + rowIntCon_, first + vals.size() );

    // let the rest of the beg array be relative to the new integral
    // condition row
    for (int i = rowIntCon_; i < (int) beg_.size(); ++i)
        beg_[i] = tmp[i - rowIntCon_] + beg_[rowIntCon_];
}

//-----------------------------------------------------------------------------
// Declaring a few LAPACK functions needed around this point

// Solve banded system stored in bandedA_
extern "C" void dgbsv_(int *N, int *KL, int *KU, int *NRHS, double *AB,
                       int *LDAB, int *IPIV, double *B,
                       int *LDB, int *INFO);

// Create LU factorization of banded system
extern "C" void dgbtrf_(int *M, int *N, int *KL, int *KU, double *AB,
                        int *LDAB, int *IPIV, int *INFO);

// Solve system using LU factorization given by dgbtrf
extern "C" void dgbtrs_(char *TRANS, int *N, int *KL, int *KU, int *NRHS,
                        double *AB, int *LDAB, int *IPIV, double *B, int *LDB,
                        int *INFO);

//-----------------------------------------------------------------------------
// In parallel setup we use Ifpack
void AtmosLocal::solve(std::shared_ptr<std::vector<double> > const &rhs)
{
    if (parallel_)
    {
        ERROR("Non-parallel solve called in parallel setting",
              __FILE__, __LINE__);
        return;
    }

    TIMER_START("AtmosLocal: solve...");

    int dim     = n_*m_*l_;
    int nrhs    = 1;
    int ldb     = dim;
    int info;

    if (buildLU_)
    {
        TIMER_START("AtmosLocal: build LU (dgbtrf)");
        dgbtrf_(&dim, &dim, &ksub_, &ksup_, &bandedA_[0],
                &ldimA_, &ipiv_[0], &info);
        buildLU_ = false; // until next request
        TIMER_STOP("AtmosLocal: build LU (dgbtrf)");
    }

    // copy rhs into sol
    *sol_  = *rhs;

    TIMER_START("AtmosLocal: solve (dgbtrs)");
    char trans  = 'N';

    // at entry sol contains the rhs, at exit it contains the solution
    dgbtrs_(&trans, &dim, &ksub_, &ksup_, &nrhs,
            &bandedA_[0], &ldimA_, &ipiv_[0],  &(*sol_)[0], &ldb, &info);

    TIMER_STOP("AtmosLocal: solve (dgbtrs)");

    TIMER_STOP("AtmosLocal: solve...");
}

//------------------------------------------------------------------
std::shared_ptr<std::vector<double> >
AtmosLocal::getCurrResVec(std::shared_ptr<std::vector<double> > const &x,
                          std::shared_ptr<std::vector<double> > const &rhs)
{
    std::shared_ptr<std::vector<double> > r =
        std::make_shared<std::vector<double> >(*x);

    applyMatrix(*x, *r);
    Utils::update(1.0, *rhs, -1.0, *r);
    return r;
}

//-----------------------------------------------------------------------------

//! Size of stencil/neighbourhood:
//! +----------++-------++----------+
//! | 12 15 18 || 3 6 9 || 21 24 27 |
//! | 11 14 17 || 2 5 8 || 20 23 26 |
//! | 10 13 16 || 1 4 7 || 19 22 25 |
//! |  below   || center||  above   |
//! +----------++-------++----------+

//!  No flow boundaries for both TT an QQ. For instance, dependencies
//!  on non-existent grid points at the western boundary are added to
//!  dependencies on the central grid point as T_{i-1} = T_{i}.

void AtmosLocal::boundaries()
{

    int west, east, north, south;
    for (int i = 1; i <= n_; ++i)
        for (int j = 1; j <= m_; ++j)
            for (int k = 1; k <= l_; ++k)
                for (int XX = 1; XX <= nun_; ++XX)
                {
                    west   = i-1;
                    east   = i+1;
                    north  = j+1;
                    south  = j-1;

                    // western boundary
                    if (west == 0 && !periodic_)
                    {
                        Al_->set(i,j,k,5,XX,XX,
                                 Al_->get(i,j,k,5,XX,XX) +
                                 Al_->get(i,j,k,2,XX,XX));
                        Al_->set(i,j,k,2,XX,XX, 0.0);
                    }

                    // eastern boundary
                    if (east == n_+1 && !periodic_)
                    {
                        Al_->set(i,j,k,5,XX,XX,
                                 Al_->get(i,j,k,5,XX,XX) +
                                 Al_->get(i,j,k,8,XX,XX));
                        Al_->set(i,j,k,8,XX,XX, 0.0);
                    }

                    // northern boundary
                    if (north == m_+1)
                    {
                        Al_->set(i,j,k,5,XX,XX,
                                 Al_->get(i,j,k,5,XX,XX) +
                                 Al_->get(i,j,k,6,XX,XX));
                        Al_->set(i,j,k,6,XX,XX, 0.0);
                    }

                    // southern boundary
                    if (south == 0)
                    {
                        Al_->set(i,j,k,5,XX,XX,
                                 Al_->get(i,j,k,5,XX,XX) +
                                 Al_->get(i,j,k,4,XX,XX));
                        Al_->set(i,j,k,4,XX,XX, 0.0);
                    }
                }
}

//-----------------------------------------------------------------------------
void AtmosLocal::write(std::vector<double> &vector, const std::string &filename)
{
    if (!vector.empty())
    {
        std::ofstream atmos_ofstream;
        atmos_ofstream.open(filename);
        for (auto &i : vector)
            atmos_ofstream << std::setprecision(12) << i << '\n';
        atmos_ofstream.close();
    }
    else
        WARNING("vector is empty", __FILE__, __LINE__);
}

//-----------------------------------------------------------------------------
int AtmosLocal::find_row(int i, int j, int k, int XX)
{
    ////////////////////////////////////////////
    // 1-based
    ////////////////////////////////////////////

    // P values are auxiliary and come after the final ordinary element
    if ( XX >= ATMOS_PP_ )
    {
        if (aux_ <= 0)
        {
            ERROR("AtmosLocal: invalid row", __FILE__, __LINE__);
            return -1;
        }
        else
            return dim_ - aux_ + (XX - ATMOS_PP_) + 1;
    }
    else // ordinary 1-based find_row
        return nun_ * ((k-1)*n_*m_ + n_*(j-1) + (i-1)) + XX;
}

//-----------------------------------------------------------------------------
// --> Many of these points are never used. For the atmosphere model
// we could use less neighbours and perhaps win a few ms.
void AtmosLocal::shift(int i, int j, int k,
                       int &i2, int &j2, int &k2, int loc)
{
    if (loc < 10)
    {
        k2 = k;
        //   +-------+    +---------+
        //   | 3 6 9 |    | 1  1  1 |
        //   | 2 5 8 | -> | 0  0  0 |
        //   | 1 4 7 |    |-1 -1 -1 |
        //   +-------+    +---------+
        j2 = j + ((loc + 2) % 3) - 1;
        //   +-------+    +---------+
        //   | 3 6 9 |    |-1  0  1 |
        //   | 2 5 8 | -> |-1  0  1 |
        //   | 1 4 7 |    |-1  0  1 |
        //   +-------+    +---------+
        i2 = i + (loc - 1) / 3 - 1;
    }
    else if (loc < 19)
    {
        k2 = k - 1;
        //   +----------+     +---------+
        //   | 12 15 18 |     | 1  1  1 |
        //   | 11 14 17 | ->  | 0  0  0 |
        //   | 10 13 16 |     |-1 -1 -1 |
        //   +----------+     +---------+
        j2 = j + ((loc + 2) % 3) - 1;
        //   +----------+     +---------+
        //   | 12 15 18 |     |-1  0  1 |
        //   | 11 14 17 | ->  |-1  0  1 |
        //   | 10 13 16 |     |-1  0  1 |
        //   +----------+     +---------+
        i2 = i + (loc - 10) / 3 - 1;
    }
    else
    {
        k2 = k + 1;
        //   +----------+     +---------+
        //   | 21 24 27 |     | 1  1  1 |
        //   | 20 23 26 | ->  | 0  0  0 |
        //   | 19 22 25 |     |-1 -1 -1 |
        //   +----------+     +---------+
        j2 = j + ((loc + 2) % 3) - 1;
        //   +----------+     +---------+
        //   | 21 24 27 |     |-1  0  1 |
        //   | 20 23 26 | ->  |-1  0  1 |
        //   | 19 22 25 |     |-1  0  1 |
        //   +----------+     +---------+
        i2 = i + (loc - 19) / 3 - 1;
    }

    if (periodic_)
    {
        if (i2 == n_+1)
            i2 = 1;
        if (i2 == 0)
            i2 = n_;
    }
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > AtmosLocal::getSolution(char mode)
{
    return Utils::getVector(mode, sol_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > AtmosLocal::getState(char mode)
{
    return Utils::getVector(mode, state_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > AtmosLocal::getRHS(char mode)
{
    return Utils::getVector(mode, rhs_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > AtmosLocal::getMassMat(char mode)
{
    return Utils::getVector(mode, diagB_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > AtmosLocal::getSST(char mode)
{
    return Utils::getVector(mode, sst_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > AtmosLocal::interfaceE(char mode)
{
    return Utils::getVector(mode, E_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > AtmosLocal::interfaceP(char mode)
{
    return Utils::getVector(mode, P_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > AtmosLocal::getPIntCoeff(char mode)
{
    return Utils::getVector(mode, pIntCoeff_);
}

//-----------------------------------------------------------------------------
void AtmosLocal::applyMatrix(std::vector<double>  const &v,
                             std::vector<double>  &out)
{
    int first;
    int last;

    TIMER_START("AtmosLocal: apply matrix");
    // Perform matrix vector product
    // 1->0 based... horrible...
    for (size_t row = 1; row <= v.size(); ++row)
    {
        first = beg_[row-1];
        last  = beg_[row] - 1;

        out[row-1] = 0;
        for (int col = first; col <= last; ++col)
            out[row-1] += co_[col-1] * v[jco_[col-1]-1];
    }
    TIMER_STOP("AtmosLocal: apply matrix");
}

// ---------------------------------------------------------------------------
// Adjust specific parameter
void AtmosLocal::setPar(std::string const &parName, double value)
{
    int ctr = 0;
    if (parName.compare(allParameters_[ctr++]) == 0)
        comb_ = value;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        sunp_ = value;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        lonf_ = value;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        humf_ = value;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        latf_ = value;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        albf_ = value;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        tdif_ = value;

    // The effects of comb_ and humf_ are combined in nuq_, so here we
    // update nuq_.
    nuq_ = comb_* humf_ * (eta_ / hdimq_ ) * (rhoo_ / rhoa_) * (r0dim_ / udim_);

    // If parameter not available we take no action
}

// ---------------------------------------------------------------------------
// Get parameter value
double AtmosLocal::getPar(std::string const &parName)
{
    int ctr = 0;
    if (parName.compare(allParameters_[ctr++]) == 0)
        return comb_;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        return sunp_;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        return lonf_;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        return humf_;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        return latf_;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        return albf_;
    else if (parName.compare(allParameters_[ctr++]) == 0)
        return tdif_;
    else // If parameter not available we return 0
        return 0;
}

//-----------------------------------------------------------------------------
std::string AtmosLocal::int2par(int ind) const
{
    return allParameters_[ind];
}

//-----------------------------------------------------------------------------
int AtmosLocal::par2int(std::string const &label) const
{
    for (size_t i = 0; i < allParameters_.size(); ++i)
    {
        if (label == allParameters_[i])
            return (int) i;
    }
    WARNING("Parameter name " << label << " is not represented in AtmosLocal",
            __FILE__, __LINE__);
    return -1;
}

//-----------------------------------------------------------------------------
void AtmosLocal::setSurfaceMask(std::shared_ptr<std::vector<int> > surfm)
{
    // clear current mask
    surfmask_->clear();

    if ((int) surfm->size() < (m_* n_))
    {
        ERROR("surfm->size() not ok:",  __FILE__, __LINE__);
    }
    else if ((int) surfm->size() > (m_ * n_))
    {
        // in this case we assume we receive an ocean landmask
        // with boundaries, which implies that the final
        // 2 * (n_+2) * (m_+2) entries are meaningful for us.
        // The final (n_+2) * (m_+2) entries contain ones.
        int maskdim = 2 * (n_+2) * (m_+2);

        surfm->erase(surfm->begin(), surfm->begin() + surfm->size() - maskdim);

        // now we put the first layer of surfm in our datamember, without borders
        for (int j = 1; j != m_+1; ++j)
            for (int i = 1; i != n_+1; ++i)
            {
                surfmask_->push_back((*surfm)[j*(n_+2) + i]);
            }

    }
    else // we trust surfm
    {
        surfmask_ = surfm;
    }

// #ifdef DEBUGGING_NEW
//     Utils::printSurfaceMask(surfmask_, "surfmask", n_);
// #endif
}
