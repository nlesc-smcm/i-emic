#include "SeaIce.H"
#include "SeaIceDefinitions.H"

SeaIce::SeaIce(Teuchos::RCP<Epetra_Comm> comm, ParameterList params)
    :
    params_          (params),
    comm_            (comm),
    n_               (params->get("Global Grid-Size n", 16)),
    m_               (params->get("Global Grid-Size m", 16)),
    l_               (params->get("Global Grid-Size l", 1)),
    periodic_        (params->get("Periodic", false)),

    taus     (0.01),    // threshold ice thickness
    epsilon  (1e2),     // Heavyside approximation steepness
    
    // background mean values
    t0o   (params->get("background ocean temp t0o", 7)),
    t0a   (params->get("background atmos temp t0a", 10)),
    t0i   (params->get("background seaice temp t0i",-15)),
    tvar  (params->get("ocean temp variation tvar", 15)),
    s0    (params->get("ocean background salinity s0", 35)),
    svar  (params->get("ocean salinity variation svar", 1)),
    q0    (params->get("atmos background humidity q0", 1e-3)),
    qvar  (params->get("atmos humidity variation qvar", 5e-4)),
    H0    (params->get("seaice background thickness H0", taus)),
    M0    (params->get("seaice background mask M0", 0)),

    // ice formation constants
    ch    (params->get("empirical constant", 0.0058)),
    utau  (params->get("skin friction velocity, ms^{-1}", 0.02)),
    rhoo  (params->get("sea water density, kg m^{-3}", 1.024e3)),
    rhoi  (params->get("ice density, kg m^{-3}", 0.913e3)),
    rhoa  (params->get("atmospheric density, kg m^{-3}", 1.25)),
    cpo   (params->get("sea water heat capacity, W s kg^{-1} K^{-1]", 4.2e3)),
    Lf    (params->get("latent heat of fusion of ice, J kg^{-1}", 3.347e5)),
    Ls    (params->get("latent heat of sublimation of ice, J kg^{-1}", 2.835e6)),
    Ic    (params->get("constant ice conductivity, W m^{-1} K^{-1}", 2.166)),

    // combined parameter
    zeta  (ch * utau * rhoo * cpo),

    // sublimation constants, parameters for saturation humidity over
    // ice
    c1    (params->get("c1", 3.8e-3)),
    c2    (params->get("c2", 21.87)),
    c3    (params->get("c3", 265.5))
    
{

    xmin_ = params->get("Global Bound xmin", 286.0) * PI_ / 180.0;
    xmax_ = params->get("Global Bound xmax", 350.0) * PI_ / 180.0;
    ymin_ = params->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
    ymax_ = params->get("Global Bound ymax", 74.0)  * PI_ / 180.0;

    dof_  = SEAICE_NUN_;

    // Create domain object
    domain_ = Teuchos::rcp(new TRIOS::Domain(n_, m_, l_, dof_,
                                             xmin_, xmax_, ymin_, ymax_,
                                             periodic_, 1.0, comm_));
    
    
}
