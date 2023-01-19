# Continuation with flux forcing
  (this is a draft, steps described here need testing)
  
    1.  Start a continuation with a non-zero salinity forcing.
    2.  Converge at some desired equilibrium. By default the ocean component saves the salinity flux (see m_probe::get_salflux). This flux can be used as a non-restoring forcing (SRES=0) and will give a zero residual for the same equilibrium state.
    3.  Now it's possible to restart a continuation at this equilibrium that has a flux forcing for salinity. For this, put SRES=0 (perhaps to be added in the xml) and Load salinity flux = true. This continuation can be in any desired parameter.
