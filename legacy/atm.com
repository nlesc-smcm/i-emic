      real, parameter :: hdima = 8400.
      real, parameter ::  rhoa = 1.25
      real, parameter ::  uatm = 0.0
      real, parameter ::    ce = 1.3e-03
      real, parameter ::    ch = 0.94 * ce  ! 
      real, parameter ::   cpa = 1000.
       real, parameter ::    uw = 8.5 ! -> mu=2.5
      real, parameter ::    d0 = 3.1e+06
      real, parameter ::    c0 = 0.43 
      real, parameter :: sigmab = 5.67e-08
      real, parameter ::   brad = 1.5
      real, parameter ::   arad = 216.0
      real, parameter ::   sun0 = 1360.

      real    Ai, Ad, As, Aa, Os, Aoa, Ooa, amua
      real    dat, davt, albe, suna, suno,upa,bmua
      common /atm/ dat(m),davt(0:m),
     >     albe(m),suna(m),suno(m),Ai, Ad, As, Os, Ooa,
     >     amua,bmua,Aa,upa(m)
*LtR 21/11/01      real, parameter ::    ch = 0.0615 * ce ! muoa=0.5
*LtR 15/11/01      real, parameter ::    ch = 0.308 * ce  ! muoa=2.5
*LtR 05/11/01      real, parameter ::    ch = 0.94 * ce   ! muoa=7.6
*      real, parameter ::    uw = 31.4 !  -> mu=48.0
*      real, parameter ::    uw = 26.2 ! -> mu=40.0
*      real, parameter ::    uw = 13.1 ! -> mu=20.0
*      real, parameter ::    uw = 6.55 ! -> mu=10.0
*      real, parameter ::    uw = 3.275 ! -> mu=5.0
*      real, parameter ::    uw = 5. !   -> mu=7.6
*      real, parameter ::    uw = 1.64 ! -> mu=2.5
*LtR 19/11/01      real, parameter ::    d0 = 0.4e+06
*LtR 16/11/01      real, parameter ::    d0 = 0.1e+06
*LtR 16/11/01      real, parameter ::    d0 = 0.7e+06
*LtR 05/11/01      real, parameter ::    d0 = 1.5e+06
*LtR 15/11/01      real, parameter ::    d0 = 1.0e+06
* LtR 30/11/01      real, parameter ::   arad = 216.0
* LtR 30/11/01      real, parameter ::   brad = 1.5
