      include 'par.com'

      real, parameter :: pi    =  3.14159265358979323846

! North Atlantic 
      real, parameter :: xmin          =  300. * pi/180.
      real, parameter :: xmax          =  340. * pi/180.
      real, parameter :: ymin          =   20. * pi/180.
      real, parameter :: ymax          =   60. * pi/180.
      real, parameter :: hdim          =  4000.
      real, parameter :: zmin          = -1.
      real, parameter :: zmax          =  0.
      real, parameter :: qz            =  2.0

      integer, parameter :: SLIP       = -1      ! -1: no slip on lateral boundaries; 1: slip 
      integer, parameter :: TRES       =  1      ! 1: restoring temp forcing; 0: prescribed flux 
      integer, parameter :: SRES       =  0      ! 1: restoring salt forcing;  0: prescribed flux  
      integer, parameter :: itopo      =  1      ! 0: data, 1: no land, 2: idealized
      integer, parameter :: ih         =  0      ! 0: no inh. mixing, 1: yes
      integer, parameter :: ifw        =  1      ! freshwater forcing 0: data, 1: idealized
      logical, parameter :: rd_mask    = .false.
      logical, parameter :: FLAT       = .true.  ! flat bottom, otherwise .false.
      integer, parameter :: iza        =  0      ! wind forcing 0: data, 1: zon ave., 2: idealized
      integer, parameter :: vmix_GLB   =  1      ! 0: no mixing and CA; 1: full mixing
      integer, parameter :: tap        =  1      ! taper 1: Gerdes et al. 2: Danabasoglu & McWilliams
                                                 ! 3: De Niet el al.
      logical, parameter :: rho_mixing = .false. ! stable(new) = .true.  
      logical, parameter :: frs        = .false. ! 0: no free surface; 1: full free surface
      
      character(len=*), parameter ::  topdir = '/orbit_usr1/data/'

!    basic parameters 
      real, parameter :: omegadim = 7.292e-05 
      real, parameter :: r0dim    = 6.37e+06
      real, parameter :: udim     = 0.1e+00
      real, parameter :: gdim     = 9.8e+00
      real, parameter :: rhodim   = 1.0e+03
      real, parameter :: t0       = 15.0
      real, parameter :: s0       = 35.0
      real, parameter :: cp0      = 4.2e+03
!    parameters equation of state  !b1 = 5.6e-05, b2 = 6.3e-06, b3 = 3.7e-08
      real, parameter :: alphaT   = 1.0e-04
      real, parameter :: alphaS   = 7.6e-04
      real, parameter :: alpt1    = 2.93       ! 2*b2*t0/b1 - 3*b3*t0*t0/b1
      real, parameter :: alpt2    = 8.3e-02    ! b2/b1 - 3*t0*b3/b1
      real, parameter :: alpt3    = 6.6e-04    ! b3/b1
!    mixing parameters
      real, parameter :: ah       = 2.5e+05
      real, parameter :: av       = 1.0e-03
      real, parameter :: kappah   = 1.0e+03
      real, parameter :: kappav   = 1.0e-04

      real    dx, dy, dz, x, y, z, xu, yv, zw, ze, zwe, dfzT, dfzW
      logical periodic
      common /grd1/ dx, dy, dz, x(n), y(m), z(l), xu(0:n), yv(0:m)
      common /grd2/ zw(0:l), ze(l), zwe(l), dfzT(l), dfzW(0:l), periodic

      integer      landm, rowintcon
      common /lan/ landm(0:n+1,0:m+1,0:l+la+1), rowintcon
      	
      real         fricum,fricvm
      common /fri/ fricum(0:n,0:m,l),fricvm(0:n,0:m,l)

      real         kapv(0:l),kaph(l),emix(n,m,0:l)
      common /kap/ kapv,kaph,emix

      integer      iout
      common /out/ iout

      real    Frc, taux, tauy, tx, ty, ft, fs, tmax, tatm, emip
      real    ftlev, fslev
      common /for1/ Frc(ndim), taux(n,m), tauy(n,m), tx(n,m), ty(n,m)
      common /for2/ ft(n,m), fs(n,m), tmax, tatm(n,m), emip(n,m)
      common /for3/ ftlev(n,m,l), fslev(n,m,l)

!      real    Al(n,m,l+la,np,nun,nun)
!      common /ope/ Al

      real    pv_adj(n,m,0:l)
      common /pv_adjust/ pv_adj
!________________________________________________________________

      integer, parameter :: UU = 1
      integer, parameter :: VV = 2
      integer, parameter :: WW = 3
      integer, parameter :: PP = 4
      integer, parameter :: TT = 5
      integer, parameter :: SS = 6

      integer, parameter :: OCEAN = 0 
      integer, parameter :: LAND  = 1
      integer, parameter :: WATER = 2
      integer, parameter :: PERIO = 3 
      integer, parameter :: ATMOS = 4 

!________________________________________________________________
!     parameter for predictor
      integer, parameter :: predictor = 2 ! 1 = Euler, 2 = Secant
!________________________________________________________________
!     parameters for corrector
      integer, parameter :: NewtonMethod  = 2
      ! 1 = Newton, 2 = Adaptive Shamanskii, 3 = Newton Chord
      integer, parameter :: OptNewtonIter = 5

!________________________________________________________________
!     parameters for preconditioning and solution method
      
      logical, parameter :: scale   = .true. ! scaling of Jacobian
      integer, parameter :: prec    = 2      ! 1 = mrilu, 2 = block-ilu preconditioner
      integer, parameter :: variant = 1      ! 1 = TS, 2 = uv, 3 = w/p - schurcomplement
      integer, parameter :: spptype = 2      ! 1 = art.compr. 2 = modified simpler
      ! droptols for mrilu
      real, parameter :: epsw_mrilu = 1D-2   ! full Jacobian (only if prec = 1)
      real, parameter :: epsw_uv    = 1D-4   ! subblock Auv (if prec = 2)
      real, parameter :: epsw_muv   = 1D-3   ! depth-averaged Auv (if spptype = 1)      
      real, parameter :: epsw_p     = 1D-4   ! modified simpler approach (if prec = 2)
      real, parameter :: epsw_TS    = 4D-3   ! subblock ATS (if prec = 2)
      integer, parameter :: maxiter = 200    ! maximum size krylov spaces
