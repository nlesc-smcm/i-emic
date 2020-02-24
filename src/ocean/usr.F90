!! A Fortran 90 replacement for the file 'usr.com'.
!!
!! Note that the dimensions n,m,l are now set in the function init
!! (implemented in usrc.F90) instead of being defined in par.com.
!! Some constants have been moved from m_usr to m_par, if you don't
!! find what you are looking for, refer to par.F90.
!!
module m_usr

  use m_par

  !===== LOCAL GRID VARIABLES =================================================
  !
  !     The grid variables on the subdomain are set by init in usrc.F90
  !

  integer :: n    = 0             ! east-west (x) direction
  integer :: m    = 0             ! north-south (y) direction
  integer :: l    = 0             ! z direction
  integer :: ndim = 0             ! total number of unknowns n*m*l*nun

  real  :: xmin,xmax            ! limits in x direction
  real  :: ymin,ymax            ! limits in y direction

  logical :: periodic             ! east-west periodicity on subdomain

  real    :: dx, dy, dz
  real,    allocatable, dimension(:)     :: x, y, z, xu, yv
  real,    allocatable, dimension(:)     :: zw, ze, zwe, dfzT, dfzW

  integer, allocatable, dimension(:,:,:) :: landm

  !     This is only used if SRES=0 (non-restoring salinity forcing),
  !     we set it to -1 to disable the integral condition
  !     Integral condition is determined by c++ code
  integer :: rowintcon = -1

  !===== GLOBAL GRID VARIABLES =================================================
  !     The grid specifications that are the same between subdomains and
  !     global domain. If adjustable, it is done in initialize in
  !     global.F90

  real, parameter :: zmin  = -1.0
  real, parameter :: zmax  =  0.0

  real    :: hdim                 ! largest depth of the ocean
  real    :: qz        = 1.0      ! vertical stretching parameter

  integer :: itopo     = 1        ! 0: data, 1: no land, 2: idealized
  logical :: FLAT      = .false.  ! flat bottom, otherwise .false.
  logical :: rd_mask   = .false.  ! read land-array from file in ./mkmask/

  !===== MIXING AND FORCING FLAGS ==============================================
  !     These are set in initialize in global.F90

  integer :: ih         = 0       ! inhomogeneous (equatorial) mixing
  integer :: vmix       = 1       ! mixing flag
  integer :: tap        = 1       ! neutral physics taper, 1: Gerdes
  ! et al., 2: Danabasoglu & McWilliams,
  ! 3: De Niet et al.
  logical :: rho_mixing = .false. ! mix density, instead of T and S

  integer :: TRES       = 1       ! restoring temp forcing
  integer :: SRES       = 1       ! restoring salt forcing
  integer :: iza        = 2       ! wind forcing, 0: data, 1: zon ave.,
                                  !               2: idealized
  integer :: its        = 1       ! salinity, 0: data, 1: idealized
  integer :: ite        = 1       ! temperature, 0: data, 1:idealized

  logical :: rd_spertm  = .false. ! read salinity perturbation mask

  integer :: coriolis_on = 1      ! Enables Coriolis force, which is disabled
                                  ! in the 2DMOC setup.

  integer :: forcing_type = 0     ! Forcing type: 0: default, 1: 2DMOC,
                                  ! 2: Northern hemisphere

  !--obsolete---
  !integer :: ifw        = 1       ! freshwater forcing 0: data, 1: idealized
                                   ! (USE its NOW)

  !integer :: CVT        = 1       ! no = 1; tanh =2; arctan = 3
  !integer :: FBT        = 0       ! fixed bottom temperature 1 , no flux 0
  !integer :: FBS        = 0       ! fixed bottom salinity +1 , no flux 0

  !===== I-EMIC COUPLING FLAGS ==============================================
  !     These are set in global.F90::initialize 
  
  integer :: coupled_T = 0         ! 0: uncoupled temperature
                                   ! 1: temperature coupled to atmosphere
  integer :: coupled_S = 0         ! 0: uncoupled salinity
                                   ! 1: salinity coupled to atmosphere
  
  !===== FORCING FIELDS ========================================================
  !
  !     The subdomain forcing fields are mostly set by forcing, which is called
  !     from init.
  !
  !   Set by an external coupled model (see inserts.F90):
  !     tatm, qatm, albe, patm, msi, qsi, qsa 

  real, allocatable, dimension(:)     :: Frc
  real, allocatable, dimension(:,:)   :: taux, tauy
  real, allocatable, dimension(:,:)   :: tatm, emip, spert, adapted_emip
  real, allocatable, dimension(:,:)   :: qatm, albe, patm
  real, allocatable, dimension(:,:)   :: msi, gsi, qsa
  real, allocatable, dimension(:,:,:) :: internal_temp, internal_salt
  real, allocatable, dimension(:,:,:) :: ftlev, fslev
  !--obsolete---
  real, allocatable, dimension(:,:)   :: tx,ty,ft,fs

  !     Nondimensionalization coefficients for temperature and
  !     salinity body forcing.
  real :: QTnd, QSnd 

  !===== OUTPUT ================================================================
  !     Output unit, formerly sent to fort.99. It is not yet clear how we can
  !     preserve THCM output on the subdomains meaningfully as the fortran files
  !     seem to get lost between calls from C++.
#ifdef THCM_STDOUT
  integer :: f99  = 6
#else
  integer :: f99  = 99
#endif
  !     This is an artefact used in forcing.f to indicate wether to create some
  !     output. It should probably be kicked out
  integer :: iout = 0

  !===== FIXED PARAMETERS ======================================================

  ! real, parameter :: omegadim = 7.272e-05  ! 2DMOC
  real, parameter :: omegadim = 7.292e-05
  ! real, parameter :: r0dim    = 6.371e+06 ! 2DMOC
  real, parameter :: r0dim    = 6.37e+06 
  real, parameter :: udim     = 0.1e+00
  real, parameter :: gdim     = 9.8e+00
  real, parameter :: rhodim   = 1.024e+03
  real, parameter :: t0       = 15

  !--> perhaps better scaling in THCM if this would actually be used...
  real, parameter :: deltat   = 1.0
  real, parameter :: deltas   = 1.0
  
  real, parameter :: s0       = 35.0
  real, parameter :: cp0      = 4.2e+03
  ! Parameters equation of state
  ! b1 = 5.6e-05, b2 = 6.3e-06, b3 = 3.7e-08
  real :: alphaT   = 1.0e-04
  real :: alphaS   = 7.6e-04
  real, parameter :: alpt1    = 2.93    ! 2*b2*t0/b1 - 3*b3*t0*t0/b1
  real, parameter :: alpt2    = 8.3e-02 ! b2/b1 - 3*t0*b3/b1
  real, parameter :: alpt3    = 6.6e-04 ! b3/b1
  ! Mixing parameters
  ! real, parameter :: ah       = 2.2e+12 !2DMOC case (zonally averaged)
  real, parameter :: ah       = 2.5e+05 !2 deg resolution
  ! real, parameter :: ah       = 1.0e+05 !1 deg resolution
  ! real, parameter :: ah       = 1.0e+04 !0.5 deg resolution
  real, parameter :: av       = 1.0e-03
  real, parameter :: kappah   = 1.0e+03
  real, parameter :: kappav   = 1.0e-04
  !--obsolete---
  ! Extra parameters for vertical profiles of k_h and k_v (England 1993)
  !real, parameter :: Ar       = 1.05e-04/pi
  !real, parameter :: Alam     = 4.5e-03
  !real, parameter :: zstar    = -2500.
  !real, parameter :: Ash      = 1.0e+03
  ! Bottom values
  !real, parameter :: bottem   = 0.0
  !real, parameter :: botsal   = 0.0

  !________________________________________________________________


  !***************************************************************************

contains

  subroutine allocate_usr(dim_n,dim_m,dim_l)

    use m_par

    implicit none

    integer dim_n,dim_m,dim_l

    ! set dimensions in m_par:
    m=dim_m
    n=dim_n
    l=dim_l
    ndim = m*n*l*nun

    allocate(x(n),y(0:m+1),z(l),xu(0:n),yv(0:m),zw(0:l),ze(l),zwe(l),&
         dfzT(l),dfzW(0:l))

    allocate(landm(0:n+1,0:m+1,0:l+1))
    landm=OCEAN;! in case no topology is read in

    allocate(Frc(ndim), taux(n,m), tauy(n,m), tx(n,m), ty(n,m))
    allocate(ft(n,m), fs(n,m), tatm(n,m), emip(n,m), spert(n,m))
    allocate(qatm(n,m), albe(n,m), patm(n,m), adapted_emip(n,m))
    allocate(msi(n,m), gsi(n,m), qsa(n,m))
    allocate(ftlev(n,m,l), fslev(n,m,l))
    allocate(internal_temp(n,m,l), internal_salt(n,m,l))
    taux   = 0.0
    tauy   = 0.0
    tatm   = 0.0
    emip   = 0.0
    spert  = 0.0
    qatm   = 0.0
    albe   = 0.0
    patm   = 0.0
    msi    = 0.0
    gsi    = 0.0
    qsa    = 0.0    
    internal_temp  = 0.0
    internal_salt  = 0.0
    adapted_emip   = 0.0

  end subroutine allocate_usr

  subroutine deallocate_usr()

    use m_par

    implicit none

    deallocate(x,y,z,xu,yv,zw,ze,zwe,&
         dfzT,dfzW)

    deallocate(landm)

    deallocate(Frc, taux, tauy, tx, ty)
    deallocate(ft, fs, tatm, emip, spert)
    deallocate(qatm, albe, patm, adapted_emip)
    deallocate(msi, gsi, qsa)
    deallocate(ftlev, fslev)
    deallocate(internal_temp, internal_salt)
  end subroutine deallocate_usr

  !! can be called from C++ to get the grid arrays:
  !! xx/yy/zz are the cell centers and range from 1:n/m/l
  !! xxu/yyz/zzw are 0-based in fortran but one-based for
  !! this subroutine since referencing C-arrays as zero-
  !! based is not a good idea (I think).
  subroutine get_grid_data(xx,yy,zz,xxu,yyv,zzw)

    implicit none

    real, dimension(n) :: xx
    real, dimension(m) :: yy
    real, dimension(l) :: zz

    real, dimension(n+1) :: xxu
    real, dimension(m+1) :: yyv
    real, dimension(l+1) :: zzw

    xx(1:n) = x(1:n)
    yy(1:m) = y(1:m)
    zz(1:l) = z(1:l)

    xxu(1:n+1) = xu(0:n)
    yyv(1:m+1) = yv(0:m)
    zzw(1:l+1) = zw(0:l)

  end subroutine get_grid_data

  subroutine set_internal_forcing(ctemp,csalt)

    implicit none

    real, dimension(m*n*l), intent(in) :: ctemp,csalt

    integer :: i,j,k,pos


    if (.not. allocated(internal_temp)) then
       stop 'internal_temp not allocated'
    end if

    if (.not. allocated(internal_salt)) then
       stop 'internal_salt not allocated'
    end if
    OPEN(51,FILE='settemp.txt',STATUS='unknown')
    REWIND(51)

    pos = 1
    do k=1,l
       do j=1,m
          do i=1,n
             internal_temp(i,j,k) = ctemp(pos)
             internal_salt(i,j,k) = csalt(pos)
             WRITE(51,*) internal_temp(i,j,k)
             pos = pos + 1
          end do
       end do
    end do
    CLOSE(51)
  end subroutine set_internal_forcing

end module m_usr
