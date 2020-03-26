#include "fdefs.h"

!**************************************************************************
!! initialize THCM: input: number of grid-points in x,y and z-direction,
!! bounds of the domain (formerly set in usr.com)
SUBROUTINE init(a_n,a_m,a_l,a_nmlglob,&
     a_xmin,a_xmax,a_ymin,a_ymax,&
     a_alphaT,a_alphaS,&
     a_ih,a_vmix,a_tap,a_rho_mixing,&
     a_coriolis_on,&
     a_periodic,a_landm,&
     a_taux,a_tauy,a_tatm,a_emip,a_spert)

  use, intrinsic :: iso_c_binding
  use m_usr
  use m_res
  use m_mix
  use m_atm
  use m_mat

  implicit none

  integer(c_int) :: a_n,a_m,a_l,a_nmlglob
  real(c_double) :: a_xmin,a_xmax,a_ymin,a_ymax
  real(c_double) :: a_alphaT, a_alphaS
  integer(c_int) :: a_ih, a_vmix, a_tap, a_rho_mixing
  integer(c_int) :: a_coriolis_on
  integer(c_int) :: a_periodic
  integer(c_int), dimension((a_n+2)*(a_m+2)*(a_l+2)) :: a_landm
  real(c_double), dimension(a_n*a_m) :: a_taux,a_tauy
  real(c_double), dimension(a_n*a_m) :: a_tatm,a_emip,a_spert

  ! LOCAL
  real    :: dzne
  integer :: i,j,k,pos

  _INFO_('THCM: init... ')

  nmlglob = a_nmlglob

  xmin    = a_xmin
  xmax    = a_xmax
  ymin    = a_ymin
  ymax    = a_ymax

  alphaT  = a_alphaT
  alphaS  = a_alphaS

  ih      = a_ih
  vmix    = a_vmix
  tap     = a_tap

  if (a_rho_mixing.ne.0) then
     rho_mixing = .true.
  else
     rho_mixing = .false.
  end if

  coriolis_on = a_coriolis_on

  !initialize atmos coefficients
  qdim = 0.01
  nuq  = 0.0
  eta  = 0.0
  dqso = 0.0

  if (a_periodic.eq.0) then
     periodic  =  .false.
  else
     periodic  =  .true.
  end if

  call allocate_usr(a_n,a_m,a_l)
  call allocate_mat()
  call allocate_atm(m)
  call allocate_res(ndim)
  call allocate_mix()

  ! Fill the local landmask array. Replace all PERIO
  ! by OCEAN if periodic is false, so that THCM does
  ! not know about global periodicity (which would
  ! mess up the matrix partitioning).
  pos = 1
  do k = 0, l+1
     do j = 0, m+1
        do i = 0, n+1
           landm(i,j,k) = a_landm(pos)
           if( (.not.periodic) .and. (landm(i,j,k).eq.PERIO)) then
              landm(i,j,k) = OCEAN
           end if
           pos = pos+1
        end do
     end do
  end do

  ! Make the dummy cells all land so the indexing
  ! scheme (assemble.F90) doesn't get confused.
  ! This means that there have to be at least two
  ! layers of overlap in the global parallel session.
  if (.not. periodic) then
     landm(0,:,:)   = LAND
     landm(n+1,:,:) = LAND
  end if
  landm(:,0,:)      = LAND
  landm(:,m+1,:)    = LAND
  landm(:,:,0)      = LAND
  landm(:,:,l+1)    = LAND

  pos = 1
  do j=1,m
     do i=1,n
        taux(i,j)  = a_taux(pos)
        tauy(i,j)  = a_tauy(pos)
        tatm(i,j)  = a_tatm(pos)
        emip(i,j)  = a_emip(pos)
        spert(i,j) = a_spert(pos)
        pos=pos+1
     end do
  end do

  call grid

  ! When the grid is known we can set nondimensionalization
  ! coefficients for the body forcing.
  dzne = dz*dfzT(l)
  QTnd = r0dim / (udim * cp0 * rhodim * hdim * dzne )
  QSnd = s0 * r0dim / ( deltas * udim * hdim * dzne )

!  write(*,*) 'THCM: nondim constants: QTnd = ', QTnd, &
!       ' QSnd = ', QSnd

  call stpnt        !
  call vmix_init    ! ATvS-Mix  USES LANDMASK
  call atmos_coef   !
  call forcing      ! USES LANDMASK
  call lin

  _INFO_('THCM: init...  done')
end subroutine init

!****************************************************************************
! deallocate all dynamically alloc'd memory
subroutine finalize

  use m_usr
  use m_res
  use m_mix
  use m_atm
  use m_mat
  implicit none

  call deallocate_usr()
  call deallocate_mat()
  call deallocate_atm()
  call deallocate_res()
  call deallocate_mix()

  close(f99)

end subroutine finalize

!*****************************************************************************
SUBROUTINE setparcs(param,value)
  !     interface for Trilinos to set the thirty continuation variables
  use, intrinsic :: iso_c_binding
  use m_usr
  implicit none
  integer(c_int) param
  real(c_double) value
  !WRITE(f99,*) 'setting par(',param,')=',value
  IF ((param>=1).AND.(param<=npar)) THEN
     PAR(param) = value
  ELSE
     WRITE(f99,*) 'error in transfer parameter to fortran'
  ENDIF
  !     ENDIF

  call forcing
  call lin

END SUBROUTINE setparcs

!*****************************************************************************
SUBROUTINE getparcs(param,value)
  !     interface for Trilinos to get the thirty continuation variables
  use, intrinsic :: iso_c_binding
  use m_usr
  implicit none
  integer(c_int) param
  real(c_double) value
  !WRITE(f99,*) 'setting par(',param,')=',value
  IF ((param>=1).AND.(param<=npar)) THEN
     value=PAR(param)
  ELSE
     WRITE(f99,*) 'error in transfer parameter from fortran'
  ENDIF
  !     ENDIF
end subroutine getparcs

!***********************************************************
SUBROUTINE getdeps(o_Ooa, o_Os, o_nus, o_eta, o_lvsc, o_qdim, o_pqsnd)
  !     interface to get Ooa and other dependencies on external model
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_atm
  implicit none
  real(c_double) o_Ooa, o_Os, o_nus, o_eta, o_lvsc, o_qdim, o_pqsnd
  real pQSnd

  pQSnd = par(COMB) * par(SALT) * QSnd

  o_Ooa   = Ooa
  o_Os    = Os
  o_nus   = nus
  o_eta   = eta
  o_lvsc  = lvsc
  o_qdim  = qdim
  o_pqsnd = pQSnd
end subroutine getdeps

!**********************************************************
SUBROUTINE get_parameters(o_r0dim, o_udim, o_hdim)
  !     interface to get a few model constants
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_atm
  implicit none
  real(c_double) o_r0dim, o_udim, o_hdim

  o_r0dim = r0dim;
  o_udim  = udim;
  o_hdim  = hdim;

end subroutine get_parameters

!**********************************************************
SUBROUTINE set_atmos_parameters(pars)
  ! Interface to set a few model parameters relevant for E-P. These
  ! parameters affect the sensitivity nus, which should be updated
  ! here.
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_atm
  implicit none

  ! This should be equal to the CommPars struct in AtmosLocal.H
  ! declarations.
  type, bind(C) :: atmos_pars
     real(c_double) :: tdim
     real(c_double) :: qdim
     real(c_double) :: nuq
     real(c_double) :: eta
     real(c_double) :: dqso
     real(c_double) :: dqsi
     real(c_double) :: dqdt
     real(c_double) :: Eo0
     real(c_double) :: Ei0
     real(c_double) :: Cs
     real(c_double) :: t0o
     real(c_double) :: t0i
     real(c_double) :: a0
     real(c_double) :: da
     real(c_double) :: tauf
     real(c_double) :: tauc
     real(c_double) :: comb
     real(c_double) :: albf
  end type atmos_pars
  
  type(atmos_pars) :: pars
  real :: dzne

  qdim  = pars%qdim
  nuq   = pars%nuq
  eta   = pars%eta
  dqso  = pars%dqso
  eo0   = pars%Eo0
  albe0 = pars%a0
  albed = pars%da

  dzne = dz*dfzT(l)

  nus  = par(COMB) * par(SALT) * eta * qdim * QSnd

  !   write(*,*) 'THCM: nondim constants: nus = ', nus

  ! --> FIXME Instead of using par(TEMP) we should have a dedicated latent heat
  ! continuation parameter.
  lvsc = par(COMB) * par(TEMP) * rhodim * lv * QTnd

  call forcing
  call lin

end subroutine set_atmos_parameters

!**********************************************************
SUBROUTINE set_seaice_parameters(pars)
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_ice
  implicit none

  ! This should be equal to the CommPars struct in SeaIce.H
  ! declarations.
  type, bind(C) :: seaice_pars
     real(c_double) :: zeta
     real(c_double) :: a0
     real(c_double) :: Lf
     real(c_double) :: s0
     real(c_double) :: rhoo
     real(c_double) :: qvar
     real(c_double) :: q0
  end type seaice_pars

  type(seaice_pars) :: pars
       
  zeta = pars%zeta ! combination of sea ice parameters
  a0   = pars%a0   ! freezing temperature S sensitivity
  Lf   = pars%Lf   ! latent heat of fusion of ice
  Qvar = pars%qvar ! typical QTsa heat flux variation
  q0   = pars%q0   ! background QTsa

  if (pars%s0.ne.s0) then
     _INFO_('WARNING conflicting reference salinity s0')
  endif

  if (pars%rhoo.ne.rhodim) then
     _INFO_('WARNING conflicting sea water density rhodim')
  endif

  call forcing
  call lin

end subroutine set_seaice_parameters

!***********************************************************
SUBROUTINE set_landmask(a_landm, a_periodic, a_reinit)
  ! interface to set a new landmask
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_mix

  implicit none

  integer(c_int) :: a_periodic
  integer(c_int) :: a_reinit
  integer(c_int), dimension((n+2)*(m+2)*(l+2)) :: a_landm

  integer :: i,j,k,pos

!  _INFO_('THCM: usrc.F90 set_landmask...')

  if (a_periodic.eq.0) then
     periodic  =  .false.
  else
     periodic  =  .true.
  end if

  ! Fill the local landmask array
  pos = 1
  do k = 0, l+1
     do j = 0, m+1
        do i = 0, n+1
           landm(i,j,k) = a_landm(pos)
           if( (.not.periodic) .and. (landm(i,j,k).eq.PERIO)) then
              landm(i,j,k) = OCEAN
           end if
           pos = pos+1
        end do
     end do
  end do

  do i = 1, n
     do j = 1, m
        do k = l, 2,-1
           if ( landm(i,j,k).eq.LAND .and. landm(i,j,k-1).eq.OCEAN ) then
              write(*,*) 'land inversion at ',i,j,k
              landm(i,j,k-1) = LAND
           endif
        enddo
     enddo
  enddo

  ! Let the dummy cells be land
  if (.not.periodic) then
     landm(0,:,:)   = LAND
     landm(n+1,:,:) = LAND
  end if
  landm(:,0,:)    = LAND
  landm(:,m+1,:)  = LAND
  landm(:,:,0)    = LAND
  landm(:,:,l+1)  = LAND

  if (a_reinit.eq.1) then
     !  A few initializations need to be repeated
     call vmix_init    ! ATvS-Mix  USES LANDMASK
     call forcing      ! USES LANDMASK
     call lin
  endif

!  _INFO_('THCM: usrc.F90 set_landmask...  done')
end subroutine set_landmask

!*****************************************************************************
SUBROUTINE writeparams()
  ! write the entire parameter list (1-30) to fort.7
  use, intrinsic :: iso_c_binding
  use m_usr
  implicit none
  integer i
  write(*,*) 'Writing to fort.7'
  write(7, *) '---------------------------'
  write(7, '(5e15.5)') (PAR(i), i=1, npar)
  write(7, *) '---------------------------'
end subroutine writeparams

!*****************************************************************************
SUBROUTINE setsres(tmp_sres)
  ! adjust sres parameter
  use, intrinsic :: iso_c_binding
  use m_usr
  implicit none
  integer(c_int) :: tmp_sres

  sres = tmp_sres

  call forcing
  call lin

end SUBROUTINE setsres

!*****************************************************************************
SUBROUTINE matrix(un)
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_mix
  use m_atm
  USE m_mat
  USE m_res
  !     construct the jacobian A and the 'mass' matrix B
  !     Put B in coB,  A-sig*B in coA
  !     sig1: u/v/T/S
  !     sig2: w/p
  implicit none
  real(c_double),dimension(ndim) :: un
  real time0, time1

  An = Al

  _DEBUG_("Build diagonal matrix B...")
  call fillcolB
#ifndef THCM_LINEAR
  _DEBUG_("Build nonlinear part of Jacobian...")
  call nlin_jac(un)
#endif
  !{ removing tons of things for eigen-analysis test...
#if 0
  !Euv
  An(UU:VV,WW,:,:,:) = 0.0
  !BTS
  An(WW,TT:SS,:,:,:) = 0.0
  !Buv, Bw
  !An(TT:SS,UU::WW,:,:,:) = 0.0
  !Gw
  !An(:,WW,PP,:,:,:) = 0.0
  !Dw
  !An(:,PP,WW,:,:,:) = 0.0
#endif
  !}


  ! ATvS-Mix ---------------------------------------------------------------------
  if (vmix_flag.ge.1) then
     call TIMER_START('mixing jac' // char(0))
     call cpu_time(time0)
     if (vmix_out.gt.0) write(99,'(a26)')"MIX| matrix...    "
     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     fix:   ', vmix_fix
     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     flag:  ', vmix_flag

     if ((vmix_fix.eq.0).and.(vmix_flag.ge.2)) call vmix_control(un)

     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     temp:  ', vmix_temp
     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     salt:  ', vmix_salt

     if (((vmix_temp.eq.1).or.(vmix_salt.eq.1)).and.(vmix_dim.gt.0)) then
        call vmix_jac(un)
     endif
     call cpu_time(time1)
     vmix_time=vmix_time+time1-time0
     if (vmix_out.gt.0) write (99,'(a26, f10.3)') 'MIX|        ...matrix done', time1-time0
     call TIMER_STOP('mixing jac' // char(0))
  endif
  ! --------------------------------------------------------------------- ATvS-Mix

  call boundaries

  call assemble
  !write(*,*) "In fortran's matrix() maxval TT =", maxval(An(:,TT,:,:,:,:))
  !write(*,*) "In fortran's matrix() maxval SS =", maxval(An(:,SS,:,:,:,:))

  ! call writematrhs(0.0)

  ! call writecsrmats
  ! stop
end SUBROUTINE matrix
!****************************************************************************
SUBROUTINE rhs(un,B)
  !     construct the right hand side B
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_mix
  use m_res
  use m_mat

  implicit none
  real(c_double),dimension(ndim) ::    un,B
  real    mix(ndim) ! ATvS-Mix
  real    Au(ndim), time0, time1
  integer i,j,k,k1,row,find_row2, mode

  !call writeparameters
  mix = 0.0
  An = Al
  ! write(*,*) 'T(n,m,l)', un(find_row2(n,m,l,TT))
#ifndef THCM_LINEAR
  call nlin_rhs(un)
#endif
  ! call forcing          !
  call boundaries       !
  call assemble
  call TIMER_START('matAvec' // char(0))
  call matAvec(un,Au)   !
  call TIMER_STOP('matAvec' // char(0))
  ! ATvS-Mix ---------------------------------------------------------------------
  if (vmix_flag.ge.1) then
     call TIMER_START('mixing rhs' // char(0))
     mode=vmix_fix
     call cpu_time(time0)
     if (vmix_out.gt.0) write(99,'(a26)')'MIX| rhs...               '
     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     fix:   ', vmix_fix

     if ((vmix_fix.eq.0).and.(vmix_flag.ge.2)) call vmix_control(un)

     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     temp:  ', vmix_temp
     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     salt:  ', vmix_salt

     if (((vmix_temp.eq.1).or.(vmix_salt.eq.1)).and.(vmix_dim.gt.0)) then
        call vmix_fun(un, mix)
     endif
     call cpu_time(time1)
     vmix_time=vmix_time+time1-time0
     if (vmix_out.gt.0) write (99,'(a26,f10.3)') 'MIX|    ...rhs done',time1-time0
     call TIMER_STOP('mixing rhs' // char(0))
  endif
  ! --------------------------------------------------------------------- ATvS-Mix

  _DEBUG2_("p0 = ", p0) ! Residue Continuation

  call TIMER_START('addition rhs' // char(0))
  B = -Au - mix + Frc - p0*(1- par(RESC))*ures
  call TIMER_STOP('addition rhs' // char(0))
#if 1
  call TIMER_START('landmask rhs' // char(0))
  if(ires == 0) then
     DO i = 1, n
        DO j = 1, m
           DO k = 1, l
              DO k1 = 1,nun
                 row = find_row2(i,j,k,k1)
                 B(row) = B(row) * (1 - landm(i,j,k))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  endif
  call TIMER_STOP('landmask rhs' // char(0))
#endif

  ! open(15,file='Frc.co')
  ! do i=1,ndim
  !    write(15,*) Frc(i)
  ! enddo
  ! close(15)

  _DEBUG2_("maxval rhs= ", maxval(abs(B)))

end SUBROUTINE rhs
!****************************************************************************
SUBROUTINE lin
  USE m_mat
  !     Thermohaline equations
  !     Produce local element matrices for linear operators
  ! +---------------------------------------------------------------------+
  ! |   Al is a list of dependencies between the unknowns, see m_mat      |
  ! |      Al(loc,A,B,i,j,k) = c                                          |
  ! |       where loc is one of the locations below:                      |
  ! |     +----------++-------++----------+                               |
  ! |     | 12 15 18 || 3 6 9 || 21 24 27 |                               |
  ! |     | 11 14 17 || 2 5 8 || 20 23 26 |                               |
  ! |     | 10 13 16 || 1 4 7 || 19 22 25 |                               |
  ! |     |  below   || center||  above   |                               |
  ! |     +----------++-------++----------+                               |
  ! |                                                                     |
  ! |     For instance, Al(14,A,B,i,j,k) = c is                           |
  ! |     d/dt A|(i,j,k) = c*B|(i,j,k-1) + ...                            |
  ! +---------------------------------------------------------------------+
  use m_usr
  use m_atm
  use m_ice
  implicit none
  !     LOCAL
  real,target :: ucsi(np,n,m,l),&
       &         uxx(np,n,m,l),uyy(np,n,m,l),uzz(np,n,m,l),&
       &         uxs(np,n,m,l),fu(np,n,m,l),px(np,n,m,l)
  real    sc(np,n,m,l),tcb(np,n,m,l)
  real    EH,EV,ph,pv,Ra,lambda, bi, dedt
  real    xes

  real    mc(np,n,m,l)
  real    QSoa(np,n,m,l), QSos(np,n,m,l)
  real    pQSnd

  ! original version:
  !      equivalence (u, v), (uy, vy), (ucsi, vcsi), (uxx, vxx, txx, uxc)
  !      equivalence (uyy, vyy, tyy, vyc), (uzz, vzz, tzz, wzc)
  !      equivalence (uxs, vxs, tbc, tyc)
  !      equivalence (fu, fv, tc), (px, py, pz)

  ! we use pointers instead:
  real,dimension(:,:,:,:),pointer :: vcsi,vxx,txx,uxc,&
       &     vyy,tyy,vyc,vzz,tzz,wzc,vxs,tbc,tyc,fv,tc,py,pz

  vcsi=>ucsi
  vxx=>uxx
  txx=>uxx
  uxc=>uxx
  vyy=>uyy
  tyy=>uyy
  vyc=>uyy
  vzz=>uzz
  tzz=>uzz
  wzc=>uzz
  vxs=>uxs
  tbc=>uxs
  tyc=>uxs
  fv=>fu
  tc=>fu
  py=>px
  pz=>px

  EV     = par(EK_V)
  EH     = par(EK_H)
  ph     = (1-par(MIXP))*par(PE_H)
  pv     = par(PE_V)
  lambda = par(LAMB)
  xes    = par(NLES)
  bi     = par(BIOT)
  Ra     = par(RAYL)
  !      rintt  = par(IFRICT)   ! ATvS-Mix

  Al = 0.0

  ! ------------------------------------------------------------------
  ! u-equation
  ! ------------------------------------------------------------------
  ! call uderiv(1,u)
  call uderiv(2,uxx)
  call uderiv(3,uyy)
  call uderiv(4,uzz)
  call uderiv(5,ucsi)
  call uderiv(6,vxs)
  call coriolis(1,fv)
  call gradp(1,px)
  Al(:,UU,UU,:,:,1:l) = -EH * (uxx+uyy+ucsi) -EV * uzz ! + rintt*u ! ATvS-Mix
  Al(:,UU,VV,:,:,1:l) = -fv - EH*vxs
  ! Al(:,UU,VV,:,:,1:l) = - EH*vxs ! for 2DMOC case
  Al(:,UU,PP,:,:,1:l) =  px

  ! ------------------------------------------------------------------
  ! v-equation
  ! ------------------------------------------------------------------
  ! call vderiv(1,v)
  call vderiv(2,vxx)
  call vderiv(3,vyy)
  call vderiv(4,vzz)
  call vderiv(5,vcsi)
  call vderiv(6,uxs)
  call coriolis(2,fu)
  call gradp(2,py)
  Al(:,VV,UU,:,:,1:l) =  fu - EH*uxs
  ! Al(:,VV,UU,:,:,1:l) =  - EH*uxs ! for 2dMOC case
  Al(:,VV,VV,:,:,1:l) = -EH*(vxx + vyy + vcsi) - EV*vzz !+ rintt*v ! ATvS-Mix
  Al(:,VV,PP,:,:,1:l) =  py

  ! ------------------------------------------------------------------
  ! w-equation
  ! ------------------------------------------------------------------
  call gradp(3,pz)
  call tderiv(6,tbc)
  Al(:,WW,PP,:,:,1:l) =  pz
  Al(:,WW,TT,:,:,1:l) = -Ra *(1. + xes*alpt1) * tbc/2.
  Al(:,WW,SS,:,:,1:l) =  lambda * Ra * tbc/2.

  ! ------------------------------------------------------------------
  ! p-equation
  ! ------------------------------------------------------------------
  call pderiv(1,uxc)
  call pderiv(2,vyc)
  call pderiv(3,wzc)
  Al(:,PP,UU,:,:,1:l) = uxc
  Al(:,PP,VV,:,:,1:l) = vyc
  Al(:,PP,WW,:,:,1:l) = wzc

  ! ------------------------------------------------------------------
  ! T-equation
  ! ------------------------------------------------------------------
  call tderiv(1,tc )
  call tderiv(2,sc )
  call tderiv(3,txx)
  call tderiv(4,tyy)
  call tderiv(5,tzz)
  call tderiv(7,tcb)

  call masksi(mc, msi); ! create sea ice mask atom

  ! dependence of TT on TT through latent heat due to evaporation
  dedt =  lvsc * eta * qdim * (deltat / qdim) * dqso
  
  if (coupled_T.eq.1) then ! coupled with external atmos
     ! FIXME is this too much mc*tc? TEM

     Al(:,TT,TT,:,:,1:l) =                &
          - ph * (txx + tyy) - pv * tzz   & ! diffusive transport
          + Ooa  * tc                     & ! sensible heat flux
          + dedt * sc                     & ! latent heat flux
          + mc * (QTnd * zeta * tc  -     & ! correction for sea ice
          Ooa * tc - dedt * sc)

     Al(:,TT,SS,:,:,1:l) = -QTnd * zeta * a0 * mc  ! salinity dependence in
                                                   ! freezing temperature
  else
     Al(:,TT,TT,:,:,1:l) = -ph * (txx + tyy) - pv * tzz + TRES*bi*tc
  endif

  ! ------------------------------------------------------------------
  ! S-equation
  ! ------------------------------------------------------------------
  ! dependence of SS on TT in evaporation term
  ! FIXME: naming is bad and confusing
  dedt  = nus * (deltat / qdim) * dqso
  pQSnd = par(COMB) * par(SALT) * QSnd

  ! FIXME: ugly
  if (coupled_S.eq.1) then ! coupled to atmosphere
     Al(:,SS,SS,:,:,1:l) = - ph * (txx + tyy) - pv * tzz &
          - mc * pQSnd * zeta * a0 / (rhodim * Lf)

     ! minus sign and nondim added (we take -Au in rhs computation)

     QSoa = -dedt * sc            ! atmosphere to ocean salinity flux
                                  ! internal component

     QSos =  pQSnd * zeta / (rhodim * Lf)    ! sea ice to ocean salinity flux
                                             ! internal component

     ! combine contributions with mask
     Al(:,SS,TT,:,:,1:l) =   QSoa + mc * (QSos - QSoa)
  else
     Al(:,SS,SS,:,:,1:l) = - ph * (txx + tyy) - pv * tzz + SRES*bi*sc
  endif


end SUBROUTINE lin

!********************************************************************
SUBROUTINE nlin_rhs(un)
  use, intrinsic :: iso_c_binding
  USE m_mat
  !     Produce local matrices for nonlinear operators for calc of Rhs
  use m_usr
  use m_mix
  use m_atm
  implicit none

  !     IMPORT/EXPORT
  real(c_double),dimension(ndim) ::    un

  !     LOCAL
  real    u(0:n  ,0:m,0:l+1), v(0:n,0:m  ,0:l+1)
  real    w(0:n+1,0:m+1,0:l  ), p(0:n+1,0:m+1,0:l+1)
  real    t(0:n+1,0:m+1,0:l+1), s(0:n+1,0:m+1,0:l+1)
  real,target ::    utx(np,n,m,l),vty(np,n,m,l),wtz(np,n,m,l)
  real    t2r(np,n,m,l),t3r(np,n,m,l)

  real    uux(np,n,m,l),uvy1(np,n,m,l),uwz(np,n,m,l),uvy2(np,n,m,l)
  real    uvx(np,n,m,l),vvy(np,n,m,l),vwz(np,n,m,l),ut2(np,n,m,l)

  real    lambda,epsr,Ra,xes

  real,dimension(:,:,:,:),pointer ::    usx,vsy,wsz

  call TIMER_START('nlin_rhs' // char(0))
  usx=>utx
  vsy=>vty
  wsz=>wtz

  lambda = par(LAMB)
  epsr   = par(ROSB)
  Ra     = par(RAYL)
  xes    = par(NLES)
  ! pvc1   = par(P_VC)
  ! pv     = par(PE_V)
  ! pvc2   = pv*(1.0 - par(ALPC))*par(ENER)
  call usol(un,u,v,w,p,t,s)
  ! rho    = lambda*s - t *( 1 + xes*alpt1) - &
  !          xes*t*t*alpt2+xes*t*t*t*alpt3

  ! ------------------------------------------------------------------
  ! u-equation
  ! ------------------------------------------------------------------
#ifndef NO_UVNLIN
  call unlin(1,uux,u,v,w)
  call unlin(3,uvy1,u,v,w)
  call unlin(5,uwz,u,v,w)
  call unlin(7,uvy2,u,v,w)
  An(:,UU,UU,:,:,1:l) = An(:,UU,UU,:,:,1:l) + epsr * (uux + uvy1 + uwz + uvy2)
#endif

  ! ------------------------------------------------------------------
  ! v-equation
  ! ------------------------------------------------------------------
#ifndef NO_UVNLIN
  call vnlin(1,uvx,u,v,w)
  call vnlin(3,vvy,u,v,w)
  call vnlin(5,vwz,u,v,w)
  call vnlin(7,ut2,u,v,w)
  An(:,VV,UU,:,:,1:l) = An(:,VV,UU,:,:,1:l) + epsr *ut2
  An(:,VV,VV,:,:,1:l) = An(:,VV,VV,:,:,1:l) + epsr*(uvx + vvy + vwz)
#endif

  ! ------------------------------------------------------------------
  ! w-equation
  ! ------------------------------------------------------------------
  call wnlin(2,t2r,t)
  call wnlin(4,t3r,t)
  An(:,WW,TT,:,:,1:l) = An(:,WW,TT,:,:,1:l) - Ra*xes*alpt2*t2r &
                                            + Ra*xes*alpt3*t3r

  ! ------------------------------------------------------------------
  ! T-equation
  ! ------------------------------------------------------------------
#ifndef NO_TSNLIN
  call tnlin(3,utx,u,v,w,t)
  call tnlin(5,vty,u,v,w,t)
  call tnlin(7,wtz,u,v,w,t)
  An(:,TT,TT,:,:,1:l) = An(:,TT,TT,:,:,1:l)+ utx+vty+wtz        ! ATvS-Mix
#endif

  ! ------------------------------------------------------------------
  ! S-equation
  ! ------------------------------------------------------------------
#ifndef NO_TSNLIN
  call tnlin(3,usx,u,v,w,s)
  call tnlin(5,vsy,u,v,w,s)
  call tnlin(7,wsz,u,v,w,s)
  An(:,SS,SS,:,:,1:l) = An(:,SS,SS,:,:,1:l)+ usx+vsy+wsz        ! ATvS-Mix
#endif

  call TIMER_STOP('nlin_rhs' // char(0))

end SUBROUTINE nlin_rhs

!****************************************************************************
SUBROUTINE nlin_jac(un)
  !     Produce local matrices for nonlinear operators for calc of Jacobian
  use, intrinsic :: iso_c_binding
  USE m_mat
  use m_usr
  use m_mix
  use m_atm
  implicit none

  !     IMPORT/EXPORT
  real(c_double), dimension(ndim) ::   un

  !     LOCAL
  real    u(0:n  ,0:m,0:l+1), v(0:n, 0:m  ,0:l+1)
  real    w(0:n+1,0:m+1,0:l  ), p(0:n+1,0:m+1,0:l+1)
  real    t(0:n+1,0:m+1,0:l+1), s(0:n+1,0:m+1,0:l+1)
  real,target :: urTx(np,n,m,l),Utrx(np,n,m,l),&
       &        vrTy(np,n,m,l),Vtry(np,n,m,l)
  real,target :: wrTz(np,n,m,l),Wtrz(np,n,m,l)
  real    t2r(np,n,m,l),t3r(np,n,m,l)

  real    lambda,epsr,Ra,xes
  real uvy1(np,n,m,l),uwz(np,n,m,l),uvy2(np,n,m,l)
  real uvx(np,n,m,l),vwz(np,n,m,l)
  real Urux(np,n,m,l),Urvy1(np,n,m,l),Urwz(np,n,m,l),Urvy2(np,n,m,l)
  real uVrx(np,n,m,l),Vrvy(np,n,m,l),Vrwz(np,n,m,l),Urt2(np,n,m,l)

  real,dimension(:,:,:,:),pointer :: urSx,Usrx,vrSy,Vsry,wrSz,Wsrz

  call TIMER_START('nlin_jac' // char(0))

  urSx=>urTx
  vrSy=>vrTy
  wrSz=>wrTz
  Usrx=>Utrx
  Vsry=>Vtry
  Wsrz=>Wtrz

  lambda = par(LAMB)
  epsr   = par(ROSB)
  Ra     = par(RAYL)
  xes    = par(NLES)
  ! pvc1   = par(P_VC)
  ! pv     = par(PE_V)
  ! pvc2   = pv*(1.0 - par(ALPC))*par(ENER)
  call usol(un,u,v,w,p,t,s)
  ! rho    = lambda*s - t *( 1 + xes*alpt1) - &
  !          xes*t*t*alpt2+xes*t*t*t*alpt3

  ! ------------------------------------------------------------------
  ! u-equation
  ! ------------------------------------------------------------------
#ifndef NO_UVNLIN
  call unlin(2,Urux,u,v,w)
  call unlin(3,uvy1,u,v,w)
  call unlin(4,Urvy1,u,v,w)
  call unlin(5,uwz,u,v,w)
  call unlin(6,Urwz,u,v,w)
  call unlin(7,uvy2,u,v,w)
  call unlin(8,Urvy2,u,v,w)
  An(:,UU,UU,:,:,1:l)  =  An(:,UU,UU,:,:,1:l) + epsr * (Urux + uvy1 + uwz + uvy2)
  An(:,UU,VV,:,:,1:l)  =  An(:,UU,VV,:,:,1:l) + epsr * (Urvy1 + Urvy2)
  An(:,UU,WW,:,:,1:l)  =  An(:,UU,WW,:,:,1:l) + epsr *  Urwz
#endif

  ! ------------------------------------------------------------------
  ! v-equation
  ! ------------------------------------------------------------------
#ifndef NO_UVNLIN
  call vnlin(1,uvx,u,v,w)
  call vnlin(2,uVrx,u,v,w)
  call vnlin(4,Vrvy,u,v,w)
  call vnlin(5,vwz,u,v,w)
  call vnlin(6,Vrwz,u,v,w)
  call vnlin(8,Urt2,u,v,w)
  An(:,VV,UU,:,:,1:l) =   An(:,VV,UU,:,:,1:l) + epsr * (Urt2 + uVrx)
  An(:,VV,VV,:,:,1:l) =   An(:,VV,VV,:,:,1:l) + epsr * (uvx + Vrvy + vwz)
  An(:,VV,WW,:,:,1:l) =   An(:,VV,WW,:,:,1:l) + epsr * Vrwz
#endif

  ! ------------------------------------------------------------------
  ! w-equation
  ! ------------------------------------------------------------------
  call wnlin(1,t2r,t)
  call wnlin(3,t3r,t)
  An(:,WW,TT,:,:,1:l) = An(:,WW,TT,:,:,1:l) - Ra*xes*alpt2*t2r &
                                            + Ra*xes*alpt3*t3r

  ! ------------------------------------------------------------------
  ! T-equation
  ! ------------------------------------------------------------------
#ifndef NO_TSNLIN
  call tnlin(2,urTx,u,v,w,t)
  call tnlin(3,Utrx,u,v,w,t)
  call tnlin(4,vrTy,u,v,w,t)
  call tnlin(5,Vtry,u,v,w,t)
  call tnlin(6,wrTz,u,v,w,t)
  call tnlin(7,Wtrz,u,v,w,t)
  An(:,TT,UU,:,:,1:l) = An(:,TT,UU,:,:,1:l) + urTx
  An(:,TT,VV,:,:,1:l) = An(:,TT,VV,:,:,1:l) + vrTy
  An(:,TT,WW,:,:,1:l) = An(:,TT,WW,:,:,1:l) + wrTz
  An(:,TT,TT,:,:,1:l) = An(:,TT,TT,:,:,1:l) + Utrx + Vtry + Wtrz        ! ATvS-Mix
#endif

  ! ------------------------------------------------------------------
  ! S-equation
  ! ------------------------------------------------------------------
#ifndef NO_TSNLIN
  call tnlin(2,urSx,u,v,w,s)
  call tnlin(3,Usrx,u,v,w,s)
  call tnlin(4,vrSy,u,v,w,s)
  call tnlin(5,Vsry,u,v,w,s)
  call tnlin(6,wrSz,u,v,w,s)
  call tnlin(7,Wsrz,u,v,w,s)
  An(:,SS,UU,:,:,1:l) = An(:,SS,UU,:,:,1:l) + urSx
  An(:,SS,VV,:,:,1:l) = An(:,SS,VV,:,:,1:l) + vrSy
  An(:,SS,WW,:,:,1:l) = An(:,SS,WW,:,:,1:l) + wrSz
  An(:,SS,SS,:,:,1:l) = An(:,SS,SS,:,:,1:l) + Usrx + Vsry + Wsrz
#endif

  call TIMER_STOP('nlin_jac' // char(0))

end SUBROUTINE nlin_jac
!****************************************************************************
SUBROUTINE usol(un,u,v,w,p,t,s)
  !     Go from un to u,v,t,h
  use m_usr
  implicit none
  !     IMPORT/EXPORT
  real    un(ndim)
  real    u(0:n  ,0:m,0:l+1), v(0:n,0:m  ,0:l+1)
  real    w(0:n+1,0:m+1,0:l  ), p(0:n+1,0:m+1,0:l+1)
  real    t(0:n+1,0:m+1,0:l+1), s(0:n+1,0:m+1,0:l+1)
  !     LOCAL
  integer i,j,k
  !     EXTERNAL
  integer find_row2

  u = 0.0
  v = 0.0
  w = 0.0
  p = 0.0
  t = 0.0
  s = 0.0
  do k = 1, l
     do j = 1, m
        do i = 1, n
           u(i,j,k) = un(find_row2(i,j,k,UU))
           v(i,j,k) = un(find_row2(i,j,k,VV))
           w(i,j,k) = un(find_row2(i,j,k,WW))
           p(i,j,k) = un(find_row2(i,j,k,PP))
           t(i,j,k) = un(find_row2(i,j,k,TT))
           s(i,j,k) = un(find_row2(i,j,k,SS))
        enddo
     enddo
  enddo
  do k = 1, l
     do j = 1, m
        if (periodic) then
           u(0,j,k)  = u(N,j,k)
           v(0,j,k)  = v(N,j,k)
           w(N+1,j,k)= w(1,j,k)
           w(0,j,k)  = w(N,j,k)
           p(N+1,j,k)= p(1,j,k)
           p(0,j,k)  = p(N,j,k)
           t(N+1,j,k)= t(1,j,k)
           t(0,j,k)  = t(N,j,k)
           s(N+1,j,k)= s(1,j,k)
           s(0,j,k)  = s(N,j,k)
        else
           u(0,j,k)  = 0.0
           u(N,j,k)  = 0.0
           v(0,j,k)  = 0.0
           v(N,j,k)  = 0.0
           p(0,j,k)  = 0.0
           p(N+1,j,k)= 0.0
           t(0,j,k)  = t(1,j,k)
           t(N+1,j,k)= t(N,j,k)
           s(0,j,k)  = s(1,j,k)
           s(N+1,j,k)= s(N,j,k)
        endif
     enddo
  enddo
  do k = 1, l
     do i = 1, n
        u(i,0,k)  = 0.0
        u(i,M,k)  = 0.0
        v(i,0,k)  = 0.0
        v(i,M,k)  = 0.0
        p(i,0,k)  = 0.0
        p(i,M+1,k)= 0.0
        t(i,0,k)  = t(i,1,k)
        t(i,M+1,k)= t(i,M,k)
        s(i,0,k)  = s(i,1,k)
        s(i,M+1,k)= s(i,M,k)
     enddo
  enddo
  do j = 1, m
     do i = 1, n
        u(i,j,0)   = u(i,j,1)
        u(i,j,l+1) = u(i,j,l)
        v(i,j,0)   = v(i,j,1)
        v(i,j,l+1) = v(i,j,l)
        w(i,j,l)   = 0.0
        w(i,j,0)   = 0.0
        p(i,j,l+1) = 0.0
        p(i,j,0)   = 0.0
        t(i,j,l+1) = t(i,j,l)
        t(i,j,0)   = t(i,j,1)
        s(i,j,l+1) = s(i,j,l)
        s(i,j,0)   = s(i,j,1)
     enddo
  enddo

  do i=1,n
     do j=1,m
        do k=1,l
           if (landm(i,j,k).eq.1) then
              u(i,j,k) = 0.0
              v(i,j,k) = 0.0
              u(i-1,j,k) = 0.0
              v(i-1,j,k) = 0.0
              u(i,j-1,k) = 0.0
              v(i,j-1,k) = 0.0
              u(i-1,j-1,k) = 0.0
              v(i-1,j-1,k) = 0.0
           endif
        enddo
     enddo
  enddo

end SUBROUTINE usol

!****************************************************************************
SUBROUTINE solu(un,u,v,w,p,t,s)
  !     Go from u,v,h,t to un
  use m_usr
  implicit none
  !     IMPORT/EXPORT
  real    un(ndim)
  real    u(0:n  ,0:m,0:l+1), v(0:n,0:m  ,0:l+1)
  real    w(0:n+1,0:m+1,0:l  ), p(0:n+1,0:m+1,0:l+1)
  real    t(0:n+1,0:m+1,0:l+1), s(0:n+1,0:m+1,0:l+1)
  integer find_row2
  !     LOCAL
  integer i,j,k

  do k = 1, l
     do j = 1, m
        do i = 1, n
           un(find_row2(i,j,k,UU)) = u(i,j,k)
           un(find_row2(i,j,k,VV)) = v(i,j,k)
           un(find_row2(i,j,k,WW)) = w(i,j,k)
           un(find_row2(i,j,k,PP)) = p(i,j,k)
           un(find_row2(i,j,k,TT)) = t(i,j,k)
           un(find_row2(i,j,k,SS)) = s(i,j,k)
        enddo
     enddo
  enddo

end SUBROUTINE solu

!****************************************************************************
SUBROUTINE stpnt!(un)
  use m_usr
  use m_atm
  implicit none

  !*******************************************************
  !     PARAMETERS:
  !*******************************************************
  ! when data are used, tmax comes from windfit
  ! otherwise tmax comes from wfun...

  par(AL_T)   =  0.1/(2*omegadim*rhodim*hdim*udim*dz*dfzT(l))
  par(RAYL)   =  alphaT*gdim*hdim/(2*omegadim*udim*r0dim)   ! Ra ~ 0.422
  par(EK_V)   =  av/(2*omegadim*hdim*hdim)                  ! E_V
  par(EK_H)   =  ah/(2*omegadim*r0dim*r0dim)                ! E_H

  ! !MdT \/
  ! par(RAYL)   = par(RAYL)*par(EK_H)            ! Rescaled Ra (2DMOC case)
  ! !MdT /\

  par(ROSB)   =  udim/(2*omegadim*r0dim)                    ! Rossby Number
  par(HMTP)   =  0.0
  par(SUNP)   =  0.0
  par(PE_H)   =  kappah/(udim*r0dim)            ! P_H0
  par(PE_V)   =  kappav*r0dim/(udim*hdim*hdim)  ! P_V0
  par(P_VC)   =  2.5e+04*par(PE_V)           ! 5.0  ! P_VC
  par(LAMB)   =  alphaS/alphaT               ! lambda
  par(SALT)   =  0.0                         ! gamma
  par(WIND)   =  0.0                         ! wind h
  par(TEMP)   =  0.0                         ! eta_T
  par(BIOT)   =  r0dim/(75.*3600.*24.*udim)  ! nonlinearity in T,S equations
  par(COMB)   =  0.0                         ! combined continuation
  par(NLES)   =  0.0                         ! nonlinear equation of state
  par(CMPR)   =  0.0
  par(ALPC)   =  1.0
  par(ENER)   =  1.0e+02
  par(MIXP)   =  0.0       ! ATvS-Mix
  par(MKAP)   =  0.0       ! Gent-Mcwilliams  ! ATvS-Mix
  par(SPL1)   =  2.0e+03   ! 1.25  !tanh      ! ATvS-Mix
  par(SPL2)   =  0.01      ! neutral physics

  ! Set mixing parameters
  call vmix_par

end SUBROUTINE stpnt

!******************************************************************
SUBROUTINE atmos_coef
  use m_usr
  use m_atm
  implicit none
  !     LOCAL
  real muoa,dzne
  integer j
  dzne = dz*dfzT(l)
  muoa = rhoa*ch*cpa*uw
  amua = (arad+brad*t0)/muoa
  bmua = brad/muoa
  Aa   = uatm*rhoa*cpa*hdima/(r0dim*muoa)
  Ai   = rhoa*hdima*cpa*udim/(r0dim*muoa)
  Ad   = rhoa*hdima*cpa*d0/(muoa*r0dim*r0dim)
  As   = sun0*(1 - c0)/(4*muoa)

  Os   = sun0*c0/4*QTnd
  Ooa  = muoa*QTnd

  nus  = 0.0 ! postpone until set_atmos_parameters is called

  ! lvsc = rhodim*lv*r0dim/(udim*cp0*rhodim*hdim*dzne)
  lvsc = 0.0 ! ...
  DO j = 1,m
     !       albe(j) = 0.15 + 0.05 * cos (y(j))
     dat(j)  = 0.9 + 1.5 * exp(-12*y(j)*y(j)/pi)
     ! suno(j) = Os*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
     ! suna(j) = As*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
     suno(j) = Os*(1-.482*(3*sin(y(j))**2-1.)/2.)
     suna(j) = As*(1-.482*(3*sin(y(j))**2-1.)/2.)
  ENDDO

  DO j = 0,m
     davt(j) = 0.9 + 1.5 * exp(-12*yv(j)*yv(j)/pi)
  ENDDO

  ! write(*,*) 'Ocean-Atmosphere pars:     dzne=', dzne,' hdim=', hdim
  ! write(*,*) '                            nus=', nus, ' lvsc=', lvsc
  ! write(*,*) '                            Ooa=', Ooa, ' muoa=', muoa

END SUBROUTINE atmos_coef

!******************************************************************
SUBROUTINE writeparameters
  use m_par
  implicit none
  integer i
  open(55,position='append')
  write(55,*) "new parameters"
  do i=1,30
     write(55,*) i,par(i)
  enddo
  close(55)
END SUBROUTINE writeparameters

!******************************************************************
subroutine write_levitus(filename)
! debugging subroutine for parallelization

  use m_usr

  implicit none

  character(len=19) :: filename !something of the form 'debug_levitusXX.txt'
  integer :: i,j

  write(*,*) "writing levitus fields to file '",rundir//filename,"'"
  open(unit=42,file=rundir//filename,status='replace')

  write(42,*) "domain bounds: "
  write(42,*) "xmin = ",xmin
  write(42,*) "xmax = ",xmax
  write(42,*) "ymin = ",ymin
  write(42,*) "ymax = ",ymax

  write(42,*) "i j landm x y taux tauy tatm emip"
  do j=1,m
     do i=1,n
        write(42,999) i,j,landm(i,j,l),x(i),y(j),taux(i,j),tauy(i,j),tatm(i,j),emip(i,j)
     end do
  end do

999 format(3I4,6E16.6)

  close(42)

end subroutine write_levitus
