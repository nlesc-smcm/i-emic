#include "fdefs.h"

!**************************************************************************
!! initialize THCM: input: number of grid-points in x,y and z-direction,
!! bounds of the domain (formerly set in usr.com)
SUBROUTINE init(a_n,a_m,a_l,a_nmlglob,&
     a_xmin,a_xmax,a_ymin,a_ymax,&
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
  integer(c_int) :: a_periodic
  integer(c_int), dimension((a_n+2)*(a_m+2)*(a_l+2)) :: a_landm
  real(c_double), dimension(a_n*a_m) :: a_taux,a_tauy
  real(c_double), dimension(a_n*a_m) :: a_tatm,a_emip,a_spert

  integer :: i,j,k,pos

  _INFO_('THCM: init... ')

  nmlglob = a_nmlglob

  xmin    = a_xmin
  xmax    = a_xmax
  ymin    = a_ymin
  ymax    = a_ymax

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
  call allocate_atm(n,m,l)
  call allocate_res(ndim)
  call allocate_mix()

  ! Fill the local landmask array. Replace all PERIO
  ! by OCEAN if periodic is false, so that THCM does
  ! not know about global periodicity (which would
  ! mess up the matrix partitioning).
  pos = 1
  do k = 0, l+la+1
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
  landm(:,:,l+la+1) = LAND

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

  call grid         !
  call stpnt        !
  call mixe         !
  call vmix_init    ! ATvS-Mix  USES LANDMASK
  call atmos_coef   !
  call forcing      ! USES LANDMASK

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

!****************************************************************************
SUBROUTINE mixe
  use m_usr
  implicit none
  integer i,j,k
  do i=1,n
     do j=1,m
        emix(i,j,0) = 0.0
        do k=1,l-1
           emix(i,j,k)= - sin(pi*(x(i)-xmin)/(xmax-xmin))*   &
                (cos((pi/2)*(y(j)-ymin)/(ymax-ymin)) + 0.2)         &
                *zw(k)
        enddo
        emix(i,j,l) = 0.0
     enddo
  enddo
  !
end SUBROUTINE mixe

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
SUBROUTINE getdeps(o_Ooa, o_Os, o_gamma, o_eta, o_lvscq)
  !     interface to get Ooa and other dependencies on external model
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_atm
  implicit none
  real(c_double) o_Ooa, o_Os, o_gamma, o_eta, o_lvscq
  o_Ooa   = Ooa
  o_Os    = Os
  o_gamma = par(COMB) * par(SALT) * nus * qdim
  o_eta   = eta
  o_lvscq = lvsc * qdim
end subroutine getdeps

!**********************************************************
SUBROUTINE get_constants(o_r0dim, o_udim, o_hdim)
  !     interface to get a few model constants
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_atm
  implicit none
  real(c_double) o_r0dim, o_udim, o_hdim

  o_r0dim = r0dim;
  o_udim  = udim;
  o_hdim  = hdim;

end subroutine get_constants

!**********************************************************
SUBROUTINE set_ep_constants(i_qdim, i_nuq, i_eta, i_dqso)
  !     interface to get a few model constants
  use, intrinsic :: iso_c_binding
  use m_usr
  use m_atm
  implicit none
  real(c_double) i_qdim, i_nuq, i_eta, i_dqso

  qdim = i_qdim 
  nuq  = i_nuq  
  eta  = i_eta  
  dqso = i_dqso

end subroutine set_EP_constants

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
SUBROUTINE matrix(un, sig1, sig2)
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
  real(c_double) :: sig1,sig2
  real time0, time1
  integer i,j,k,k1,row,ii,jj,kk,find_row2
  integer ix, iy, iz,ie

  !     clean old arrays:
  _DEBUG_("Zero out matrix arrays...")
  coB  = 0
  Al   = 0
  begA(1:ndim+1) = 0
  coA(1:maxnnz)  = 0.D0
  jcoA(1:maxnnz) = 0

  _DEBUG_("Build diagonal matrix B...")
  call fillcolB
  _DEBUG_("Build linear part of Jacobian...")
  call lin
#ifndef THCM_LINEAR
  _DEBUG_("Build nonlinear part of Jacobian...")
  call nlin_jac(un)
#endif
  !{ removing tons of things for eigen-analysis test...
#if 0
  !Euv
  Al(:,:,:,:,UU:VV,WW) = 0.0
  !BTS
  Al(:,:,:,:,WW,TT:SS) = 0.0
  !Buv, Bw
  !Al(:,:,:,:,TT:SS,UU::WW) = 0.0
  !Gw
  !Al(:,:,:,:,WW,PP) = 0.0
  !Dw
  !Al(:,:,:,:,PP,WW) = 0.0
#endif
  !}


  ! ATvS-Mix ---------------------------------------------------------------------
  if (vmix_flag.ge.1) then
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
  endif
  ! --------------------------------------------------------------------- ATvS-Mix

  do k = 1, l+la
     do j = 1, m
        do i = 1, n
           do ii = 1,nun
              row = find_row2(i,j,k,ii)
              Al(i,j,k,5,ii,ii) = Al(i,j,k,5,ii,ii) - sig1*coB(row)
           end do

           ! note that B is 0 in W/P points, but we add something there, too
           ! to make sure the diagonal gets included in the matrix (sig2~mach.eps)
           row = find_row2(i,j,k,WW)
           Al(i,j,k,5,WW,WW) = Al(i,j,k,5,WW,WW) - sig2
           row = find_row2(i,j,k,PP)
           Al(i,j,k,5,PP,PP) = Al(i,j,k,5,PP,PP) - sig2
        enddo
     enddo
  enddo
  
  call boundaries

  call assemble
  !write(*,*) "In fortran's matrix() maxval TT =", maxval(Al(:,:,:,:,TT,:))
  !write(*,*) "In fortran's matrix() maxval SS =", maxval(Al(:,:,:,:,SS,:))

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
  real mix(ndim) ! ATvS-Mix
  real    Au(ndim), time0, time1
  integer i,j,k,k1,row,find_row2, mode, iter

  !call writeparameters
  mix  = 0.0
  Al   = 0
  begA = 0
  coA  = 0
  jcoA = 0
  call lin              !
  write(*,*) 'T(n,m,l)', un(find_row2(n,m,l,TT))
#ifndef THCM_LINEAR
  call nlin_rhs(un)
#endif
  call forcing          !
  call boundaries       !
  call assemble
  call matAvec(un,Au)   !
  ! ATvS-Mix ---------------------------------------------------------------------
  if (vmix_flag.ge.1) then
     mode=vmix_fix
     call cpu_time(time0)
     if (vmix_out.gt.0) write(99,'(a26)')'MIX| rhs...               '
     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     fix:   ', vmix_fix

     if ((vmix_fix.eq.0).and.(vmix_flag.ge.2)) call vmix_control(un)

     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     temp:  ', vmix_temp
     if (vmix_out.gt.0) write(99,'(a16,i10)') 'MIX|     salt:  ', vmix_salt

     if (((vmix_temp.eq.1).or.(vmix_salt.eq.1)).and.(vmix_dim.gt.0)) then
        call vmix_fun(un,mix,vmix_fix)
     endif
     call cpu_time(time1)
     vmix_time=vmix_time+time1-time0
     if (vmix_out.gt.0) write (99,'(a26,f10.3)') 'MIX|    ...rhs done',time1-time0
  endif
  ! --------------------------------------------------------------------- ATvS-Mix

  _DEBUG2_("p0 = ", p0) ! Residue Continuation

  B = -Au - mix + Frc - p0*(1- par(RESC))*ures

#if 1
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
  ! |      Al(i,j,k,loc,A,B) = c                                          |
  ! |       where loc is one of the locations below:                      |
  ! |     +----------++-------++----------+                               |
  ! |     | 12 15 18 || 3 6 9 || 21 24 27 |                               |
  ! |     | 11 14 17 || 2 5 8 || 20 23 26 |                               |
  ! |     | 10 13 16 || 1 4 7 || 19 22 25 |                               |
  ! |     |  below   || center||  above   |                               |
  ! |     +----------++-------++----------+                               |
  ! |                                                                     |
  ! |     For instance, Al(i,j,k,14,A,B) = c is                           |
  ! |     d/dt A|(i,j,k) = c*B|(i,j,k-1) + ...                            |
  ! +---------------------------------------------------------------------+
  use m_usr
  use m_atm
  implicit none
  !     LOCAL
  real,target :: u(n,m,l,np),uy(n,m,l,np),ucsi(n,m,l,np),&
       &      uxx(n,m,l,np),uyy(n,m,l,np),uzz(n,m,l,np),&
       &      uxs(n,m,l,np),fu(n,m,l,np),px(n,m,l,np)
  real ub(n,m,l,np),vb(n,m,l,np),sc(n,m,l,np),tcb(n,m,l,np)
  real    yc(n,m,la,np),yc2(n,m,la,np),yxx(n,m,la,np),yyy(n,m,la,np)
  real    EH,EV,ph,pv,Ra,lambda, bi, ahcor, dedt
  real    hv(n,m,l,np),yadv(n,m,la,np)
  real    uxxc(n,m,l,np),uyyc(n,m,l,np),vxxc(n,m,l,np),vyyc(n,m,l,np)
  real    xes,rintb,rwint !, rintt ! ATvS-Mix

  ! original version:
  !      equivalence (u, v), (uy, vy), (ucsi, vcsi), (uxx, vxx, txx, uxc)
  !      equivalence (uyy, vyy, tyy, vyc), (uzz, vzz, tzz, wzc)
  !      equivalence (uxs, vxs, tbc, tyc)
  !      equivalence (fu, fv, tc), (px, py, pz)

  ! we use pointers instead:
  real,dimension(:,:,:,:),pointer :: v,vy,vcsi,vxx,txx,uxc,&
       &     vyy,tyy,vyc,vzz,tzz,wzc,vxs,tbc,tyc,fv,tc,py,pz

  v=>u
  vy=>uy
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
  call uderiv(1,ub)     
  call uderiv(2,uxx)
  call uderiv(3,uyy)
  call uderiv(4,uzz)
  call uderiv(5,ucsi)
  call uderiv(6,vxs)
  call uderiv(7,u)      
  call coriolis(1,fv)
  call gradp(1,px)
  Al(:,:,1:l,:,UU,UU) = -EH * (uxx+uyy+ucsi) -EV * uzz ! + rintt*u ! ATvS-Mix
  Al(:,:,1:l,:,UU,VV) = -fv - EH*vxs
  ! Al(:,:,1:l,:,UU,VV) = - EH*vxs ! for 2DMOC case
  Al(:,:,1:l,:,UU,PP) =  px

  ! ------------------------------------------------------------------
  ! v-equation
  ! ------------------------------------------------------------------
  call vderiv(1,vb )    
  call vderiv(2,vxx)
  call vderiv(3,vyy)
  call vderiv(4,vzz)
  call vderiv(5,vcsi)
  call vderiv(6,uxs)
  call vderiv(7,v)      
  call coriolis(2,fu)
  call gradp(2,py)
  Al(:,:,1:l,:,VV,UU) =  fu - EH*uxs
  ! Al(:,:,1:l,:,VV,UU) =  - EH*uxs ! for 2dMOC case
  Al(:,:,1:l,:,VV,VV) = -EH*(vxx + vyy + vcsi) - EV*vzz !+ rintt*v ! ATvS-Mix
  Al(:,:,1:l,:,VV,PP) =  py

  ! ------------------------------------------------------------------
  ! w-equation
  ! ------------------------------------------------------------------
  call gradp(3,pz)
  call tderiv(6,tbc)
  Al(:,:,1:l,:,WW,PP) =  pz
  Al(:,:,1:l,:,WW,TT) = -Ra *(1. + xes*alpt1) * tbc/2.
  Al(:,:,1:l,:,WW,SS) =  lambda * Ra * tbc/2.

  ! ------------------------------------------------------------------
  ! p-equation
  ! ------------------------------------------------------------------
  call pderiv(1,uxc)
  call pderiv(2,vyc)
  call pderiv(3,wzc)
  Al(:,:,1:l,:,PP,UU) = uxc
  Al(:,:,1:l,:,PP,VV) = vyc
  Al(:,:,1:l,:,PP,WW) = wzc

  ! ------------------------------------------------------------------
  ! T-equation
  ! ------------------------------------------------------------------
  call tderiv(1,tc )
  call tderiv(2,sc )
  call tderiv(3,txx)
  call tderiv(4,tyy)
  call tderiv(5,tzz)
  call tderiv(7,tcb)

  ! dependence of TT on TT through latent heat due to evaporation
  dedt = par(COMB) * par(TEMP) * lvsc * qdim * & 
       eta * (deltat / qdim) * dqso

  write(*,*) 'dedt=', dedt, ' eta=', eta, ' dqso=',dqso

  if (la > 0) then ! deprecated local atmosphere
     Al(:,:,1:l,:,TT,TT) = - ph * (txx + tyy) - pv * tzz + Ooa*tc
  else if (coupled_atm.eq.1) then ! coupled with external atmos
     Al(:,:,1:l,:,TT,TT) = - ph * (txx + tyy) - pv * tzz + &
          Ooa*tc! + dedt*sc
  else
     Al(:,:,1:l,:,TT,TT) = - ph * (txx + tyy) - pv * tzz + TRES*bi*tc
  endif

  ! ------------------------------------------------------------------
  ! S-equation
  ! ------------------------------------------------------------------
  ! dependence of SS on TT through evaporation
  dedt = par(COMB) * par(SALT) * nus * qdim * &
       eta * (deltat / qdim) * dqso
  
  if (coupled_atm.eq.1) then
     Al(:,:,1:l,:,SS,SS) = - ph * (txx + tyy) - pv * tzz 

     ! minus sign and nondim added (we take -Au in rhs computation)
     Al(:,:,1:l,:,SS,TT) = - dedt * sc !  * (r0dim / udim) 
  else
     Al(:,:,1:l,:,SS,SS) = - ph * (txx + tyy) - pv * tzz + SRES*bi*sc
  endif
  
  ! ------------------------------------------------------------------
  ! atmosphere layer
  ! ------------------------------------------------------------------
  if (la > 0) then
     call yderiv(1,yc )
     call yderiv(4,yc2 )
     call yderiv(2,yxx)
     call yderiv(3,yyy)
     call yderiv(5,yadv)
     Al(:,:,l+1:l+la,5,UU,UU)  =  1.
     Al(:,:,l+1:l+la,5,VV,VV)  =  1.
     Al(:,:,l+1:l+la,5,WW,WW)  =  1.
     Al(:,:,l+1:l+la,5,PP,PP)  =  1.
     Al(:,:,l+1:l+la,5,SS,SS)  =  1.
     Al(:,:,l+1:l+la,:,TT,TT)  = - Ad * (yxx + yyy) - yc + &
                                   Aa * yadv + bmua*yc2
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
  real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
  real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
  real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
  real    rho(0:n+1,0:m+1,0:l+la+1)
  real,target ::    utx(n,m,l,np),vty(n,m,l,np),wtz(n,m,l,np)
  real    t2r(n,m,l,np),t3r(n,m,l,np)
  real    cat1(n,m,l,np),cas1(n,m,l,np)
  real    cat2(n,m,l,np),cas2(n,m,l,np)
  real    uux(n,m,l,np),uvy1(n,m,l,np),uwz(n,m,l,np),uvy2(n,m,l,np)
  real    uvx(n,m,l,np),vvy(n,m,l,np),vwz(n,m,l,np),ut2(n,m,l,np)
  real    bolt(n,m,la,np)
  real    lambda,epsr,Ra,xes,pvc1,pvc2, pv

  real,dimension(:,:,:,:),pointer ::    usx,vsy,wsz

  usx=>utx
  vsy=>vty
  wsz=>wtz

  lambda = par(LAMB)
  epsr   = par(ROSB)
  Ra     = par(RAYL)
  xes    = par(NLES)
  pvc1   = par(P_VC)
  pv     = par(PE_V)
  pvc2   = pv*(1.0 - par(ALPC))*par(ENER)
  call usol(un,u,v,w,p,t,s)
  rho    = lambda*s - t *( 1 + xes*alpt1) - &
           xes*t*t*alpt2+xes*t*t*t*alpt3

  ! ------------------------------------------------------------------
  ! u-equation
  ! ------------------------------------------------------------------
#ifndef NO_UVNLIN
  call unlin(1,uux,u,v,w)
  call unlin(3,uvy1,u,v,w)
  call unlin(5,uwz,u,v,w)
  call unlin(7,uvy2,u,v,w)
  Al(:,:,1:l,:,UU,UU) = Al(:,:,1:l,:,UU,UU) + epsr * (uux + uvy1 + uwz + uvy2)
#endif

  ! ------------------------------------------------------------------
  ! v-equation
  ! ------------------------------------------------------------------
#ifndef NO_UVNLIN
  call vnlin(1,uvx,u,v,w)
  call vnlin(3,vvy,u,v,w)
  call vnlin(5,vwz,u,v,w)
  call vnlin(7,ut2,u,v,w)
  Al(:,:,1:l,:,VV,UU) = Al(:,:,1:l,:,VV,UU) + epsr *ut2
  Al(:,:,1:l,:,VV,VV) = Al(:,:,1:l,:,VV,VV) + epsr*(uvx + vvy + vwz)
#endif

  ! ------------------------------------------------------------------
  ! w-equation
  ! ------------------------------------------------------------------
  call wnlin(2,t2r,t)
  call wnlin(4,t3r,t)
  Al(:,:,1:l,:,WW,TT) = Al(:,:,1:l,:,WW,TT) - Ra*xes*alpt2*t2r &
                                            + Ra*xes*alpt3*t3r

  ! ------------------------------------------------------------------
  ! T-equation
  ! ------------------------------------------------------------------
#ifndef NO_TSNLIN
  call tnlin(3,utx,u,v,w,t,rho)
  call tnlin(5,vty,u,v,w,t,rho)
  call tnlin(7,wtz,u,v,w,t,rho)
  Al(:,:,1:l,:,TT,TT) = Al(:,:,1:l,:,TT,TT)+ utx+vty+wtz        ! ATvS-Mix
#endif

  ! ------------------------------------------------------------------
  ! S-equation
  ! ------------------------------------------------------------------
#ifndef NO_TSNLIN
  call tnlin(3,usx,u,v,w,s,rho)
  call tnlin(5,vsy,u,v,w,s,rho)
  call tnlin(7,wsz,u,v,w,s,rho)
  Al(:,:,1:l,:,SS,SS) = Al(:,:,1:l,:,SS,SS)+ usx+vsy+wsz        ! ATvS-Mix
#endif

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
  real    u(0:n  ,0:m,0:l+la+1), v(0:n, 0:m  ,0:l+la+1)
  real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
  real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
  real    rho(0:n+1,0:m+1,0:l+la+1)
  real,target :: urTx(n,m,l,np),Utrx(n,m,l,np),&
       &        vrTy(n,m,l,np),Vtry(n,m,l,np)
  real,target :: wrTz(n,m,l,np),Wtrz(n,m,l,np)
  real    cat1(n,m,l,np),cat2(n,m,l,np),&
       &        cat3(n,m,l,np),cat4(n,m,l,np)
  real    t2r(n,m,l,np),t3r(n,m,l,np)

  real    cas1(n,m,l,np),cas2(n,m,l,np),&
       &        cas3(n,m,l,np),cas4(n,m,l,np)
  real    bolt(n,m,la,np)
  real    lambda,epsr,Ra,xes,pvc1,pvc2,pv
  real uux(n,m,l,np),uvy1(n,m,l,np),uwz(n,m,l,np),uvy2(n,m,l,np)
  real uvx(n,m,l,np),vvry(n,m,l,np),vwz(n,m,l,np),wvrz(n,m,l,np)
  real Urux(n,m,l,np),Urvy1(n,m,l,np),Urwz(n,m,l,np),Urvy2(n,m,l,np)
  real uVrx(n,m,l,np),Vrvy(n,m,l,np),Vrwz(n,m,l,np),Urt2(n,m,l,np)

  real,dimension(:,:,:,:),pointer :: urSx,Usrx,vrSy,Vsry,wrSz,Wsrz
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
  pvc1   = par(P_VC)
  pv     = par(PE_V)
  pvc2   = pv*(1.0 - par(ALPC))*par(ENER)
  call usol(un,u,v,w,p,t,s)
  rho    = lambda*s - t *( 1 + xes*alpt1) - &
       &            xes*t*t*alpt2+xes*t*t*t*alpt3

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
  Al(:,:,1:l,:,UU,UU)  =  Al(:,:,1:l,:,UU,UU) + epsr * (Urux + uvy1 + uwz + uvy2)
  Al(:,:,1:l,:,UU,VV)  =  Al(:,:,1:l,:,UU,VV) + epsr * (Urvy1 + Urvy2)
  Al(:,:,1:l,:,UU,WW)  =  Al(:,:,1:l,:,UU,WW) + epsr *  Urwz
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
  Al(:,:,1:l,:,VV,UU) =   Al(:,:,1:l,:,VV,UU) + epsr * (Urt2 + uVrx)
  Al(:,:,1:l,:,VV,VV) =   Al(:,:,1:l,:,VV,VV) + epsr * (uvx + Vrvy + vwz)
  Al(:,:,1:l,:,VV,WW) =   Al(:,:,1:l,:,VV,WW) + epsr * Vrwz
#endif

  ! ------------------------------------------------------------------
  ! w-equation
  ! ------------------------------------------------------------------
  call wnlin(1,t2r,t)
  call wnlin(3,t3r,t)
  Al(:,:,1:l,:,WW,TT) = Al(:,:,1:l,:,WW,TT) - Ra*xes*alpt2*t2r &
                                            + Ra*xes*alpt3*t3r

  ! ------------------------------------------------------------------
  ! T-equation
  ! ------------------------------------------------------------------
#ifndef NO_TSNLIN
  call tnlin(2,urTx,u,v,w,t,rho)
  call tnlin(3,Utrx,u,v,w,t,rho)
  call tnlin(4,vrTy,u,v,w,t,rho)
  call tnlin(5,Vtry,u,v,w,t,rho)
  call tnlin(6,wrTz,u,v,w,t,rho)
  call tnlin(7,Wtrz,u,v,w,t,rho)
  Al(:,:,1:l,:,TT,UU) = Al(:,:,1:l,:,TT,UU) + urTx
  Al(:,:,1:l,:,TT,VV) = Al(:,:,1:l,:,TT,VV) + vrTy
  Al(:,:,1:l,:,TT,WW) = Al(:,:,1:l,:,TT,WW) + wrTz
  Al(:,:,1:l,:,TT,TT) = Al(:,:,1:l,:,TT,TT) + Utrx + Vtry + Wtrz        ! ATvS-Mix
#endif

  ! ------------------------------------------------------------------
  ! S-equation
  ! ------------------------------------------------------------------
#ifndef NO_TSNLIN
  call tnlin(2,urSx,u,v,w,s,rho)
  call tnlin(3,Usrx,u,v,w,s,rho)
  call tnlin(4,vrSy,u,v,w,s,rho)
  call tnlin(5,Vsry,u,v,w,s,rho)
  call tnlin(6,wrSz,u,v,w,s,rho)
  call tnlin(7,Wsrz,u,v,w,s,rho)
  Al(:,:,1:l,:,SS,UU) = Al(:,:,1:l,:,SS,UU) + urSx
  Al(:,:,1:l,:,SS,VV) = Al(:,:,1:l,:,SS,VV) + vrSy
  Al(:,:,1:l,:,SS,WW) = Al(:,:,1:l,:,SS,WW) + wrSz
  Al(:,:,1:l,:,SS,SS) = Al(:,:,1:l,:,SS,SS) + Usrx + Vsry + Wsrz
#endif

end SUBROUTINE nlin_jac
!****************************************************************************
SUBROUTINE usol(un,u,v,w,p,t,s)
  !     Go from un to u,v,t,h
  use m_usr
  implicit none
  !     IMPORT/EXPORT
  real    un(ndim)
  real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
  real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
  real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
  !     LOCAL
  integer i,j,k,row
  !     EXTERNAL
  integer find_row2

  u = 0.0
  v = 0.0
  w = 0.0
  p = 0.0
  t = 0.0
  s = 0.0
  do k = 1, l+la
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
  do k = 1, l+la
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
           v(N,j,k)= 0.0
           p(0,j,k)  = 0.0
           p(N+1,j,k)= 0.0
           t(0,j,k)  = t(1,j,k)
           t(N+1,j,k)= t(N,j,k)
           s(0,j,k)  = s(1,j,k)
           s(N+1,j,k)= s(N,j,k)
        endif
     enddo
  enddo
  do k = 1, l+la
     do i = 1, n
        u(i,0,k)  = 0.0
        u(i,M,k)= 0.0
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
        IF (la == 0) t(i,j,l+1) = t(i,j,l)
        t(i,j,0)   = t(i,j,1)
        s(i,j,l+1) = s(i,j,l)
        s(i,j,0)   = s(i,j,1)
     enddo
  enddo
18 format(i4,i4,4(g12.4))
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
  real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
  real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
  real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
  integer find_row2
  !     LOCAL
  integer i,j,k,row

  do k = 1, l+la
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
  integer i,j,k,row,find_row2
  !     real    un(ndim)

  !*******************************************************
  !     PARAMETERS:
  !*******************************************************
  ! when data are used, tmax comes from windfit
  ! otherwise tmax comes from wfun...

  par(AL_T)   =  0.1/(2*omegadim*rhodim*hdim*udim*dz*dfzT(l))
  par(RAYL)   =  alphaT*gdim*hdim/(2*omegadim*udim*r0dim)   ! Ra
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
  par(P_VC)   =  2.5e+04*par(PE_V)   ! 5.0  ! P_VC
  par(LAMB)   =  alphaS/alphaT     ! lambda
  par(SALT)   =  0.0               ! gamma
  par(WIND)   =  0.0               ! wind h
  par(TEMP)   =  0.0               ! eta_T
  par(BIOT)   =  r0dim/(75.*3600.*24.*udim) ! nonlinearity in T,S equations
  par(COMB)   =  0.0               ! combined continuation
  par(NLES)   =  0.0       ! nonlinear equation of state
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
  integer i,j
  dzne = dz*dfzT(l)
  muoa = rhoa*ch*cpa*uw
  amua = (arad+brad*t0)/muoa
  bmua = brad/muoa
  Aa   = uatm*rhoa*cpa*hdima/(r0dim*muoa)
  Ai   = rhoa*hdima*cpa*udim/(r0dim*muoa)
  Ad   = rhoa*hdima*cpa*d0/(muoa*r0dim*r0dim)
  As   = sun0*(1 - c0)/(4*muoa)
  Os   = sun0*c0*r0dim/(4*udim*hdim*dzne*rhodim*cp0)
  Ooa  = muoa*r0dim/(udim*cp0*rhodim*hdim*dzne)
  nus  = s0 * r0dim / (udim*deltas*hdim*dzne  ) ! without qdim!
  lvsc = lv * r0dim / (udim*cp0*hdim*dzne ) ! without qdim!
  DO j = 1,m
     !       albe(j) = 0.15 + 0.05 * cos (y(j))
     albe(j) = 0.3
     dat(j)  = 0.9 + 1.5 * exp(-12*y(j)*y(j)/pi)
     suno(j) = Os*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
     suna(j) = As*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
  ENDDO

  DO j = 0,m
     davt(j) = 0.9 + 1.5 * exp(-12*yv(j)*yv(j)/pi)
  ENDDO

  open(8, file = rundir//'suno.txt')
  write(8, *) suno
  close(8)
  
  write(*,*) 'Ocean-Atmosphere pars:     dzne=', dzne,' hdim=', hdim
  write(*,*) '                            nus=', nus, ' lvsc=', lvsc
  write(*,*) '                            Ooa=', Ooa, ' muoa=', muoa    

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
