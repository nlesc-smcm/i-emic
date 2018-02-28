#include "fdefs.h"
! /**********************************************************************
! * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
! * Permission to use, copy, modify, redistribute is granted           *
! * as long as this header remains intact.                             *
! * contact: jonas@math.rug.nl                                         *
! **********************************************************************/

! This module contains variables from usr.F90 which are adjusted to
! represent the global domain in a parallel computation. This informa-
! tion is used for I/O purposes, i.e. to read land-cells etc. m_global
! should be used sparsely and with care to avoid confusions with the
! subdomain (usr.F90).

module m_global

  ! zmin,zmax and qz can be taken from usr.F90 because
  ! they are not used in the C++ part (at least not yet!)
  use m_par

  ! these things are shared between global and subdomain.
  use m_usr, only :                          &
       zmin, zmax, la,                       &
       hdim, qz,                             &
       alphaT, alphaS,                       &
       ih, vmix_GLB, tap, rho_mixing,        &
       itopo, flat, rd_mask,                 &
       coupled_atm,                          &
       TRES, SRES, iza, ite, its, rd_spertm, &
       rowintcon, f99, t0, s0

  use m_atm, only : Ooa, suno

  implicit none

  integer, parameter :: nf   =  4

  integer :: n = 0
  integer :: m = 0
  integer :: l = 0  !! these may be different from the ones in usr.F90!

  integer icp
  real,dimension(:),allocatable ::   u, up
  real,dimension(:,:),allocatable :: w
  real,dimension(nf,2) :: sig
  real :: xl, xlp, det, tval

  integer :: ndim = 0

  ! here are some things we cannot use from m_usr as they may (and are
  ! likely to) differ between the computational and global domains
  real :: xmin, xmax, ymin, ymax
  real, dimension(:),allocatable :: x,y,z,xu,yv,zw,ze,zwe
  real :: dx, dy, dz
  integer, dimension(:,:,:),allocatable :: landm
  real, dimension(:,:), allocatable :: taux, tauy, tatm, emip, spert
  logical ::  periodic
  real, dimension(:,:,:), allocatable :: internal_temp,internal_salt

  character(len=999) :: maskfile, spertmaskfile

contains

  !! allocate memory for the global data structures. You are allowed
  !! to specify m,n,l to be 0 so that no memory is allocated (most
  !! subdomains will only use xmin,xmax,ymin and ymax)
  subroutine initialize(a_n,a_m,a_l,&
       a_xmin,a_xmax,a_ymin,a_ymax,a_hdim,a_qz,&
       a_alphaT,a_alphaS,&
       a_ih,a_vmix_GLB,a_tap,a_rho_mixing,&
       a_periodic,a_itopo,a_flat,a_rd_mask,&
       a_TRES,a_SRES,a_iza,a_ite,a_its,a_rd_spertm,&
       a_coupled_atm)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int) :: a_n,a_m,a_l
    real(c_double) :: a_xmin,a_xmax,a_ymin,a_ymax,a_hdim,a_qz
    real(c_double) :: a_alphaT, a_alphaS
    integer(c_int) :: a_ih,a_vmix_GLB,a_tap,a_rho_mixing
    integer(c_int) :: a_periodic,a_itopo,a_flat,a_rd_mask
    integer(c_int) :: a_TRES,a_SRES,a_iza,a_ite,a_its,a_rd_spertm
    integer(c_int) :: a_coupled_atm

    xmin  = a_xmin
    xmax  = a_xmax
    ymin  = a_ymin
    ymax  = a_ymax
    hdim  = a_hdim
    qz    = a_qz

    alphaT   = a_alphaT
    alphaS   = a_alphaS

    ih       = a_ih
    vmix_GLB = a_vmix_GLB
    tap      = a_tap
    if (a_rho_mixing.ne.0) then
       rho_mixing = .true.
    else
       rho_mixing = .false.
    end if

    if (a_periodic .ne. 0) then
       periodic   = .true.
    else
       periodic   = .false.
    end if
    itopo    = a_itopo
    if (a_flat .ne. 0) then
       flat      = .true.
    else
       flat      = .false.
    end if
    if (a_rd_mask .ne. 0) then
       rd_mask   = .true.
    else
       rd_mask   = .false.
    end if

    TRES     = a_TRES
    SRES     = a_SRES
    iza      = a_iza
    ite      = a_ite
    its      = a_its
    if (a_rd_spertm .ne. 0) then
       rd_spertm = .true.
    else
       rd_spertm = .false.
    end if

    !========= I-EMIC coupling  ================================
    coupled_atm = a_coupled_atm
    !===========================================================

    ! if (n/=a_n .or. m/=a_m .or. l/=a_l) then
    
    n = a_n
    m = a_m
    l = a_l
    ndim = n*m*l*nun

    call deallocate_global

    write(*,*) 'allocating fortran arrays'
    allocate(u(ndim), up(ndim), w(ndim,nf))       
    allocate(landm(0:n+1,0:m+1,0:l+la+1))
    allocate(taux(n,m),tauy(n,m))
    allocate(tatm(n,m),emip(n,m),spert(n,m))
    allocate(internal_temp(n,m,l),internal_salt(n,m,l))
    internal_salt=0.0
    internal_temp=0.0
    landm = 0

    allocate(x(n),y(m),z(l),xu(0:n),yv(0:m),zw(0:l),ze(l),zwe(l))
    
    call g_grid
    ! end if ! new dimensions?

  end subroutine initialize

  subroutine deallocate_global

    implicit none

    write(*,*) 'deallocating fortran arrays'

    if (allocated(u)) deallocate(u)
    if (allocated(up)) deallocate(up)
    if (allocated(w)) deallocate(w)

    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(z)) deallocate(z)
    
    if (allocated(xu)) deallocate(xu)
    if (allocated(yv)) deallocate(yv)
    if (allocated(zw)) deallocate(zw)
    if (allocated(ze)) deallocate(ze)
    if (allocated(zwe)) deallocate(zwe)
    
    if (allocated(taux)) deallocate(taux)
    if (allocated(tauy)) deallocate(tauy)
    
    if (allocated(landm)) deallocate(landm)

    if (allocated(tatm))  deallocate(tatm)
    if (allocated(emip))  deallocate(emip)
    if (allocated(spert)) deallocate(spert)
    
    if (allocated(internal_temp)) deallocate(internal_temp)
    if (allocated(internal_salt)) deallocate(internal_salt)

  end subroutine deallocate_global

  !! this is meant as an interface for the C++ code, it just calls
  !! deallocate_global, really:
  subroutine finalize()

    implicit none

    call deallocate_global

  end subroutine finalize

  !! this is a copy of the corresponding function in matetc.f
  !! using the redefined dimensinos m,n,l from this module.
  !! It has to be kept consistent with the ordering used there !!!
  integer function find_row2(i,j,k,XX)
    ! find row in matrix A of variable XX at grid point (i,j,k)

    use m_par

    implicit none

    integer i,j,k,XX
    find_row2 = nun*((k-1)*n*m+ n*(j-1) + i-1)+XX

  end function find_row2

  !! this is a part of the subroutine grid from grid.f. We need it to
  !! initialize the vectors x,y,z...
  !! It should be kept consistent with that subroutine
  SUBROUTINE g_grid
    use m_par
    implicit none
    integer i,j,k
    !     EXTERNAL
    real  dfdz,fz

    write(*,*) '============GRID==========='
    write(*,10) xmin*180/pi,xmax*180/pi,ymin*180/pi,ymax*180/pi
10  format(1x,'configuration:[',f6.1,1x,f6.1,'] x [',f6.1,1x,f6.1,']')
    write(*,20) n,m,l
20  format(1x,'resolution:',i8,'x',i8,'x',i8)
    write(*,30) qz
30  format(1x,'stretching parameter',g12.4)
    write(*,*) '============GRID==========='

    dx = (xmax-xmin)/N
    dy = (ymax-ymin)/M
    dz = (zmax-zmin)/L
    DO i=1,n
       x(i) = (real(i)-0.5)*dx + xmin
       xu(i)= (real(i)    )*dx + xmin
    ENDDO

    xu(0) = xmin
    DO j=1,m
       y(j) = (real(j)-0.5)*dy + ymin
       yv(j)= (real(j)    )*dy + ymin
    ENDDO
    yv(0) = ymin
    DO k=1,l
       ze(k)  = (real(k)-0.5)*dz + zmin
       zwe(k) = (real(k)    )*dz + zmin
       z(k)   = fz(ze(k) ,qz)
       zw(k)  = fz(zwe(k),qz)
       ! compute derivatives of mapping at T- points
       !        dfzT(k) = dfdz(ze(k),qz)
       ! compute derivatives of mapping at w- points
       !        dfzW(k) = dfdz(zwe(k),qz)
    ENDDO
    zw(0) = zmin
    !      dfzw(0) = dfdz(zmin,qz)

  end subroutine g_grid

  !! this is supposed to be called once by the root
  !! proc with a standard 0:n+1,0:m+1,0:l+la+1 1D C
  !! array which is then distributed to all the subdomains
  !! and re-inserted into the subdomains when calling usrc:init
  subroutine get_landm(cland)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), dimension((n+2)*(m+2)*(l+la+2)) :: cland

    integer :: i,j,k,pos


    _INFO_('THCM: global.F90 get_landm...')
    call topofit

    pos = 1
    do k=0,l+la+1
       do j=0,m+1
          do i=0,n+1
             cland(pos) = landm(i,j,k)
             pos=pos+1
          end do
       end do
    end do
    _INFO_('THCM: global.F90 get_landm... done')
  end subroutine get_landm

  ! Obtain the landmask as it is right now in THCM,
  ! without reading from a file
  subroutine get_current_landm(cland)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), dimension((n+2)*(m+2)*(l+la+2)) :: cland

    integer :: i,j,k,pos

    _INFO_('THCM: global.F90 get_current_landm...')
    pos = 1
    do k=0,l+la+1
       do j=0,m+1
          do i=0,n+1
             cland(pos) = landm(i,j,k)
             pos=pos+1
          end do
       end do
    end do
    _INFO_('THCM: global.F90 get_current_landm... done')
  end subroutine get_current_landm

  !! Similar but vice versa, set landm from c array
  subroutine set_landm(cland)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), dimension((n+2)*(m+2)*(l+la+2)) :: cland
    integer :: i,j,k,pos


!    _INFO_('THCM: global.F90 set_landm...')

    landm = 0;

    pos = 1
    do k=0,l+la+1
       do j=0,m+1
          do i=0,n+1
             landm(i,j,k) = cland(pos)
             pos=pos+1
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

!    write(*,*) '_______ set_landm ____________', z(l)*hdim
!
!    do j = m+1, 0, -1
!       write(*,'(  92i1)') landm(:,j,l)
!    enddo


!    _INFO_('THCM: global.F90 set_landm... done')
  end subroutine set_landm

  !! this is supposed to be called once by the root
  !! proc with standard 1:n+1,0:m 1D C
  !!arrays which are then distributed to all the subdomains
  !! and re-inserted into the subdomains when calling usrc:init
  !  subroutine get_levitus(sres,ctatm,cemip)
  !
  !  implicit none
  !
  !  integer :: SRES
  !  real, dimension(n*m) :: ctatm, cemip
  !
  !  integer :: i,j,pos
  !
  !  if (SRES.eq.1) then
  !    call levitus_sal
  !  else
  !    ! read salt flux from fort.15
  !    call read_forcing(emip,15)
  !  end if
  !  call levitus_sst
  !
  !  pos = 1
  !  do j=1,m
  !    do i=1,n
  !      ctatm(pos) = tatm(i,j)
  !      cemip(pos) = emip(i,j)
  !      pos = pos+1
  !    end do
  !  end do
  !
  !  end subroutine get_levitus

  subroutine get_windfield(ctaux,ctauy)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), dimension(n*m) :: ctaux, ctauy
    real, dimension(n) :: tauz
    integer :: i,j,pos

    if(iza.ne.2) then
       call windfit        ! read data with subroutine from forcing.F90
       if(iza.eq.1) then   ! average zonally
          do j=1,m
             tauz(j) = 0.0
             do i=1,n
                tauz(j) = tauz(j) + taux(i,j)
             end do
             taux(1:n,j) = tauz(j)/n
             tauy(1:n,j) = 0.0
          end do
       end if
    else                  ! idealized wind, set in forcing.F90
       taux(1:n,1:m) = 0.0
       tauy(1:n,1:m) = 0.0
    end if

    pos = 1
    do j=1,m
       do i=1,n
          ctaux(pos) = taux(i,j)
          ctauy(pos) = tauy(i,j)
          pos = pos+1
       end do
    end do

  end subroutine get_windfield

  subroutine get_temforcing(ctatm)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), dimension(n*m) :: ctatm
    integer :: i,j,pos

    if (coupled_atm.ne.1) then
       if (ite.eq.0) then
          if (TRES.eq.0) then ! read heat flux and store it in tatm
             _INFO_("Expecting heat flux to be provided by Ocean hdf5 interface...")
          else                ! read sst
             _INFO_("Read sst from levitus...")
             call levitus_sst
          end if
       else if (ite.eq.1) then  ! idealized forcing, set in forcing.F90
          _INFO_("Using idealized temperature forcing, see forcing.F90")
          _INFO_(" For now we put tatm = 0")
          tatm(1:n,1:m) = 0.0
       else
          _INFO2_("Invalid option assigned to ite: ", ite)
       end if
    else ! accepting external atmosphere within I-EMIC
       _INFO_("Accepting external atmosphere, we are now part of the I-EMIC")
       _INFO_(" For now we put tatm = 0")
       tatm(1:n,1:m) = 0.0
    end if

    pos = 1
    do j=1,m
       do i=1,n
          ctatm(pos) = tatm(i,j)
          pos = pos+1
       end do
    end do

  end subroutine get_temforcing

  subroutine get_internal_temforcing(ctemp)

    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), dimension(n*m*l) :: ctemp
    integer :: i,j,k,pos

    !choose the followings

    write(f99,*) 'reading internal temperature forcing from "'//topdir//'levitus/new/t00an1'//'"'
    call levitus_internal(topdir//'levitus/new/t00an1',internal_temp,.false.,'TEMP')
    !call levitus_internal(topdir//'levitus/new/avtemp',internal_temp,.false.,'TEMP')

    pos = 1
    do k=1,l
       do j=1,m
          do i=1,n
             ctemp(pos) =internal_temp(i,j,k)
             pos = pos+1
          end do
       end do
    end do

  end  subroutine   get_internal_temforcing


  subroutine get_salforcing(cemip)

    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), dimension(n*m) :: cemip
    integer :: i,j,pos

    if (its.ne.1 .and. coupled_atm.eq.0) then
       ! read (virtual) salt flux and store it in emip
       if (SRES.eq.0) then 
          _INFO_("Expecting salinity flux to be provided by Ocean hdf5 interface...")
       else                ! read sss
          call levitus_sal
       end if
    else   ! idealized or coupled forcing, set in forcing.F90
       emip(1:n,1:m) = 0.0
    end if

    pos = 1
    do j=1,m
       do i=1,n
          cemip(pos) = emip(i,j)
          pos = pos+1
       end do
    end do

  end subroutine get_salforcing

  subroutine get_internal_salforcing(csalt)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), dimension(n*m*l) :: csalt
    integer :: i,j,k,pos

    !choose the followings
    write(f99,*) 'reading internal salt forcing from "'//topdir//'levitus/new/s00an1'//'"'
    call levitus_internal(topdir//'levitus/new/s00an1',internal_salt,.false.,'SALT')
    !call levitus_internal(topdir//'levitus/new/avsalt',internal_salt,.false.,'SALT')
    !chose the followings

    pos = 1
    do k=1,l
       do j=1,m
          do i=1,n
             csalt(pos) = internal_salt(i,j,k)
             pos = pos+1
          end do
       end do
    end do
  end subroutine get_internal_salforcing

  subroutine get_spert(cspert)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), dimension(n*m) ::  cspert
    integer :: i,j,pos
    if (rd_spertm) then
       call read_spertm
    else
       spert(1:n,1:m) = real(SRES)
    end if

    pos = 1
    do j=1,m
       do i=1,n
          cspert(pos) = spert(i,j)
          pos=pos+1
       end do
    end do

  end subroutine get_spert

  !! get levitus T and S for a certain month. These
  !! will be distributed by the C++ code and passed
  !! back into m_monthly by calling set_levitus
  !! 0<=month<=11 indicates which month you want.
  subroutine get_monthly_forcing(ctatm,cemip,ctaux,ctauy,month)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), dimension(n*m) :: ctatm, cemip, ctaux, ctauy
    integer(c_int), intent(in) :: month
    character(len=2) :: ibuf
    character(len=1024) :: fname

    real, dimension(n,m) :: taux_bak, tauy_bak, tatm_bak, emip_bak

    real :: temfun, salfun ! from forcing.F90

    integer :: i,j,pos

    ! if idealized forcing is desired, construct annual mean:
    if (ite==1) then
       do j=1,m
          do i=1,n
             tatm(i,j) = temfun(x(i),y(j))
          end do
       end do
    end if

    if (its==1) then
       do j=1,m
          do i=1,n
             emip(i,j) = salfun(x(i),y(j))
          end do
       end do
    end if

    ! use these arrays from now on:
    its=0
    ite=0

    ! keep annual mean data in backup locations
    taux_bak(:,:) = taux
    tauy_bak(:,:) = tauy
    tatm_bak(:,:) = tatm
    emip_bak(:,:) = emip

    write(ibuf,'(1I2.2)') month

    if (ite.eq.0) then
       fname = topdir//'levitus/monthly/t'//ibuf//'an1'
       write(f99,*) 'read levitus SST from file"'//trim(fname)//'"'
       call levitus_interpol(trim(fname),tatm,-5.,50.,l,1)
    else if (ite.eq.1) then
       ! idealized - just add a constant sin(t) term (this is only for testing...)
       tatm=sin(month/12*pi)
       ite=0 ! use tatm from now on
    end if

    if (its.eq.0 .and. coupled_atm.eq.0) then
       fname = topdir//'levitus/monthly/s'//ibuf//'an1'
       write(f99,*) 'read levitus S from file"'//trim(fname)//'"'
       call levitus_interpol(trim(fname),emip,30.,40.,l,1)
    else
       ! idealized - just add a constant sin(t) term (this is only for testing...)
       emip=sin(month/12*pi)
       its=0 ! use the emip array from now on!
    end if

    if (iza==0) then
       !wind (cf. forcing.F90)
       call windfit_monthly(month)
    else
       ! idealized - just add a constant sin(t) term (this is only for testing...)
       taux=sin(month/12*pi)
       tauy=0;
    end if

    pos = 1
    do j=1,m
       do i=1,n
          ctatm(pos) = tatm(i,j)
          cemip(pos) = emip(i,j)
          ctaux(pos) = taux(i,j)
          ctauy(pos) = tauy(i,j)
          pos = pos+1
       end do
    end do

    ! restore annual mean data in original arrays
    taux(:,:) = taux_bak
    tauy(:,:) = tauy_bak
    tatm(:,:) = tatm_bak
    emip(:,:) = emip_bak


  end subroutine get_monthly_forcing

  !! get levitus 3D  T and S for a certain month. These
  !! will be distributed by the C++ code and passed
  !! back into m_monthly by calling set_levitus
  !! 0<=month<=11 indicates which month you want.
  subroutine get_monthly_internal_forcing(ctemp,csalt,month)

    use, intrinsic :: iso_c_binding
    use m_lev

    implicit none

    real(c_double), dimension(n*m*l) :: ctemp,csalt

    integer(c_int), intent(in) :: month
    character(len=2) :: ibuf
    character(len=1024) :: fname
    real, dimension(n,m,l) :: temp_bak, salt_bak




    integer :: i,j,k,pos
    real :: dep

    ! keep annual mean data in backup locations
    temp_bak(:,:,:) = internal_temp
    salt_bak(:,:,:) = internal_salt

    write(ibuf,'(1I2.2)') month
    write(f99,*) 'reading internal seasonal temperature forcing from "'//topdir//'levitus/monthly/t'//ibuf//'an1'//'"'
    call levitus_internal(topdir//'levitus/monthly/t'//ibuf//'an1',internal_temp,.true.,'TEMP')
    write(f99,*) 'reading internal seasonal salinity forcing from "'//topdir//'levitus/monthly/s'//ibuf//'an1'//'"'
    call levitus_internal(topdir//'levitus/monthly/s'//ibuf//'an1',internal_salt,.true.,'SALT')

    ! note: the monthly mean data sets only contain temperature and salt up to a depth of
    ! 1500m. For any layers deeper than that we use the annual mean value instead.
    do k=l,1,-1
       dep=-z(k)*hdim
       if (dep.lt.depth(nlev_monthly)) then
          internal_temp(:,:,k)=temp_bak(:,:,k)
          internal_salt(:,:,k)=salt_bak(:,:,k)
       end if
    end do

    ! copy the 3D arrays to the 1D C++ arrays
    pos = 1
    do k=1,l
       do j=1,m
          do i=1,n
             ctemp(pos) = internal_temp(i,j,k)

             csalt(pos) = internal_salt(i,j,k)

             pos = pos+1
          end do
       end do
    end do

    ! restore annual mean data in original arrays
    internal_temp(:,:,:) = temp_bak
    internal_salt(:,:,:) = salt_bak

  end subroutine get_monthly_internal_forcing

  subroutine compute_flux(sol)

    implicit none

    real, dimension(n*m*l*nun), intent(in) :: sol
    real, dimension(n,m,l+la) :: T,S
    real, dimension(n,m) :: hf,fwf
    integer :: i,j,k,row
    real    :: temfun, salfun

    do k = 1, (l+la)
       do j = 1, m
          do i = 1, n
             row = find_row2(i,j,k,TT)
             T(i,j,k) = sol(row)
             row = find_row2(i,j,k,SS)
             S(i,j,k) = sol(row)
          enddo
       enddo
    enddo

    do j = 1,m
       do i = 1,n
          if (landm(i,j,l).eq.OCEAN) then
             hf(i,j)  = par(BIOT)*&
                  ( tatm(i,j) - T(i,j,l)/(par(TEMP)*par(COMB)) )
             fwf(i,j) = par(BIOT)*&
                  ( (1-its)*emip(i,j) + its*salfun(x(i),y(j))&
                  - S(i,j,l)/(par(SALT)*par(COMB)) )
          endif
       enddo
    enddo

    call write_forcing("flux.heat",  hf,14)
    call write_forcing("tatm.out", tatm,14)
    call write_forcing("flux.salt", fwf,15)
    call write_forcing("emip.out", emip,15)

  end subroutine compute_flux

end module m_global
