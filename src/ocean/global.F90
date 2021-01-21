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
       zmin, zmax,                           &
       hdim, qz,                             &
       itopo, flat, rd_mask,                 &
       coupled_T, coupled_S,                 &
       TRES, SRES, iza, ite, its, rd_spertm, &
       forcing_type,                         &
       rowintcon, f99, t0, s0

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
  real, dimension(:), pointer :: x,y,z,xu,yv,zw
  real :: dx, dy, dz
  integer, dimension(:,:,:),allocatable :: landm
  real, dimension(:,:), allocatable :: taux, tauy, tatm, emip, spert
  logical ::  periodic
  real, dimension(:,:,:), allocatable :: internal_temp,internal_salt

  character(len=999) :: maskfile, spertmaskfile
  character(len=999) :: windfile, sstfile, sssfile

contains

  !! allocate memory for the global data structures. You are allowed
  !! to specify m,n,l to be 0 so that no memory is allocated (most
  !! subdomains will only use xmin,xmax,ymin and ymax)
  subroutine initialize(a_n,a_m,a_l,&
       a_xmin,a_xmax,a_ymin,a_ymax,a_hdim,a_qz,&
       a_periodic,a_itopo,a_flat,a_rd_mask,&
       a_TRES,a_SRES,a_iza,a_ite,a_its,a_rd_spertm,&
       a_coupled_T, a_coupled_S,&
       a_forcing_type,&
       a_maskfile, a_spertmaskfile,&
       a_windfile, a_sstfile, a_sssfile)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int) :: a_n,a_m,a_l
    real(c_double) :: a_xmin,a_xmax,a_ymin,a_ymax,a_hdim,a_qz
    integer(c_int) :: a_periodic,a_itopo,a_flat,a_rd_mask
    integer(c_int) :: a_TRES,a_SRES,a_iza,a_ite,a_its,a_rd_spertm
    integer(c_int) :: a_coupled_T, a_coupled_S
    integer(c_int) :: a_forcing_type

    character(c_char), dimension(*) :: a_maskfile, a_spertmaskfile
    character(c_char), dimension(*) :: a_windfile, a_sstfile, a_sssfile

    xmin  = a_xmin
    xmax  = a_xmax
    ymin  = a_ymin
    ymax  = a_ymax
    dx = (xmax-xmin)/N
    dy = (ymax-ymin)/M
    dz = (zmax-zmin)/L
    hdim  = a_hdim
    qz    = a_qz

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
    forcing_type = a_forcing_type

    !===================== file names  =========================
    call set_filename(a_maskfile, maskfile)
    call set_filename(a_spertmaskfile, spertmaskfile)
    call set_filename(a_windfile, windfile)
    call set_filename(a_sstfile, sstfile)
    call set_filename(a_sssfile, sssfile)
    !===========================================================

    !========= I-EMIC coupling  ================================
    coupled_T = a_coupled_T
    coupled_S = a_coupled_S
    !===========================================================

    ! if (n/=a_n .or. m/=a_m .or. l/=a_l) then

    n = a_n
    m = a_m
    l = a_l
    ndim = n*m*l*nun

    call deallocate_global

    write(*,*) 'allocating fortran arrays'
    allocate(u(ndim), up(ndim), w(ndim,nf))
    allocate(landm(0:n+1,0:m+1,0:l+1))
    allocate(taux(n,m),tauy(n,m))
    allocate(tatm(n,m),emip(n,m),spert(n,m))
    allocate(internal_temp(n,m,l),internal_salt(n,m,l))
    internal_salt=0.0
    internal_temp=0.0
    landm = 0
  end subroutine initialize

  subroutine deallocate_global

    implicit none

    write(*,*) 'deallocating fortran arrays'

    if (allocated(u)) deallocate(u)
    if (allocated(up)) deallocate(up)
    if (allocated(w)) deallocate(w)

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

  subroutine set_filename(in_file, out_file)

    use, intrinsic :: iso_c_binding
    implicit none

    character(c_char), dimension(*) :: in_file
    character(*) out_file

    integer i, len, maxlen

    len = 1
    do while (in_file(len).ne.c_null_char)
       out_file(len:len) = in_file(len)
       len = len + 1
    end do
    maxlen = LEN_TRIM(out_file)
    do i=len,maxlen
       out_file(i:i) = ' '
    end do

  end subroutine set_filename

  subroutine set_maskfile(a_maskfile)

    use, intrinsic :: iso_c_binding
    implicit none

    character(c_char), dimension(*) :: a_maskfile

    call set_filename(a_maskfile, maskfile)

  end subroutine set_maskfile

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

  subroutine set_x(asize,ptr) bind(C,name="set_global_x")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), intent(in) :: asize
    real(c_double), dimension(asize), intent(in), target :: ptr

    x(1:asize)=>ptr
  end subroutine set_x

  subroutine set_y(asize,ptr) bind(C,name="set_global_y")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), intent(in) :: asize
    real(c_double), dimension(asize), intent(in), target :: ptr

    y(1:asize)=>ptr
  end subroutine set_y

  subroutine set_z(asize,ptr) bind(C,name="set_global_z")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), intent(in) :: asize
    real(c_double), dimension(asize), intent(in), target :: ptr

    z(1:asize)=>ptr
  end subroutine set_z

  subroutine set_xu(asize,ptr) bind(C,name="set_global_xu")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), intent(in) :: asize
    real(c_double), dimension(asize), intent(in), target :: ptr

    xu(0:asize-1)=>ptr
  end subroutine set_xu

  subroutine set_yv(asize,ptr) bind(C,name="set_global_yv")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), intent(in) :: asize
    real(c_double), dimension(asize), intent(in), target :: ptr

    yv(0:asize-1)=>ptr
  end subroutine set_yv

  subroutine set_zw(asize,ptr) bind(C,name="set_global_zw")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), intent(in) :: asize
    real(c_double), dimension(asize), intent(in), target :: ptr

    zw(0:asize-1)=>ptr
  end subroutine set_zw

  !! this is supposed to be called once by the root
  !! proc with a standard 0:n+1,0:m+1,0:l+1 1D C
  !! array which is then distributed to all the subdomains
  !! and re-inserted into the subdomains when calling usrc:init
  subroutine get_landm(cland)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), dimension((n+2)*(m+2)*(l+2)) :: cland

    integer :: i,j,k,pos


    _INFO_('THCM: global.F90 get_landm...')
    call topofit

    pos = 1
    do k=0,l+1
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

    integer(c_int), dimension((n+2)*(m+2)*(l+2)) :: cland

    integer :: i,j,k,pos

    _INFO_('THCM: global.F90 get_current_landm...')
    pos = 1
    do k=0,l+1
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

    integer(c_int), dimension((n+2)*(m+2)*(l+2)) :: cland
    integer :: i,j,k,pos


!    _INFO_('THCM: global.F90 set_landm...')

    landm = 0;

    pos = 1
    do k=0,l+1
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


    write(*,*) ' THCM: iza = ', iza

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

    if (coupled_T.eq.0) then
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
    else ! accepting  models within I-EMIC
       _INFO_(" Ocean: accepting external temperature, we are now part of the I-EMIC")
       _INFO_("   For now we put tatm = 0")
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

    write(f99,*) 'reading internal temperature forcing from "'//locate_file('levitus/new/t00an1')//'"'
    call levitus_internal(locate_file('levitus/new/t00an1'),internal_temp,'TEMP')
    !call levitus_internal(locate_file('levitus/new/avtemp'),internal_temp,'TEMP')

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

    if ((its.eq.0).and.(coupled_S.eq.0)) then
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
    write(f99,*) 'reading internal salt forcing from "'//locate_file('levitus/new/s00an1')//'"'
    call levitus_internal(locate_file('levitus/new/s00an1'),internal_salt,'SALT')
    !call levitus_internal(locate_file('levitus/new/avsalt'),internal_salt,'SALT')
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

  subroutine compute_flux(sol)

    implicit none

    real, dimension(n*m*l*nun), intent(in) :: sol
    real, dimension(n,m,l) :: T,S
    real, dimension(n,m) :: hf,fwf
    integer :: i,j,k,row
    real    :: salfun

    do k = 1, l
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
                  ( (1-its)*emip(i,j) + its*salfun(y(j))&
                  - S(i,j,l)/(par(SALT)*par(COMB)) )
          endif
       enddo
    enddo

    call write_forcing("flux.heat",  hf,14)
    call write_forcing("tatm.out", tatm,14)
    call write_forcing("flux.salt", fwf,15)
    call write_forcing("emip.out", emip,15)

  end subroutine compute_flux

  !------------------------------------------------------------------
  ! locate_file: finding files in rundir and topdir, return full path
  
  function locate_file(file) result(full)
    implicit none
    character(len=*)          :: file
    character(:), allocatable :: full
    integer                   :: status
    character(len=256)        :: msg
    
    open(unit=3511, file=file, status='old', iostat=status, iomsg=msg)
    close(3511)
    
    if (status.ne.0) then
       open(unit=3511, file=topdir//file, status='old', err=123, &
            iostat=status, iomsg=msg)
       close(3511)
       write(*,*) '  found ', file, ' in ', topdir
       full = topdir//file
    else
       write(*,*) '  found ', file, ' in rundir'
       full = file
    endif
    
    return
    
123 write(*,*) 'WARNING: failed to open: ', file
    write(*,*) ' topdir: ', trim(topdir)
    write(*,*) '    msg: ', trim(msg)
    write(*,*) ' iostat: ', status
    call throw_error('failed to open file')
    
  end function locate_file
  
  !-----------------------------------------------------------------------------
  ! let the cpp code throw an error for us
  subroutine throw_error(msg)
    implicit none
    character(len=*), intent(in) :: msg
    external thcm_throw_error

    call thcm_throw_error(msg)
  end subroutine throw_error

end module m_global
