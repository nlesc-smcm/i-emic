#include "fdebug.h"
! this module has been made 'global' i.e. it uses m_global instead
! of m_global and should be used only by the root process in case
! Trilinos-THCM is run in parallel
!***************************************************************
SUBROUTINE topofit
  use m_global
  implicit none

  real       depth(n,m),dav,hbs(n,m)

  integer,parameter :: nd=6*360, md=6*180
  integer :: lwrk, liwrk

  integer*2 d10(nd,md)
  real      rd(md,nd),xd(nd),yd(md),dxd,dyd,wrk2(4*n+nd+4)
  integer   ifail,px,py,iwrk(n+nd)
  real      lambda(nd+4),mu(md+4),cspl(nd*md),wrk((nd+6)*(md+6))
  real      dumd(n*m), xi(n*m), yi(n*m)

  integer    i,j,io,sfend,sfstart

  _DEBUG_('entering topofit...')
  
  lwrk=4*n+nd+4
  liwrk=n+nd


  depth = 0

  IF (rd_mask) THEN
     _DEBUG_('call readmask...')
     call readmask
  ELSE
     _DEBUG_('call depth3land...')
     call depth3land(depth)
  ENDIF
  _DEBUG_('leave topofit')
END SUBROUTINE topofit
!*****************************************************************
SUBROUTINE readmask
  use m_global
  implicit none
  integer i,j,k
  integer nw,nnw,nsw,ns,ne,nne,nse,nn,nsum,status

  write(*,*) '===========TOPOGRAPHY==================================='
  write(*,*) 'land mask is read in from file: mkmask/'//trim(maskfile)

  open(unit=500, file=locate_file('mkmask/'//trim(maskfile)), &
       status='old', err=123)

  write(*,*) '===========TOPOGRAPHY==================================='

  landm = LAND
  
  do k = 0, l+1
     read(500,*)
     do j = m+1, 0, -1
        read(500,'(362i1)',iostat=status) (landm(i,j,k),i=0,n+1)
     enddo
  enddo

  close(500)


  ! do k = 1, l
  !    do i = 3, n-2
  !       do j = 2, m-1
  !          if ( landm(i,j,k).eq.OCEAN ) then
  !             nw  =  landm(i-1,j,k)
  !             nnw =  landm(i-1,j+1,k)
  !             nsw =  landm(i-1,j-1,k)
  !             ne  =  landm(i+1,j,k)
  !             nne =  landm(i+1,j+1,k)
  !             nse =  landm(i+1,j-1,k)
  !             ns  =  landm(i,j-1,k)
  !             nn  =  landm(i,j+1,k)

  !             nsum = nw+ne+ns+nn
  !             if (nsum.gt.2) then
  !                landm(i,j,k) = LAND ! disabled this hack?
  !             endif

  !             if (((nw + ns) .eq. 2).and. (nsw.eq.OCEAN)) then
  !                landm(i,j,k) = LAND ! disabled this hack?
  !             endif

  !          endif
  !       enddo
  !    enddo
  ! enddo
  !
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

  if( FLAT ) then ! remove topography
     do k = 1, l-1
        landm(:,:,k) = landm(:,:,l)
     enddo
  endif

  i = n/2
  j = 6*m/8
  k = l
  ! integral-condition not implemented in parallel, -1 disables it completely:
  !rowintcon = find_row2(i,j,k,SS)
  rowintcon=-1

  ! write(*,*) '____',i,j,k, rowintcon,'____'
  ! do k = 1, l
  !    write(*,*) '______________',z(k)*hdim,'__________________'
  !    do j = m+1, 0, -1
  !       write(*,'(92i1)') landm(:,j,k)
  !    enddo
  ! enddo

  return

123 write(*,*) 'WARNING: failed to read land mask from file'
  write(*,*) '         mkmask/', trim(maskfile), '.'
  write(*,*) '         Continuing with no land.'
  
  call throw_error('failed to read land mask')

END SUBROUTINE readmask
!*****************************************************************
SUBROUTINE depth3land(depth)
  use m_global
  implicit none
  real    depth(n,m)
  !w
  logical atl,iatl
  logical mam,imam
  !w
  integer i, j, k

  real     ph1,ph2,ph3,ph4,thn,tha,thsa,thd

  _DEBUG_('entering depth3land...')
  
  depth = depth / hdim

  landm = LAND
  do k = 1, l
     do j = 1, m
        do i = 1, n
           if ( z(k) > depth(i,j) ) landm(i,j,k) = WATER
        enddo
     enddo
  enddo

  SELECT CASE(itopo)
  CASE(0) ! topography from data
     i = n/2
     j = m/2
     k = l
     do while ( landm(i,j,k) /= WATER )
        i = i + 1
        if(i > n) stop 'in depth3land: cannot find ocean point'
     enddo
     call fillbays_old
     call flood(i,j,k,WATER,OCEAN)
     where( landm == WATER ) landm = LAND
     !rowintcon = find_row2(i,j,k,SS)
     rowintcon = -1
     write(*,*) 'integral condition ',i,j,k,landm(i,j,k),rowintcon
     ! END CASE 0
  CASE(1) ! no continents
     _DEBUG_('case 1: no continents')
     i = n
     j = m
     k = l
     !rowintcon = find_row2(i,j,k,SS)
     rowintcon = -1
     write(*,*) 'integral condition ',i,j,k,landm(i,j,k),rowintcon
     landm(1:n,1:m,1:l) = OCEAN
     !         landm(7:9,1:m,1:5) = LAND
     ! END CASE 1
  CASE(2) ! Miocene
     i = n-1
     j = m/2
     k = l-1
     !rowintcon = find_row2(i,j,k,SS)
     rowintcon = -1
     write(*,*) 'integral condition ',i,j,k,landm(i,j,k),rowintcon
     landm(1:n,1:m,1:l) = OCEAN
     !
     ph1 = 250*pi/180.
     ph2 = 315*pi/180.
     ph3 = 10*pi/180.
     ph4 = 65.*pi/180.
     thd = -60*pi/180.
     thsa = -35*pi/180.
     thn = 10.*pi/180.
     tha =  30*pi/180.
     !
     !     south america
     !
     DO i=1,n
        IF ((x(i).lt.ph2).and.(x(i).gt.ph1)) THEN
           DO j = 1,m
              IF ((y(j).lt.0.).and.(y(j).gt.thd)) THEN
                 DO k = 1,l
                    landm(i,j,k) = LAND
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     !
     !     south africa
     !
     DO i=1,n
        IF ((x(i).lt.ph4).and.(x(i).gt.ph3)) THEN
           DO j = 1,m
              IF ((y(j).lt.thn).and.(y(j).gt.thsa)) THEN
                 !              IF ((y(j).lt.0.).and.(y(j).gt.thd)) THEN
                 DO k = 1,l
                    landm(i,j,k) = LAND
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     !
     !    north america
     !
     DO i=1,n
        IF ((x(i).lt.ph2).and.(x(i).gt.ph1)) THEN
           DO j = 1,m
              IF ((y(j).lt.ymax).and.(y(j).gt.tha)) THEN
                 DO k = 1,l
                    landm(i,j,k) = LAND
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     !
     !    asia
     !
     DO i=1,n
        IF ((x(i).lt.ph4).and.(x(i).gt.ph3)) THEN
           DO j = 1,m
              IF ((y(j).lt.ymax).and.(y(j).gt.tha)) THEN
                 DO k = 1,l
                    landm(i,j,k) = LAND
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     ! END CASE 2
  CASE(3) !DB - SH case
     i = n
     j = m
     k = l
     !rowintcon = find_row2(i,j,k,SS)
     rowintcon = -1
     write(*,*) 'integral condition ',i,j,k,landm(i,j,k),rowintcon
     landm(1:n,1:m,1:l) = OCEAN
     !     landm(18:20,3:14,1:l) = LAND
     landm(18:20,1:16,1:l) = LAND
     ! END CASE 3
  CASE(4) !DB - DH case
     i = n
     j = m
     k = l
     !rowintcon = find_row2(i,j,k,SS)
     rowintcon=-1
     write(*,*) 'integral condition ',i,j,k,landm(i,j,k),rowintcon
     landm(1:n,1:m,1:l) = OCEAN
     !     landm(1:1,1:m,1:l) = LAND
     !     landm(n:n,1:m,1:l) = LAND
     landm(22:24,6:m,1:l) = LAND
     ! resolution 40 x 30 x 16 , domain 100-260 x -60 - 60
     ! END CASE 4
  END SELECT

  if( FLAT ) then ! remove topography
     do k = 1, l-1
        landm(:,:,k) = landm(:,:,l)
     enddo
  endif

  !w
  if (.false.) then
     do j = 1, m
        do i = 1, n
           iatl = atl(x(i),y(j))
           if (.not.iatl) landm(i,j,:) = LAND
        enddo
     enddo
  endif
  if (.false.) then
     do j = 1, m
        do i = 1, n
           imam = mam(x(i),y(j))
           if (imam) landm(i,j,:) = LAND
        enddo
     enddo
  endif
  !w

  if (periodic) then
     where( landm(1,:,:) == OCEAN .and. landm(n,:,:) == OCEAN )
        landm(n+1,:,:) = PERIO
        landm(  0,:,:) = PERIO
     end where
  endif
  
  write(*,*) '===========TOPOGRAPHY=================='
  write(*,10) itopo
  write(*,*) 'land mask is written to unit 77'
  write(*,*) '===========TOPOGRAPHY=================='

  open(77,FILE=rundir//'fort.77')
  do k = 0, l+1
     write(77,*) "-------------------------------------"
     do j = m+1, 0, -1
        write(77,'(362i1)') landm(:,j,k)
     enddo
  enddo
  close(77)

10 format(1x,'you have chosen topography option',1x,i8)
end SUBROUTINE depth3land

!*****************************************************************
recursive subroutine flood( i, j, k, old, new )
  use m_global
  implicit none
  integer i, j, k, old, new
  if( landm(i,j,k) == old ) then
     landm(i,j,k) = new
     call flood(i+1,j  ,k  , old, new)
     call flood(i-1,j  ,k  , old, new)
     call flood(i  ,j+1,k  , old, new)
     call flood(i  ,j-1,k  , old, new)
     call flood(i  ,j  ,k+1, old, new)
     call flood(i  ,j  ,k-1, old, new)
  endif
  _DEBUG_('leaving depth3land')
end subroutine flood

!*****************************************************************
subroutine fillbays
  use m_global
  implicit none

  integer it, i, j, k, ns, ew
  integer nedit,oldland(0:n+1,0:m+1,0:l+1)

  do it = 1, 15
     oldland = landm
     do k = 1, l
        do j = 1, m
           do i = 1, n
              ns = 0
              ew = 0
              if ( landm(i+1,j  ,k) == LAND ) ew = ew + 1
              if ( landm(i-1,j  ,k) == LAND ) ew = ew + 1
              if ( landm(i  ,j+1,k) == LAND ) ns = ns + 1
              if ( landm(i  ,j-1,k) == LAND ) ns = ns + 1
              if ( (ew+ns).ge.3 ) landm(i,j,k) = LAND
              if ( (ew+ns).le.1 ) landm(i,j,k) = OCEAN
           enddo
        enddo
     enddo
     do j = 1, m
        do i = 1, n
           if ( landm(i,j,l-1) == LAND ) landm(i,j,l) = LAND
        enddo
     enddo
     !
     nedit = 0
     do k = 1, l
        do j = 1, m
           do i = 1, n
              nedit = nedit + abs(landm(i,j,k)-oldland(i,j,k))
           enddo
        enddo
     enddo
     write(*,*)'fillbays, it : ',it,' number of edits: ',nedit
     if (nedit.eq.0) goto 100
  enddo
100 continue
  write(*,*)'fillbays: ',it-1,' iterations'

end subroutine fillbays
!*****************************************************************
subroutine fillbays_old
  use m_global
  implicit none

  integer it, i, j, k, ns, ew, ne, se
  integer nedit,oldland(0:n+1,0:m+1,0:l+1)

  do it = 1, 40
     oldland = landm
     do k = 1, l
        do j = 1, m
           do i = 1, n
              ns = 0
              ew = 0
              ne = 0
              se = 0
              if ( landm(i+1,j  ,k) == LAND ) ew = ew + 1
              if ( landm(i-1,j  ,k) == LAND ) ew = ew + 1
              if ( landm(i  ,j+1,k) == LAND ) ns = ns + 1
              if ( landm(i  ,j-1,k) == LAND ) ns = ns + 1
              if ( landm(i+1,j+1,k) == LAND ) ne = ne + 1
              if ( landm(i-1,j-1,k) == LAND ) ne = ne + 1
              if ( landm(i-1,j+1,k) == LAND ) se = se + 1
              if ( landm(i+1,j-1,k) == LAND ) se = se + 1
              if ( ew==2.or.ns==2.or.ne==2.or.se==2 ) landm(i,j,k) = LAND
              if ( ew==0.or.ns==0.or.ne==0.or.se==0 ) landm(i,j,k) = OCEAN
           enddo
        enddo
     enddo
     do j = 1, m
        do i = 1, n
           if ( landm(i,j,l-1) == LAND ) landm(i,j,l) = LAND
        enddo
     enddo
     !
     nedit = 0
     do k = 1, l
        do j = 1, m
           do i = 1, n
              nedit = nedit + abs(landm(i,j,k)-oldland(i,j,k))
           enddo
        enddo
     enddo
     write(*,*)'fillbays, it : ',it,' number of edits: ',nedit
     if (nedit.eq.0) goto 100
  enddo
100 continue
  write(*,*)'fillbays: ',it-1,' iterations'

end subroutine fillbays_old
!******************************************************************
FUNCTION atl(xin,yin)
  implicit none
  !
  !     INPUT/OUTPUT
  logical atl
  real xin,yin
  !     LOCAL
  real x,y,pi
  real x_hud,x_nam,x_sam,x_eur,x_med,x_afr,x_mam
  !
  pi =  3.14159265358979323846
  x  = xin * 180/pi  ! x E [0,360]
  y  = yin * 180/pi  ! y E [-90,90]
  if (x.lt.0.)   x = x+360.
  if (x.gt.360.) x = x-360.
  !
  x_hud = 285.         ! Hudson Strait          (75W)
  x_nam = 260.         ! North America         (100W)
  x_sam = 290.         ! South America          (70W)
  x_eur = 30.          ! Europe                 (30E)
  x_med = 40.          ! Mediterranean          (40E)
  x_afr = 20.          ! Africa                 (20E)
  x_mam = 290.-1.5*y   ! Middle American land bridge
  !
  atl = .true.
  if ((y.ge.-55).and.(y.lt.0)) then
     !     if (y.lt.0) then
     if  ((x.ge.x_afr).and.(x.le.x_sam)) atl = .false.
  elseif ((y.ge.0).and.(y.lt.20)) then
     if  ((x.ge.x_afr).and.(x.le.x_mam)) atl = .false.
  elseif ((y.ge.20).and.(y.lt.30)) then
     if  ((x.ge.x_afr).and.(x.le.x_nam)) atl = .false.
  elseif ((y.ge.30).and.(y.lt.50)) then
     if  ((x.ge.x_med).and.(x.le.x_nam)) atl = .false.
  elseif (y.ge.50) then
     if  ((x.ge.x_eur).and.(x.le.x_hud)) atl = .false.
  endif
  !
  return
end FUNCTION atl
!******************************************************************
FUNCTION mam(xin,yin)
  implicit none
  !
  !     INPUT/OUTPUT
  logical mam
  real xin,yin
  !     LOCAL
  real x,y,pi
  real x_mam_w,x_mam_e
  !
  pi =  3.14159265358979323846
  x  = xin * 180/pi  ! x E [0,360]
  y  = yin * 180/pi  ! y E [-90,90]
  if (x.lt.0.)   x = x+360.
  if (x.gt.360.) x = x-360.
  !
  x_mam_w = 285.-1.5*y   ! Middle American land bridge
  x_mam_e = 295.-1.5*y   ! Middle American land bridge
  !
  mam = .false.
  if ((y.ge.0).and.(y.lt.30)) then
     if  ((x.ge.x_mam_w).and.(x.le.x_mam_e)) mam = .true.
  endif
  !
  return
end FUNCTION mam

