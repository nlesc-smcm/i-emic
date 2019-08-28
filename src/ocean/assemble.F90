!****************************************************************************
SUBROUTINE assemble
  !     assemble the global matrix A from the local matrices
  use m_mat
  use m_usr
  implicit none
  call TIMER_START('assemble' // char(0))
  call fillcolA
  call intcond

  call TIMER_STOP('assemble' // char(0))

end SUBROUTINE assemble

!****************************************************************************
! JT: I have un-commented the lines where the T and S points are set to -1,
! I think that they are required for the time integration scheme
SUBROUTINE fillcolB
  !     fill the columns of B
  USE m_mat
  use m_usr
  use m_atm
  implicit none
  integer find_row2
  integer i,j,k

  call TIMER_START('fillcolB' // char(0))

  ! Put B in coB,  B is a diagonal matrix.
  coB = 0.0
  do k = 1, l
     do j = 1, m
        do i = 1, n
           if ( landm(i,j,k) == OCEAN ) then
              if ( landm(i+1,j,k) /= LAND ) coB(find_row2(i,j,k,UU)) = -par(ROSB)
              if ( landm(i,j+1,k) /= LAND ) coB(find_row2(i,j,k,VV)) = -par(ROSB)
              !for testing periodic orbit solver, we re-scale this temporarily
              !if ( landm(i+1,j,k) /= LAND ) coB(find_row2(i,j,k,UU)) = -1.0
              !if ( landm(i,j+1,k) /= LAND ) coB(find_row2(i,j,k,VV)) = -1.0
              coB(find_row2(i,j,k,TT)) = -1.0
              coB(find_row2(i,j,k,SS)) = -1.0
           endif
        enddo
     enddo
  enddo
  
  if (rowintcon>0) then
     ! if(SRES == 0) coB(rowintcon + SS) = 0.0 !zero in B for integral condition
     if(SRES == 0) coB(rowintcon) = 0.0 !zero in B for integral condition
  end if

  call TIMER_STOP('fillcolB' // char(0))

end SUBROUTINE fillcolB

!****************************************************************************
SUBROUTINE fillcolA
  !     Fill the columns of A
  USE m_mat
  use m_usr
  implicit none
  integer find_row2
  integer i,j,k,ii,jj,kk,v,row,i2,j2,k2


  call TIMER_START('fillcolA' // char(0))
  !  +------------------------------------+
  !  |  stencil/neighbourhood             |
  !  | +----------++-------++----------+  |
  !  | | 12 15 18 || 3 6 9 || 21 24 27 |  |
  !  | | 11 14 17 || 2 5 8 || 20 23 26 |  |
  !  | | 10 13 16 || 1 4 7 || 19 22 25 |  |
  !  | |  below   || center||  above   |  |
  !  | +----------++-------++----------+  |
  !  |                                    |
  !  | Compressed sparse row format:      |
  !  |   co{.} : values                   |
  !  |  jco{.} : column indices           |
  !  |  beg{.} : row pointer              |
  !  +------------------------------------+
  !
  ! +-----------------------------------------------------------------------------+
  ! | Obtaining the matrix:                                                       |
  ! | 1) The top iteration travels through the grid points:                       |
  ! |     vertically       k = 1, l                                               |
  ! |     meridionally     j = 1, m                                               |
  ! |     zonally          i = 1, n                                               |
  ! |                                                                             |
  ! | 2) For each grid point we iterate over its neighbours: kk = 1, np           |
  ! |     Every neighbour's coordinate set is stored: (i2, j2, k2)                |
  ! |                                                                             |
  ! | 3) For each neighbour we iterate over every unknown: ii = 1, nun            |
  ! |     For every unknown we obtain its row and set a value in beg{.}.          |
  ! |     Assume an unknown is connected to itself and every other unknown in     |
  ! |     every point in the stencil, then the number of elements in each row is  |
  ! |     at most nun*np = 6*27. The row pointer for row r is therefore set to    |
  ! |     nun*np*(r-1) + 1. Later on this will be corrected... (?)                |
  ! |                                                                             |
  ! | 4) Then we iterate over every unknown including ii itself: jj = 1, nun      |
  ! |     The coefficient c in d/dt ii|(i,j,k) = c jj|(i2,j2,k2) + ...            |
  ! |     is stored in the row corresponding to ii|(i,j,k) and the column         |
  ! |     corresponding to jj|(i2,j2,k2).                                         |
  ! +-----------------------------------------------------------------------------+
  begA = 0
  v = 1
  row = 1
  do k = 1, l
     do j = 1, m
        do i = 1, n
           Alocal = An(:,:,:,i,j,k)
           do ii = 1, nun
              begA(row) = v
              do kk = 1,np
                 do jj = 1, nun
                    if (abs(Alocal(kk,ii,jj)).gt.1.0e-10) then
                       coA(v) = Alocal(kk,ii,jj)
                       ! shift(i,j,k,i2,j2,k2,kk) returns the neighbour at location kk
                       !  w.r.t. the center of the stencil (5) defined above.
                       !  it is faster to do this in here than outside of the loop.
                       call shift(i,j,k,i2,j2,k2,kk)
                       ! find_row2(i,j,k,ii) returns the row in the matrix for variable
                       !  ii at grid point (i,j,k) (matetc.F90)
                       jcoA(v) = find_row2(i2,j2,k2,jj)
                       v = v + 1
                    end if
                 end do
              end do
              row = row + 1
           end do

        end do
     end do
  end do

  ! final element of beg{.} array should be final row + 1
  begA(ndim + 1) = v

  call TIMER_STOP('fillcolA' // char(0))
end SUBROUTINE fillcolA

!****************************************************************************
SUBROUTINE shift(i,j,k,i2,j2,k2,kk)
  ! Defines location of neighbouring grid points
  ! +----------++-------++----------+
  ! | 12 15 18 || 3 6 9 || 21 24 27 |
  ! | 11 14 17 || 2 5 8 || 20 23 26 |
  ! | 10 13 16 || 1 4 7 || 19 22 25 |
  ! |  below   || center||  above   |
  ! +----------++-------++----------+
  ! shift(i,j,k,i2,j2,k2,kk) returns the neighbour at location kk
  !  w.r.t. the center of the stencil (5) defined above

  use m_usr
  implicit none
  integer i,j,k,i2,j2,k2,kk

  if (kk < 10) then
     k2 = k
     j2 = j - 1 + mod(kk+2,3)
     i2 = i - 1 + int((kk-1)/3)
  else if (kk < 19) then
     k2 = k - 1
     j2 = j - 1 + mod(kk+2,3)
     i2 = i - 1 + int((kk-10)/3)
  else
     k2 = k + 1
     j2 = j - 1 + mod(kk+2,3)
     i2 = i - 1 + int((kk-19)/3)
  end if

  if (periodic) then
     if (i2.eq.0) then
        i2 = n
     elseif (i2.eq.(n+1)) then
        i2 = 1
     endif
  endif

end SUBROUTINE shift

!****************************************************************************
SUBROUTINE intcond
  !     Impose integral condition
  USE m_mat
  use m_usr
  implicit none
  !      include 'mat.com'
  integer find_row2
  integer i, j, k, v, ic, shift, jcoIC(ndim)
  real    coIC(ndim)

  !     Replace P equation with a 'normalization' condition for p
  !      ic = rowintcon + PP
  !     call dirset(ic,0.0)

  !     Replace S equation with an 'integral' condition for s

  ! this is now done in Trilinos, so we simply !
  return

  if( SRES == 1 ) return
  if (rowintcon<0) return
  ic = rowintcon ! + SS
  call dirset(ic,0.0)

  v = 0
  do k = 1, l
     do j = 1, m
        do i = 1, n
           if ( landm(i,j,k) == OCEAN ) then
              v = v + 1
              jcoIC(v) = find_row2(i,j,k,SS)
              coIC(v) = cos(y(j)) * dfzT(k)
           endif
        enddo
     enddo
  enddo

  shift = v - begA(ic+1) + begA(ic)
  if ( shift > 0 ) then
     do i = ic + 1, ndim + 1
        begA(i) = begA(i) + shift
     enddo
     do i = begA(ndim+1) - 1, begA(ic+1), -1
        coA(i) =  coA(i - shift)
        jcoA(i) = jcoA(i - shift)
     enddo
  endif

  do i = 1, v
     j = begA(ic) + i - 1
     coA(j) =  coIC(i)
     jcoA(j) = jcoIC(i)
  enddo


end SUBROUTINE intcond

!****************************************************************************
SUBROUTINE dirset(row,Fd)
  !     inforce dirichlet condition
  USE m_mat
  use m_usr
  implicit none

  integer row, v
  real    Fd

  Frc(row)  =  Fd
  do v=begA(row),begA(row+1)-1
     coA(v)= 0.0
     if (jcoA(v).eq.row) coA(v) = 1.0
  enddo

end SUBROUTINE dirset

!****************************************************************************
