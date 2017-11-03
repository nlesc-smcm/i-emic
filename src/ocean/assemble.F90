!****************************************************************************
SUBROUTINE assemble
  !     assemble the global matrix A from the local matrices
  use m_mat
  use m_usr
  implicit none
  call preprocessA_old !--> is this always necessary?
  call fillcolA
  call intcond

  call packA

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
  integer i,j,k,row

  !     Put B in coB,  B is a diagonal matrix.
  !     PRINT *,'parROSB = ',par(ROSB)
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
  if(la == 1) then
     do j = 1, m
        do i = 1, n
           coB(find_row2(i,j,l+la,TT)) = -Ai
        enddo
     enddo
  endif
  if (rowintcon>0) then
     ! if(SRES == 0) coB(rowintcon + SS) = 0.0 !zero in B for integral condition
     if(SRES == 0) coB(rowintcon) = 0.0 !zero in B for integral condition
  end if

end SUBROUTINE fillcolB

!****************************************************************************
SUBROUTINE preprocessA
  USE m_mat
  use m_usr
  implicit none
  integer i,j,k,ii,jj,kk, nnz, locnnz, count

  active = .false.

  i = int(n/2)
  j = int(m/2)
  k = int(l/2)

  if (landm(i,j,k) == OCEAN) then
     do ii = 1,nun
        do jj = 1,nun
           do kk = 1,np
              !                  if (abs(al(i,j,k,kk,ii,jj)) > 1.0E-14) then
              if (al(i,j,k,kk,ii,jj) /= 0D0) then
                 active(kk,ii,jj) = .true.
              end if
           end do
        end do
        active(5,ii,ii) = .true.
     end do
  else
     call preprocessA_old
  end if
end SUBROUTINE preprocessA

!****************************************************************************
SUBROUTINE preprocessA_old
  USE m_mat
  use m_usr
  implicit none
  integer i,j,k,ii,jj,kk, nnz, locnnz, count

  active = .false.
  count = 0
  do ii = 1,nun
     do jj = 1,nun
        do kk = 1,np
           loop: do i = 1,n
              do j = 1,m
                 do k = 1,l+la
                    ! --> 1.0E-10 small enough?
                    if (abs(al(i,j,k,kk,ii,jj)) > 1.0E-10) then
                       active(kk,ii,jj) = .true.
                       cycle loop
                    end if
                 end do
              end do
           end do loop
           !         if (active(kk,ii,jj)) then
           !            count = count +1
           !         end if
        end do
     end do
  end do
  !      write(f99,*) "total nnz in A:", nnz
  !      write(f99,*) "with ", count, " active connections out of ", nun*nun*np
#ifdef DEBUGGING
  ! this can be useXul for determining the maximal matrix graph
  open(42,file='active.txt')
  do ii=1,nun
     do jj=1,nun
        do kk=1,np
           if (active(kk,ii,jj)) then
              write(42,*) ii,jj,kk
           end if
        end do
     end do
  end do
  close(42)
#endif
end SUBROUTINE preprocessA_old


!****************************************************************************
SUBROUTINE fillcolA
  !     Fill the columns of A
  USE m_mat
  use m_usr
  implicit none
  integer find_row2
  integer i,j,k,ii,jj,kk,v,w,is,js,ks,bis,row,col,i2,j2,k2

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
  ! |     vertically       k = 1, l+la                                            |
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
  coA = 0.0
  is  = nun        !--> any use?
  js  = nun*n      !--> any use?
  ks  = nun*n*m    !--> any use?
  do k = 1, l+la
     do j = 1, m
        do i = 1, n

           do kk = 1,np

              ! shift(i,j,k,i2,j2,k2,kk) returns the neighbour at location kk
              !  w.r.t. the center of the stencil (5) defined above
              call shift(i,j,k,i2,j2,k2,kk)
              do ii = 1, nun

                 ! find_row2(i,j,k,ii) returns the row in the matrix for variable
                 !  ii at grid point (i,j,k) (matetc.F90)
                 row = find_row2(i,j,k,ii)
                 v   = nun*np*(row-1)
                 begA(row) = v + 1   ! assuming every row will have nun*np entries
                 do jj = 1, nun
                    if (active(kk,ii,jj)) then
                       w = v + (jj-1)*np ! ? grouped ordering?                       
                       coA(w+kk)  = al(i,j,k,kk,ii,jj)
                       jcoA(w+kk) = find_row2(i2,j2,k2,jj)
                    end if
                 end do
              end do
           end do

        end do
     end do
  end do

  ! final element of beg{.} array should be final row + 1
  begA(ndim + 1) = nun * np * ndim + 1

  ! --> what happens here?
  if (periodic) then
     bis= nun*(n-1)
     do k = 1, l+la
        do j = 1, m
           do ii = 1, nun
              row = find_row2(1,j,k,ii)
              v = np * nun * (row - 1)
              do jj = 1, nun
                 w = v + (jj-1) * np
                 jcoA(w+1) = find_row2(n,j-1,k  ,jj)
                 jcoA(w+2) = find_row2(n,j  ,k  ,jj)
                 jcoA(w+3) = find_row2(n,j+1,k  ,jj)
                 jcoA(w+10)= find_row2(n,j-1,k-1,jj)
                 jcoA(w+11)= find_row2(n,j  ,k-1,jj)
                 jcoA(w+12)= find_row2(n,j+1,k-1,jj)
                 jcoA(w+19)= find_row2(n,j-1,k+1,jj)
                 jcoA(w+20)= find_row2(n,j  ,k+1,jj)
                 jcoA(w+21)= find_row2(n,j+1,k+1,jj)
              enddo
           enddo
           do ii = 1, nun
              row = find_row2(n,j,k,ii)
              v = np * nun * (row - 1)
              do jj = 1, nun
                 w = v + (jj - 1) * np
                 jcoA(w+7) = find_row2(1,j-1,k  ,jj)
                 jcoA(w+8) = find_row2(1,j  ,k  ,jj)
                 jcoA(w+9) = find_row2(1,j+1,k  ,jj)
                 jcoA(w+16)= find_row2(1,j-1,k-1,jj)
                 jcoA(w+17)= find_row2(1,j  ,k-1,jj)
                 jcoA(w+18)= find_row2(1,j+1,k-1,jj)
                 jcoA(w+25)= find_row2(1,j-1,k+1,jj)
                 jcoA(w+26)= find_row2(1,j  ,k+1,jj)
                 jcoA(w+27)= find_row2(1,j+1,k+1,jj)
              enddo
           enddo
        enddo
     enddo
  endif

end SUBROUTINE fillcolA

!****************************************************************************
SUBROUTINE shift(i,j,k,i2,j2,k2,np)
  ! Defines location of neighbouring grid points
  ! +----------++-------++----------+
  ! | 12 15 18 || 3 6 9 || 21 24 27 |
  ! | 11 14 17 || 2 5 8 || 20 23 26 |
  ! | 10 13 16 || 1 4 7 || 19 22 25 |
  ! |  below   || center||  above   |
  ! +----------++-------++----------+
  ! shift(i,j,k,i2,j2,k2,np) returns the neighbour at location np
  !  w.r.t. the center of the stencil (5) defined above

  implicit none
  integer i,j,k,i2,j2,k2,np

  if (np < 10) then
     k2 = k
     j2 = j - 1 +  mod(np+2,3)
     i2 = i - 1 + int((np-1)/3)
  else if (np < 19) then
     k2 = k - 1
     j2 = j - 1 +  mod(np+2,3)
     i2 = i - 1 + int((np-10)/3)
  else
     k2 = k + 1
     j2 = j - 1 +  mod(np+2,3)
     i2 = i - 1 + int((np-19)/3)
  end if
end SUBROUTINE shift

!****************************************************************************
SUBROUTINE intcond_old
  !     Impose integral condition
  USE m_mat
  use m_usr
  implicit none
  !      include 'mat.com'
  integer find_row2
  integer i, j, k, v, L1, row
  !     Replace p equation at part. point with
  !     a 'normalization' condition for p
  !     L1 = nun*((L/2-1)*N*M+ N*(M/2-1) + N/2-1) + PP
  L1 = find_row2(n,m,l,PP)
  do v = begA(L1),begA(L1+1)-1
     coA(v)= 0.0
     if (jcoA(v).eq.L1) then
        coA(v) = 1.0
     endif
  enddo
  Frc(L1) = 0.0
  !     Replace S equation at ndim with an 'integral' condition for s
  do v = begA(ndim),begA(ndim+1)-1
     coA(v) = 0.0
  enddo

  v = begA(ndim)
  do k = 1, l
     do j = 1, m
        do i = 1, n
           if ( landm(i,j,k) == OCEAN ) then
              jcoA(v)= find_row2(i,j,k,SS)
              coA(v) = cos(y(j)) * dfzT(k)
              v = v+1
           endif
        enddo
     enddo
  enddo
  begA(ndim+1) = v
  Frc(ndim) =  0.0

end SUBROUTINE intcond_old
!****************************************************************************
SUBROUTINE intcond
  !     Impose integral condition
  USE m_mat
  use m_usr
  implicit none
  !      include 'mat.com'
  integer find_row2
  integer i, j, k, v, ic, row, shift, jcoIC(ndim)
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
SUBROUTINE packA
  !     remove zero entries in A
  USE m_mat
  use m_usr
  implicit none

  integer vn,v,begin,i,oldsize
  integer ii,jj,kk,XX,i2,j2,k2

  vn = 1
  oldsize = begA(ndim+1) - 1
  do i = 1, ndim
     begin = vn
     do v = begA(i), begA(i+1) - 1
        if (abs(coA(v)).gt.1.0e-15) then
           coA(vn)  =  coA(v)
           jcoA(vn) = jcoA(v)
           vn = vn + 1
           if ((jcoA(v).gt.ndim).or.(jcoA(v).lt.0)) then
              write(*,*) 'row = ', i, 'col = ', jcoA(v), 'entry = ', coA(v)
              call findex(i,ii,jj,kk,XX)
              write(*,*) 'row grid point (i,j,k):',ii,jj,kk
              write(*,*) 'variable',XX
              vn = mod(v-1,27)+1
              call shift(ii,jj,kk,i2,j2,k2,vn)
              write(*,*) 'col grid point (i,j,k):',i2,j2,k2
              !                  call findex(jcoA(v),ii,jj,kk,XX)
              !                  write(*,*) ii,jj,kk,XX
              stop 'in packA: index out of range'
           endif
        endif
     enddo
     begA(i) = begin
  enddo
  begA(ndim+1) = vn
  !     write(f99,999) oldsize,begA(ndim+1)-1

999 format('size before and after packing: ',2i8)

end SUBROUTINE packA

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
