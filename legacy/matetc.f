********************************************************************
      SUBROUTINE initw(w,ndim,nf)
      implicit none
*     IMPORT/EXPORT
      integer  ndim,nf
      real     w(ndim,nf)
*     LOCAL
      integer  i,l
*     FUNCTIONS
      real    g05caf
      real    rand
*
      DO i=1,ndim
       DO l=1,nf
!        w(i,l) = g05caf()
 !        w(i,l) = rand()
       ENDDO
      ENDDO
*
      END
*****************************************************************************
      real FUNCTION l2nrm(f,ndim)
      implicit none
      integer i,ndim
      real    f(ndim),nrm
*
      nrm = 0.0
      DO i=1,ndim
       nrm = f(i)*f(i) + nrm
      ENDDO
      l2nrm = sqrt(nrm)
*
      END
*****************************************************************************
      real FUNCTION linrm(f,ndim)
      implicit none
      integer i,ndim
      real    f(ndim),nrm
*
      nrm = 0.0
      DO i=1,ndim
       IF ( abs(f(i)).GT.nrm ) nrm = abs(f(i))
      ENDDO
      linrm = nrm
*
      END
****************************************
      subroutine findex(row,i,j,k,XX)
      ! find variable XX and grid point (i,j,k) that
      ! correspond to row 'row' in matrix A
      implicit none
      include 'usr.com'
      integer row, i, j, k, XX
      IF (prec.eq.1) THEN
         call findex_interleaved_ordering(row,i,j,k,XX)
      ELSE
         call findex_block_ordering(row,i,j,k,XX)
      END IF
      end subroutine findex

****************************************
      subroutine findex_block_ordering(row,i,j,k,XX)
      ! find variable XX and grid point (i,j,k) that
      ! correspond to row 'row' in matrix A
      implicit none
      include 'usr.com'
      integer t, r, XX, row, i, j, k, nn
      nn = n*m*(l+la)
      if (row.LE.2*nn) then
         XX = 2 - mod(row,2) ! UU or VV
         r = (row + 2 - XX)/2
      elseif (row.LE.3*nn) then
         XX = WW
         r = mod(row-1,nn)+1
      elseif (row.LE.4*nn) then
         XX = PP
         r = mod(row-1,nn)+1
      elseif (row.LE.6*nn) then
         XX = 6 - mod(row,2) ! TT or SS
         r = (mod(row-1,2*nn)+1 +6 - XX)/2
      else
         STOP "error in findex (matetc.f): row > ndim " 
      end if
      k = int((r-1)/(n*m))+1
      r = mod(r-1,n*m)+1
      j = int((r-1)/n)+1
      i = mod(r-1,n)+1
      end subroutine findex_block_ordering

****************************************
      subroutine findex_interleaved_ordering(row,i,j,k,XX)
      ! find variable XX and grid point (i,j,k) that
      ! correspond to row 'row' in matrix A
      implicit none
      include 'usr.com'
      integer t, r, XX, row, i, j, k

      XX = mod(row - 1, nun) + 1
      r = (row - XX)/nun
      t = mod(r, n*m)
      k = (r - t)/(n*m) + 1
      i = mod(t, n) + 1
      j = (t - i + 1)/n + 1
      end subroutine findex_interleaved_ordering

****************************************
      integer function find_row(i,j,k)
      implicit none
      include 'par.com'
      integer i,j,k
      find_row = nun*((k-1)*n*m+ n*(j-1) + i-1)
      WRITE(99,*) "use of find_row is prohibited, use find_row2 instead"
      STOP
      end function find_row

****************************************
      integer function find_row2(i,j,k,XX)
      ! find row in matrix A of variable XX at grid point (i,j,k)
      implicit none
!      include 'par.com'
      include 'usr.com'
      integer i,j,k,XX
!     functions
      integer find_row_interleaved_ordering, find_row_block_ordering
      IF (prec.eq.1) THEN
         find_row2 = find_row_interleaved_ordering(i,j,k,XX)
      ELSE
         find_row2 = find_row_block_ordering(i,j,k,XX)
      END IF
      end function find_row2

****************************************
      integer function find_row_block_ordering(i,j,k,XX)
      ! find row in matrix A of variable XX at grid point (i,j,k)
      implicit none
      include 'usr.com'
      integer i,j,k,XX
      select case (XX)
      case (UU:VV)
         find_row_block_ordering = 2*((k-1)*n*m + (j-1)*n + i-1) + XX
      case (WW)
         find_row_block_ordering = (k-1)*n*m + (j-1)*n + i + (WW-1)*n*m*(l+la)
      case (PP)
         find_row_block_ordering = (k-1)*n*m + (j-1)*n + i + (PP-1)*n*m*(l+la)
      case (TT:SS)
         find_row_block_ordering = 2*((k-1)*n*m + (j-1)*n + i-1) +(TT-1)*n*m*(l+la)+(XX-4)
      case default
         find_row_block_ordering = -1
      end select
      ! check dimensions
      if ((i>n).or.(i<1).or.(j>m).or.(j<1).or.(k>l+la).or.(k<1)) then
         find_row_block_ordering = -1
      end if
      end function find_row_block_ordering

****************************************
      integer function find_row_interleaved_ordering(i,j,k,XX)
      ! find row in matrix A of variable XX at grid point (i,j,k)
      implicit none
      include 'usr.com'
      integer i,j,k,XX
      find_row_interleaved_ordering = nun*((k-1)*n*m+ n*(j-1) + i-1)+XX
      ! check dimensions
      if ((i>n).or.(i<1).or.(j>m).or.(j<1).or.(k>l+la).or.(k<1)) then
         find_row_interleaved_ordering = -1
      end if
      end function find_row_interleaved_ordering

*******************************************************************************
      SUBROUTINE matAvec(v1,v2)
*     This multiplies sparse matrix A and vector v1 to vector v2
      USE m_mat
      implicit none
      include 'par.com'
!      include 'mat.com'
      real     v1(ndim),v2(ndim)
*     LOCAL
      integer  i,v
*
      v2 = 0.0
      DO i = 1,ndim
       DO v = begA(i),begA(i+1)-1
         v2(i) = coA(v)*v1(jcoA(v)) + v2(i)
       ENDDO
      ENDDO
*
      END
*******************************************************************************
      SUBROUTINE matBvec(v1,v2)
*     This multiplies sparse matrix B and vector v1 to vector v2
*     B is a diagonal matrix
      USE m_mat
      implicit none
      include 'par.com'
!      include 'mat.com'
      real     v1(ndim),v2(ndim)
*     LOCAL
      integer  i
*
      DO i = 1,ndim
       v2(i) = coB(i)*v1(i)
      ENDDO
*
      END
************************************************************************
      SUBROUTINE writematrhs(rl) 
*     write the matrix A + rhs to a file
      USE m_mat
      implicit none
      include 'par.com'
!      include 'mat.com'
*     INPUT/OUTPUT
      real     rl(ndim) 
*     LOCAL
      integer i
*
      open( 9,file='B.co')     
      open(10,file='A.info')     
      open(11,file='A.beg')
      open(12,file='A.jco')
      open(13,file='A.co')
      open(14,file='A.rl')     
      write(10,*) ndim,begA(ndim+1)-1      
      DO i=1,ndim+1
        write(11,*) begA(i)
      ENDDO
      DO i=1,begA(ndim+1)-1
        write(12,*) jcoA(i)
      ENDDO
      DO i=1,begA(ndim+1)-1
        write(13,*) coA(i)
      ENDDO
      DO i=1,ndim
        write( 9,*) coB(i)
        write(14,*) rl(i)
      ENDDO
      close( 9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
*     STOP
*
      END
*******************************************************
      subroutine writevector(ndim,x)
      implicit none
      integer ndim, i
      real x(1)
      open(10,file='A.sol', status='UNKNOWN')
*      print*, 'Vector'
      do i = 1, ndim
         write(10,*) x(i)
*         print '(A,I4,A,G10.10)', '[', i, '] = ', x(i)
      enddo
      close(10)
*      stop
      end
********************************************************
      subroutine writevectori(ndim,x)
      implicit none
      integer ndim, i
      integer x(1)
*      open(10,file='A.sol', status='UNKNOWN')
      print*, 'Vector'
      do i = 1, ndim
*         write(10,*) x(i)
         print '(A,I4,A,I8)', '[', i, '] = ', x(i)
      enddo
      close(10)
      end
********************************************************
      subroutine printmat( a, n)
      integer n, i, j
      real a(n,n)
      write(*,*) ' '
      do i = 1, n
         write(*,'(10(f10.4,x))') (a(i,j), j = 1, n)
      enddo
      write(*,*) ' '
      end
