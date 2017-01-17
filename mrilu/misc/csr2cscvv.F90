!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix
#define partmatrix  anymatrix
#define csrdmatrix  anymatrix

#endif

MODULE m_csr2cscvv

CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE csr2cscvv (a, l)

USE m_build
USE m_wacsr

TYPE (csrmatrix)			, POINTER		:: a
TYPE (csrmatrix)			, POINTER		:: l

!     Transpose CSR to CSC format and vice versa


!     Store the sub-matrix A(1:nr,1:nc) in the CSR matrix  A,
!     into the existing CSC sub-matrix  L, the Lower triangular block,
!     in a partition.
!     The CSR matrix  A  is stored in
!        [Nupp, a%beg, a%jco, a%co],
!     with row and column numbers in  1:Nlast.
!     The CSC matrix  L, = A21,  is stored into
!        [Nupp, l%beg, l%jco, l%co],
!     with the original row/column numbers in  1:NeqDon+Nlast.

!     Arguments:
!     ==========
!     NeqDon   i   Number of equations (rows/columns) done.
!                  Nlast+NeqDon is the size of the root matrix.
!     Nlast    i   Number of rows and columns of the matrix, containing
!                  the sub-matrix  A.
!     Nupp     i   Number of columns in sub-matrix  A.  (Nupp < Nlast)
!     Permut   i   Permut(NeqDon+1:NeqDon+Nlast): mapping from
!                  1:Nlast  into  1:NeqDon+Nlast, the original
!                  row/column numbers.

!#enddoc

!     Local Variables:
!     ================
INTEGER :: i, j, ij, next, nr, nc

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'csr2cscvv'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Compute the length of the columns of A:

call wacsr(a%n,a%nnz,l)
nr=a%n
nc=nr
l%beg(1:nc+1) = 0

DO   i = 1, nr
  DO   ij = a%beg(i), a%beg(i+1) - 1
    j       = a%jco(ij) + 1
    l%beg(j) = l%beg(j) + 1
  END DO
END DO

!     {  A(j: 1<=j<=Nupp: l%beg(j+1) = "length j-th column in A")

!     Compute indices to beginning of column segments from the length

l%beg(1) = 1
DO   j = 1, nc
  l%beg(j+1) = l%beg(j+1) + l%beg(j)
END DO

!     Store the elements and adjust the column numbers:

!     {  A(j: 1<=j<=Nupp: l%beg(j) = "position next element in column j")  }

DO   i = 1, nr
  DO   ij = a%beg(i), a%beg(i+1) - 1
    j          = a%jco(ij)
    next       = l%beg(j)
    l%co(next)  = a%co(ij)
    l%jco(next) = i
    l%beg(j)    = next + 1
  END DO
END DO

!     Reshift l%beg:

DO   j = nc, 1, -1
  l%beg(j+1) = l%beg(j)
END DO
l%beg(1) = 1

!     End of  csr2cscvv
END SUBROUTINE csr2cscvv


END MODULE