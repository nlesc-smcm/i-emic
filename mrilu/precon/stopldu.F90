!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix
#define partmatrix  anymatrix
#define csrdmatrix  anymatrix

#endif

MODULE m_stopldu

CONTAINS

SUBROUTINE stopldu (neqdon, nlast, nupp, a, permut, xpart)

USE m_build
USE m_wapldu
USE m_dump

INTEGER						, INTENT(IN)		:: neqdon
INTEGER						, INTENT(IN)            :: nlast
INTEGER						, INTENT(IN)            :: nupp
TYPE (csrdmatrix)				, POINTER		:: a
INTEGER, DIMENSION(1:neqdon+nlast)		, INTENT(IN)            :: permut
TYPE (partmatrix)				, POINTER		:: xpart

!     Store the left and upper border of the (sub-)matrix A into
!     a new, PLDU type, partition.

!     The matrix  A  has been partitioned into:
!        A11 | A12
!        ----+----
!        A21 | A22
!     by means of the permutation given in Permut.

!     A11 is a block diagonal matrix stored in the structure:
!        [Nupp, a%offd%dia%com]

!     The off-diagonal part of matrix  A  is stored, in CSR format, in
!     the structure:
!        [Nlast, a%offd%beg, a%offd%jco, a%offd%co]

!     Block A11 is a diagonal matrix.
!     Block A12 is Upper triangular.
!     Block A21 is Lower triangular.

!     Arguments:
!     ==========
!     NeqDon   i   Number of equations (rows/columns) done.
!     Nlast    i   The size of the partitioned matrix A.
!     Nupp     i   The size of the upper partition A11 in matrix A.
!     A%dia%blksiz   i   Number of rows/columns in the diagonal blocks of A.
!                  Nlast+NeqDon is the size of the root matrix.
!                  New row and column numbers are incremented with
!                  NeqDon.
!     a%offd%beg     i   a%offd%beg(r), 1<=r<=Nlast: Index in 'a%offd%jco' and 'a%offd%co' of
!                  the first non-zero off-diagonal element in row  r
!                  of matrix  A.
!                  a%offd%beg(Nlast+1): Index in 'a%offd%jco' and 'a%offd%co of last
!                  nonzero element + 1.
!     a%offd%jco     i   a%offd%jco(nz), a%offd%beg(r)<=nz<a%offd%beg(r+1), 1<=r<=Nlast: Column
!                  number of the non-zero off-diagonal element  a%offd%co(nz)
!                  in row  r  of matrix  A.
!     a%offd%co      i   a%offd%co(a%offd%beg(r):a%offd%beg(r+1)-1), 1<=r<=Nlast: Values of the
!                  nonzero off-diagonal elements in row  r  of matrix  A.
!     a%dia%com     i   Inverse of main (block-)diagonal of left upper
!                  submatrix, A_11, of matrix A.
!                  Each block is stored in column major order.
!     A%lotr   i   A%lotr(i): Index, in a%offd%jco and a%offd%co, of last nonzero
!                  element in Lower Triangular block of A?
!     Permut   i   Vector defines the permutation which partitions the
!                  matrix into the 4 blocks A11, A12, A21 and A22.
!     xpart   i   Location of descriptor of the new partition
!                  with the L, D and U factors.

!#enddoc

!     Local Variables:
!     ================

INTEGER 					:: nnza12, nnza21

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'stopldu'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Number of nonzeros in A12:
nnza12 = SUM(a%offd%beg(2:nupp+1) - a%offd%beg(1:nupp))

!     Number of nonzeros in A21:
nnza21 = SUM( a%lotr(nupp+1:nlast) ) - SUM( a%offd%beg(nupp+1:nlast) )

nnza21 = nnza21 + (nlast-nupp)+1

!     Allocate a PLDU partition [xpart].
!     (Descriptor has been initialised!)
CALL wapldu (nupp, a%dia%blksiz, nnza12, nnza21, neqdon, xpart)

!     Store the Block diagonal matrix [A11] into the Block Diagonal,
!     stored by Column, part of partition [xpart]:

#ifdef DEBUG
  IF (nupp /= xpart%dia%n  ) THEN
    PRINT '(A, 2X, A, /, 3X, A)', 'Internal error in', rounam, 'Inconsistent data in PLDU or DIAtp representation!'
    CALL dump(__FILE__,__LINE__,'Inconsistent data')
  END IF
#endif

!     Store the block diagonal [A11] into the Diagonal part of
!     partition [xpart]:

xpart%dia%com = a%dia%com(1:a%dia%blksiz,1:nupp)

!     Store the the Upper triangular block, a CSR matrix, [A12] into the
!     CSR sub-matrix of partition [xpart]:

CALL stocsrp (neqdon, nlast, nupp, permut, a%offd,  xpart%utr)

!     Store the the Lower triangular block, a CSC matrix, [A21] into the
!     CSC sub-matrix of partition [xpart]:

CALL stocscp (neqdon, nlast, nupp, permut, a%offd, xpart%ltr)

END SUBROUTINE stopldu

!=======================================================================

!#begindoc

!-----------------------------------------------------------------------

SUBROUTINE stocscp (neqdon, nlast, nupp, permut, a, l)

USE m_build

INTEGER					, INTENT(IN)		:: neqdon
INTEGER					, INTENT(IN)    	:: nlast
INTEGER					, INTENT(IN)            :: nupp
INTEGER, DIMENSION(1:neqdon+nlast)	, INTENT(IN)          	:: permut
TYPE (csrmatrix)			, POINTER		:: a
TYPE (cscmatrix)			, POINTER		:: l

!     Store the Lower triangular block in CSC format:


!     Store the sub-matrix A(Nupp+1:Nlast,1:Nupp) in the CSR matrix  A,
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
INTEGER :: i, j, ij, next

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'stocscp'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Compute the length of the columns of A:

l%beg(1:nupp+1) = 0

DO   i = nupp+1, nlast
  DO   ij = a%beg(i), a%beg(i+1) - 1
    j       = a%jco(ij) + 1
    IF (j <= nupp) THEN
      l%beg(j) = l%beg(j) + 1
    END IF
  END DO
END DO

!     {  A(j: 1<=j<=Nupp: l%beg(j+1) = "length j-th column in A")

!     Compute indices to beginning of column segments from the legnth

l%beg(1) = 1
DO   j = 1, nupp
  l%beg(j+1) = l%beg(j+1) + l%beg(j)
END DO

!     Store the elements and adjust the column numbers:

!     {  A(j: 1<=j<=Nupp: l%beg(j) = "position next element in column j")  }

DO   i = nupp+1, nlast
  DO   ij = a%beg(i), a%beg(i+1) - 1
    j          = a%jco(ij)
    IF (j <= nupp) THEN
      next       = l%beg(j)
      l%co(next)  = a%co(ij)
      l%jco(next) = permut(neqdon+i)
      l%beg(j)    = next + 1
    END IF
  END DO
END DO

!     Reshift l%beg:

DO   j = nupp, 1, -1
  l%beg(j+1) = l%beg(j)
END DO
l%beg(1) = 1

!     End of  stocscp
END SUBROUTINE stocscp

!=======================================================================

!#begindoc

!-----------------------------------------------------------------------

SUBROUTINE stocsrp (neqdon, nlast, nupp, permut, a, u)

USE m_build

INTEGER					, INTENT(IN)		:: neqdon
INTEGER					, INTENT(IN)            :: nlast
INTEGER					, INTENT(IN)            :: nupp
INTEGER	, DIMENSION(1:neqdon+nlast)	, INTENT(IN)            :: permut
TYPE (csrmatrix)			, POINTER		:: a
TYPE (csrmatrix)			, POINTER		:: u

!     Store the CSR sub-matrix  A,  stored in
!        [Nupp, a%beg, a%jco, a%co],
!     with column numbers in  1:Nlast,
!     into the existing CSR sub-matrix  U, stored in
!        [Nupp, u%beg, u%jco, u%co],
!     with the original row/column numbers in  1:NeqDon+Nlast.

!     Arguments:
!     ==========
!     NeqDon   i   Number of equations (rows/columns) done.
!                  Nlast+NeqDon is the size of the root matrix.
!     Nlast    i   Number of rows and columns of the matrix, containing
!                  the sub-matrix  A.
!     Nupp     i   Number of rows in matrix  A.  (Nupp < Nlast)
!     Permut   i   Permut(NeqDon+1:NeqDon+Nlast): mapping from
!                  1:Nlast  into  1:NeqDon+Nlast, the original
!                  row/column numbers.

!#enddoc

!     Local Variables:
INTEGER :: nz, rnr, unnz

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'stocsrp'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     u%beg(.) index of next nonzero element in u%jco and u%co:

unnz     = 0
u%beg(1) = 1

DO rnr = 1, nupp
  DO nz = a%beg(rnr), a%beg(rnr+1) - 1
    unnz      = unnz+1
    u%co(unnz) = a%co(nz)
!   u%jco(.) original column number of nonzero element.
    u%jco(unnz) = permut(neqdon+a%jco(nz))
  END DO
  u%beg(rnr+1) = unnz+1
END DO


END SUBROUTINE stocsrp

END MODULE