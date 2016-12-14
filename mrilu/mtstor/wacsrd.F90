!#begindoc

#ifndef WITH_UNION

#define csrdmatrix  anymatrix

#endif

MODULE m_wacsrd

CONTAINS
 
SUBROUTINE wacsrd (n, blksiz, MaxNnz, x)

USE m_dump
USE m_build
USE m_wacsr
USE m_wadia

INTEGER, INTENT(IN)                     :: n
INTEGER, INTENT(IN)                     :: blksiz
INTEGER, INTENT(IN)                     :: MaxNnz
TYPE (csrdmatrix), POINTER              :: x

!     Allocate a matrix descriptor and the segments for a CSRD,
!     Compact Sparse Row Diagonal extracted, matrix.
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at  'x'.

!     Arguments:
!     ==========
!     N        i   Number of rows/columns in matrix
!     BlkSiz   i   Number of rows/columns in diagonal blocks.
!                  'N' should be an integer multiple of 'BlkSiz'!
!     MaxNnz      i   Number of non-zeros in the off-diagonal matrix.
!     x      o   Location of Matrix descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER 	:: rounam = 'wacsrd'
INTEGER				:: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam

!     Check input sensible:
IF (n < 0)   CALL dump(__FILE__,__LINE__,'Illegal n < 0')
IF (MaxNnz < 0) CALL dump(__FILE__,__LINE__,'Illegal MaxNnz < 0')

PRINT '(A, I0, A, I0, A, I0)', 'Allocating ', n, 'x', n, ' CSRD-matrix with ', MaxNnz, ' non-zeros and blocksize ', blksiz
#endif

!     Request segment for the storage descriptor
ALLOCATE ( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Request segment for the vector with the indices of the last
!     elements in the Lower Triangular part of the matrix:
ALLOCATE( x%lotr(1:n), STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate storage for a sparse (block-) diagonal matrix (DIAtp)
!     Request segment for a diagonal matrix:

CALL wadia(n, blksiz, x%dia)

!     Allocate storage for the off-diagonal matrix (type: CSR)
!     Request segments for a sparse CSR matrix:

CALL wacsr(n, MaxNnz, x%offd)

!     Save information in storage descriptor:
!     Storage Type, Matrix Dimension and the matrices with the diagonal
!     and off-diagonal  values.

x%typ  = csrdtp
x%n  = n

END SUBROUTINE wacsrd

END MODULE










