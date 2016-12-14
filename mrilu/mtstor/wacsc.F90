!#begindoc
 
#ifndef WITH_UNION

#define cscmatrix  anymatrix

#endif

MODULE m_wacsc

CONTAINS

SUBROUTINE wacsc (n, MaxNnz, x)

USE m_dump
USE m_build

INTEGER			, INTENT(IN)                :: n
INTEGER			, INTENT(IN)                :: MaxNnz
TYPE (cscmatrix)	, POINTER                   :: x

!     Allocate a matrix descriptor and the segments for a sparse matrix
!     in CSC format.
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at 'x'.

!     Arguments:
!     ==========
!     N        i   Number of rows and columns in matrix
!     MaxNnz      i   Number of non-zeros in matrix
!     x   o   Location of Matrix descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'wacsc'
INTEGER                      :: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam

!     Check input sensible:
IF (n < 0)   CALL dump(__FILE__,__LINE__,'Illegal n < 0')
IF (MaxNnz < 0) CALL dump(__FILE__,__LINE__,'Illegal MaxNnz < 0')

PRINT '(A, I0, A, I0, A, I0, A)', 'Allocating ', n, 'x', n, ' CSS-matrix with ', MaxNnz, ' non-zeros'
#endif

!     Request segment for the matrix descriptor:
ALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Request segments for a sparse CSC matrix:

ALLOCATE( x%beg(1:n+1), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( x%jco(1:MaxNnz), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( x%co(1:MaxNnz), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Save information in descriptor:
!     Storage Type, Matrix Dimension and the locations of the vectors
!     for the CSC representation of a matrix.

x%typ 	= csctp
x%n	= n
x%nnz	= MaxNnz
x%beg	= 1

END SUBROUTINE wacsc

END MODULE













