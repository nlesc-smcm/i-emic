!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_wacsr

CONTAINS

SUBROUTINE wacsr (n, MaxNnz, x)

USE m_build
USE m_dump

INTEGER				, INTENT(IN)	:: n
INTEGER				, INTENT(IN)    :: MaxNnz
TYPE (csrmatrix)		, POINTER 	:: x

!     Allocate a matrix descriptor and the segments for a sparse matrix
!     in CSR format.
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at 'x'.

!     Arguments:
!     ==========
!     N        i   Number of rows and columns in matrix
!     MaxNnz      i   Number of non-zeros in matrix
!     x   o   Location of Matrix descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'wacsr'
INTEGER                      :: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam

!     Check input sensible:
IF (n < 0)   CALL dump(__FILE__,__LINE__,'Illegal n < 0')
IF (MaxNnz < 0) CALL dump(__FILE__,__LINE__,'Illegal MaxNnz < 0')

PRINT '(A, I0, A, I0, A)', 'Allocating a CSR-matrix with ', n, ' rows and ', MaxNnz, ' non-zeros'
#endif

!     Request segment for the matrix descriptor:
ALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Request segments for a sparse CSR matrix:

ALLOCATE( x%beg(1:n+1), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( x%jco(1:MaxNnz), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( x%co(1:MaxNnz), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Save information in descriptor:
!     Storage Type, Matrix Dimension and the locations of the vectors
!     for the CSR representation of a matrix.

x%typ 	= csrtp
x%n	= n
x%nnz	= MaxNnz
x%beg	= 1

END SUBROUTINE wacsr

END MODULE
