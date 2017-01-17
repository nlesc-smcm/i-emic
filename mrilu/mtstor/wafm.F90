!#begindoc

#ifndef WITH_UNION

#define fmmatrix  anymatrix

#endif

MODULE m_wafm

CONTAINS
 
SUBROUTINE wafm (n, x)

USE m_dump
USE m_build

INTEGER			, INTENT(IN)                  :: n
TYPE (fmmatrix)		, POINTER                     :: x

!     Allocate a matrix descriptor and the segment for a Full Matrix,
!     FMtp, of size  N x N.
!     The descriptor is filled so that the requested segment can be
!     referenced through the descriptor at 'x'

!     Arguments:
!     ==========
!     N        i   Number of rows and columns in matrix.
!     x   o   Location of Matrix descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER 	:: rounam = 'wafm'
INTEGER				:: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam

!     Check input sensible:
IF (n < 0)   CALL dump(__FILE__,__LINE__,'Illegal n < 0')

PRINT '(A, I0, A, I0, A)', 'Allocating ', n, 'x', n, ' FN-matrix'
#endif

!     Request segment for the matrix descriptor:

ALLOCATE( x, STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Request segment for a full matrix:

ALLOCATE( x%com(1:n,1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Save information in descriptor:
!     Storage Type, Matrix Dimension and the locations of an NxN matrix.

x%typ 	= fmtp
x%n 	= n

!     End of  wafm
END SUBROUTINE wafm

END MODULE
