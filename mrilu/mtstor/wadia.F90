!#begindoc

#ifndef WITH_UNION

#define diamatrix  anymatrix

#endif

MODULE m_wadia

CONTAINS

SUBROUTINE wadia (n, blksiz, x)

USE m_dump
USE m_build

INTEGER			, INTENT(IN)                     :: n
INTEGER			, INTENT(IN)                     :: blksiz
TYPE (diamatrix)	, POINTER 	  		 :: x

!     Allocate a matrix descriptor and the segment for a sparse
!     (block-) DIAgonal matrix, DIAtp.
!     The descriptor is filled so that the requested segment can be
!     referenced through the descriptor at 'x'

!     Arguments:
!     ==========
!     N        	i   Number of rows and columns in matrix.
!     BlkSiz   	i   Number of rows/colums in a block  ('BlkSiz' >= 1).
!                   'N' should be an integer multiple of 'BlkSiz'!
!     x   	o   Location of Matrix descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER 	:: rounam = 'wadia'
INTEGER				:: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam

!     Check input sensible:
IF (n < 0)        CALL dump(__FILE__,__LINE__,'Illegal n < 0')
IF (blksiz < 0)   CALL dump(__FILE__,__LINE__,'Illegal blksiz < 0')

PRINT '(A, I0, A, I0, A, I0)', 'Allocating ', n, 'x', n, ' DIA-matrix with blocksize', blksiz
#endif

IF (MOD(n, blksiz) /= 0)   CALL dump(__FILE__,__LINE__,'Internal error: blksiz not compatible with number rows n')

!     Request segment for the matrix descriptor:

ALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Request segment for a diagonal matrix:

ALLOCATE( x%com(1:blksiz,1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Save information in descriptor:
!     Storage Type, Matrix Dimension, the number of blocks, the block
!     size and the locations of a (block-) Diagonal matrix.

x%typ   	= diatp
x%n   		= n
x%blksiz 	= blksiz

!     End of  wadia
END SUBROUTINE wadia

END MODULE







