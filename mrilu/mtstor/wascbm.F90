!#begindoc
 
#ifndef WITH_UNION

#define scbmmatrix  anymatrix

#endif

MODULE m_wascbm

CONTAINS

SUBROUTINE wascbm (n, nupp, blksiz, nz12, nz21, nz22, x)

USE m_dump
USE m_build
USE m_wacsr
USE m_wadia

INTEGER			, INTENT(IN)                    :: n
INTEGER			, INTENT(IN)                    :: nupp
INTEGER			, INTENT(IN)                    :: blksiz
INTEGER			, INTENT(IN)                    :: nz12
INTEGER			, INTENT(IN)                    :: nz21
INTEGER			, INTENT(IN)                    :: nz22
TYPE (scbmmatrix)	, POINTER                     	:: x

!     Allocate a matrix descriptor and the segments for a Schur
!     Complement from Block Matrix.
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at 'x'.

!     Arguments:
!     ==========
!     N        	i   Number of rows/columns in the matrix, A.
!     Nupp     	i   Number of rows/columns in the left upper part, A11.
!                   The following relation should hold:  0 < Nupp < N
!     BlkSiz   	i   Number of rows/columns in diagonal blocks.
!                   'N' and 'Nupp' should be an integer multiple of
!                   'BlkSiz'!
!     nz12     	i   Number of non-zeros in the right upper block, A12, of
!                   the matrix.
!     nz21     	i   Number of non-zeros in the lower left block, A21, of
!                   the matrix.
!     nz22     	i   Number of non-zeros in the lower right block, A22, of
!                  the matrix.
!     x      	o   Location of Matrix descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER 	:: rounam = 'wascbm'
INTEGER				:: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam

!     Check if input sensible:
IF (n < 0)       CALL dump(__FILE__,__LINE__,'Illegal n < 0')
IF (nupp < 0)    CALL dump(__FILE__,__LINE__,'Illegal nupp < 0')
IF (blksiz <= 0) CALL dump(__FILE__,__LINE__,'Illegal blksiz <= 0')
IF (nz12 < 0)    CALL dump(__FILE__,__LINE__,'Illegal nz12 < 0')
IF (nz21 < 0)    CALL dump(__FILE__,__LINE__,'Illegal nz21 < 0')
IF (nz22 < 0)    CALL dump(__FILE__,__LINE__,'Illegal nz22 < 0')

PRINT '(A, I0, A, I0, A, I0, A, I0)', 'Allocating ', n, 'x', n, &
 ' SCBM-matrix with ', nupp, &
 ' rows/columns in left upper part and blocksize', blksiz
#endif

!     Request segment for the storage descriptor
ALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate storage for a sparse (block-) diagonal matrix (DIAtp)
CALL wadia (nupp, blksiz, x%a11d)

!     Allocate storage for the left upper submatrix (type: CSR)
CALL wacsr (nupp, nz12, x%a12)

!     Allocate storage for the lower left submatrix (type: CSR)
CALL wacsr (n-nupp, nz21, x%a21)

!     Allocate storage for the lower right submatrix (type: CSR)
CALL wacsr (n-nupp, nz22, x%a22)

!     Save information in storage descriptor:
!     Storage Type, Matrix Dimensions and the descriptors of the
!     submatrices.

x%typ  	= scbmtp
x%n  	= n
x%g 	= nupp
x%nschur= n-nupp

!     End of  wascbm
END SUBROUTINE wascbm

END MODULE
