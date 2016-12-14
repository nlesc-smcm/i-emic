!#begindoc
 
#ifndef WITH_UNION

#define prcmatrix  anymatrix

#endif

MODULE m_waprc

CONTAINS

SUBROUTINE waprc (n, x)

USE m_dump
USE m_build

INTEGER			, INTENT(IN)		:: n
TYPE (prcmatrix)	, POINTER               :: x

!     Allocates a descriptor for the complete preconditioner of an  NxN
!     matrix, and allocates the segments for:
!     . the scale factors for the matrix and right hand side  ScaFct
!     . the permutation vector for the Red-Black reordering   PerRB
!     . the column indices of the last non-zero in left part  LasLeft
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at  'x'.

!     Arguments:
!     ==========
!     N    i   Number of rows/columns in the matrix.
!     x    o   Location of Matrix descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER 	:: rounam = 'waprc'
INTEGER				:: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam

IF (n < 0)   CALL dump(__FILE__,__LINE__,'Illegal n < 0')

PRINT '(A, I0, A, I0, A)', 'Allocating ', n, 'x', n, ' PRC-matrix'
#endif

!     Request segment for the storage descriptor:
ALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate storage for the scale factors  ScaFct:
ALLOCATE( x%scale(1:n), STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate storage for the permutation vector  PerRB:
ALLOCATE( x%perrb(1:n), STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')


!     Save information in storage descriptor:
!     Locations of the unknown descriptors are initialised with NULLINX.

x%typ   = prctp
x%n  	= n
x%g  	= 0
x%nschur= n
NULLIFY(x%aro)
NULLIFY(x%mlp)

!     End of  waprc
END SUBROUTINE waprc

END MODULE
