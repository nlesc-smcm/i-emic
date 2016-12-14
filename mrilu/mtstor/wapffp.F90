!#begindoc
 
#ifndef WITH_UNION

#define partmatrix  anymatrix

#endif

MODULE m_wapffp

CONTAINS

SUBROUTINE wapffp (n, offset, x)

USE m_dump
USE m_build
USE m_wafm

INTEGER			, INTENT(IN)                  :: n
INTEGER			, INTENT(IN)                  :: offset
TYPE (partmatrix)	, POINTER                     :: x

!     Allocate Partition with Full Factors and Pivots.

!     Allocate a partition
!     descriptor and the segments for the last, a PFFP type, partition
!     in a MLP type matrix.
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at  'x'.

!     Arguments:
!     ==========
!     N        	i   Number of unknowns and equations represented by the partition.
!     Offset   	i   Offset of partition in MLP matrix  P.
!     x   	o   Location of the Partition descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER 	:: rounam = 'wapffp'
INTEGER				:: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam

PRINT '(A, I0, A, I0, A, I0)', 'Allocating ', n, 'x', n, ' PFFP-matrix with offset ', offset
#endif

!     Request workspace for partition descriptor:
ALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate workspace for pivot elements:
ALLOCATE( x%piv(1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate storage for a Full Matrix (type: FMtp):
CALL wafm (n, x%fm)

!     Save information in partition descriptor:
!     Partition Type, Offset for block partition in root matrix and the
!     locations of the diagonal and off-diagonal elements:

x%typ   = pffptp
x%off  = offset
x%n = n  
NULLIFY(x%prev)
NULLIFY(x%next)


END SUBROUTINE wapffp

END MODULE




