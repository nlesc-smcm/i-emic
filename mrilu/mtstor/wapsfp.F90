!#begindoc
 
#ifndef WITH_UNION

#define partmatrix  anymatrix

#endif

MODULE m_wapsfp

CONTAINS

SUBROUTINE wapsfp (n, MaxNnz, offset, x)

USE m_dump
USE m_build
USE m_wacsr
USE m_wadia

INTEGER			, INTENT(IN)		:: n
INTEGER			, INTENT(IN)     	:: MaxNnz
INTEGER			, INTENT(IN)            :: offset
TYPE (partmatrix)	, POINTER               :: x

!     Allocate Partition with Sparse Factors and Pivots.

!     Allocate a partition
!     descriptor and the segments for the last, a PSFP type, partition
!     in a MLP type matrix.
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at 'x'.

!     Arguments:
!     ==========
!     N        	i   Number of unknowns and equations represented by the
!                   partition.
!     MaxNnz      	i   Number of non-zeros in the off-diagonal matrix.
!     Offset   	i   Offset of partition in MLP matrix  P.
!     x   	o   Location of the Partition descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER 	:: rounam = 'wapsfp'
INTEGER				:: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam

PRINT '(A, I0, A, I0, A, I0)', 'Allocating ', n, 'x', n, ' PSFP-matrix with ', MaxNnz, ' non-zeros '
#endif

!     Request workspace for partition descriptor:
ALLOCATE( x, STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate workspace for pivot elements.
ALLOCATE( x%piv(1:n), STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate workspace for last lower triangular elements.
ALLOCATE( x%lnzl(1:n), STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate storage for a sparse diagonal matrix (type: DIAtp)
CALL wadia (n, 1, x%dia)

!     Allocate a storage for a CSR matrix (type: CSR)
CALL wacsr (n, MaxNnz, x%offd)

!     Save information in partition descriptor:
!     Partition Type, Offset for block partition in root matrix and the
!     locations of the diagonal and off-diagonal elements:

x%n    = n
x%typ  = psfptp
x%off  = offset
NULLIFY(x%prev)
NULLIFY(x%next)

END SUBROUTINE wapsfp

END MODULE
