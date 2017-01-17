!#begindoc

#ifndef WITH_UNION

#define mlpmatrix  anymatrix

#endif

MODULE m_wamlp

CONTAINS
 
SUBROUTINE wamlp (n, x)

USE m_dump
USE m_build

INTEGER			, INTENT(IN)                  :: n
TYPE (mlpmatrix)	, POINTER                     :: x

!     Allocate:
!     . The matrix descriptor for a Multilevel Partitioned LDU
!       matrix.  The list of partitions is empty.
!     . The segment for the permutation vector.
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at  X.

!     Arguments:
!     ==========
!     N   	i   Number of unknowns and equations in matrix
!     x   	o   Location of Matrix descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'wamlp'
INTEGER			     :: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam

!     Check input sensible
IF (n < 0)   CALL dump(__FILE__,__LINE__,'Illegal n < 0')

PRINT '(A, I0, A, I0, A)', 'Allocating ', n, 'x', n, ' MLP-matrix'
#endif

!     Request segment for the matrix descriptor:
ALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Request segment for the permutation vector:
ALLOCATE( x%perm(1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Save information in storage descriptor:
!     Storage Type, Matrix Dimension and the location of the
!     Preconditioner permutation.

x%typ  	= mlptp
x%n  	= n

!     Initialise locations of first and last partition to the head of
!     the doubly linked list:

NULLIFY(x%first)
NULLIFY(x%last)

!     End of wamlp
END SUBROUTINE wamlp

END MODULE
