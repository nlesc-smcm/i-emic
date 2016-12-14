!#begindoc
 
#ifndef WITH_UNION

#define prcmatrix  anymatrix

#endif

MODULE m_matvecp

CONTAINS

SUBROUTINE matvecp (s, x, Sx)

USE m_dump
USE m_build
USE m_matvec
USE m_dperv

TYPE (prcmatrix)			, POINTER		:: S
DOUBLE PRECISION, DIMENSION(1:S%nschur)	, INTENT(IN)         	:: x
DOUBLE PRECISION, DIMENSION(1:S%nschur)	, INTENT(OUT)         	:: Sx
 
!     Matrix vector product with permutation.
!     Computes the matrix vector product
!        Prc S Prc' x
!     and stores the result into  'Sx'.

!     Arguments:
!     ==========
!     S%NSCHUR		i   Number of rows/columns in matrix  S.
!     S%aro     	i   Location of descriptor of the
!                           Schur-complement matrix  S.
!     S%mlp%perm   	i   Permutation vector which characterizes the
!                           preconditioner permutation matrix  Prc.
!     x        		i   Input vector
!     Sx       		o   Matrix vector product:  Prc S Prc' x

!#enddoc

!     Local Parameters:
!     =================


CHARACTER (LEN=*), PARAMETER :: rounam = 'matvecp'

!     Local variables:
!     ================

INTEGER						:: ier
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: temp

ALLOCATE( temp(1:S%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Calculate:  y(1:S%NSCHUR) := Prc' x (= inv(Prc) x):
CALL dperv (.true., S%mlp%perm, x, temp)


!     Calculate:  Sx(1:S%NSCHUR) := S y (= S Prc' x):
CALL matvec (S%nschur, 1.0D0, S%aro, temp, Sx)

!     Permute solution into place
!     Calculate:  Sx := Prc Sx
!                    = Prc S Prc' x
CALL dperv (.false., S%mlp%perm, Sx, Sx )

!     End of  matvecp

DEALLOCATE( temp, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

END SUBROUTINE matvecp

END MODULE















