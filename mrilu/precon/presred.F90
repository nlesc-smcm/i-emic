!#begindoc
 
#ifndef WITH_UNION

#define scbmmatrix  anymatrix

#endif

MODULE m_presred

CONTAINS

SUBROUTINE presred (A, b, x, b2, x2)

USE m_build
USE m_diavec
USE m_csrvec
USE m_dump

TYPE (scbmmatrix)			, POINTER	:: A
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)    :: b
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)    :: x
DOUBLE PRECISION, DIMENSION(1:A%nschur)	, INTENT(OUT)   :: b2
DOUBLE PRECISION, DIMENSION(1:A%nschur)	, INTENT(OUT)   :: x2

!     Preprocessing before Solving the Reduced system.

!     Construct both right-hand side 'b2' and starting value for the
!     solution 'x2' of the reduced system
!        S x2  =  b2
!     from the right-hand side 'b' and starting value for the solution
!     'x' of the original system
!        A x = b.

!     with
!        S    = A_22 - A_21 (1/A_11) A_12
!               is the Schur complement of  A_11  in  A.
!        b2 = b_2 - A_21 (1/A_11) b_1
!     and
!        A_11  is the  GxG  left upper block-diagonal sub-matrix of  A.

!     The original equation
!        A x = b
!     has beeen partitioned as:

!        ( A_11 | A_12 ) ( x_1 )     ( b_1 )
!        (------+------) (-----)  =  (-----)
!        ( A_21 | A_22 ) ( x_2 )     ( b_2 )

!     == {  equivalent with:  }

!        A_11 x_1  =  b_1 - A_12 x_2      ,
!        A_22 x_2  =  b_2 - A_21 x_1

!     ==

!        x_1       =  (1/A_11) (b_1 - A_12 x_2)  ,
!        A_22 x_2  =  b_2 - A_21 (1/A_11) b_1 + A_21 (1/A_11) A_12 x_2

!     ==

!        x_1                              =  (1/A_11) (b_1 - A_12 x_2)    ,
!        (A_22 - A_21 (1/A_11) A_12) x_2  =  (b_2 - A_21 (1/A_11) b_1)

!     The reduced problem
!        S x_2 =  b_2 - A_21 (1/A_11) b_1
!     or
!        S x2  = b2
!     can be solved using a CG-like method with a preconditioner.

!     Arguments:
!     ==========
!     A%N	i   Number of rows/columns in matrix  A.
!     A%G	i   Number of rows/columns in left upper partition A_11.
!     A%a11d   	i   Location of descriptor for the inverse of
!                   left upper partition A_11 of matrix A, inv(A11).
!     A%a21    	i   Location of descriptor for the lower left partition A_21 of matrix A.
!     b        	i   Right-hand side of equation  (b_1 | b_2)'.
!     x        	i   Starting value of solution  (x_1 | x_2)'.
!     b2     	o   Right-hand side of the reduced problem.
!     x2       	o   Starting value of solution of the reduced system.

!#enddoc

!     Local Variables
!     ===============

INTEGER						:: ier
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: wrk

CHARACTER (LEN=*), PARAMETER :: rounam = 'presred'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Construct the starting vector for the reduced problem:

x2 = x(A%g+1:A%n)

!     Construct Right-Hand side for reduced problem, using wrk(1:G) as
!     scratch space.

!     Construct:  wrk(1:A%G) = (1/A_11) * b_1

ALLOCATE( wrk(1:A%g), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

CALL diavec (.false., 1.0D0, A%a11d, b, wrk)

!     Construct:  b2 = b_2 - A_21*(1/A_11)*b_1

b2 = b(A%g+1:A%n)

CALL csrvec (-1.0D0, A%a21, wrk, b2)

DEALLOCATE( wrk, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

END SUBROUTINE presred

END MODULE