!#begindoc

#ifndef WITH_UNION

#define scbmmatrix  anymatrix

#endif

MODULE m_possred

CONTAINS
 
SUBROUTINE possred (A, x2, b, x)

USE m_build
USE m_diavec
USE m_csrvec

TYPE (scbmmatrix)			, POINTER		:: A
DOUBLE PRECISION, DIMENSION(1:A%nschur)	, INTENT(IN)		:: x2
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN OUT)	:: b
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(OUT)           :: x

!     Postprocessing after Solving the Reduced system.

!     Compute segment  x_1  of  x = ( x_1 | x_2 )', with:
!        x_2 := x2
!        x_1 := (1/A_11) (b_1 - A_12 x_2)
!     and store  ( x_1 | x_2 )', the solution of  A x = b, into  x.

!     The matrix  A  of the original system is partitioned as:
!            ( A_11 | A_12 )
!        A = (------+------)
!            ( A_21 | A_22 )
!     where  A_11  is a block-diagonal sub-matrix.

!     Arguments:
!     ==========
!     A%N	i   Number of rows/columns in matrix  A.
!     A%G	i   Number of rows/columns in left upper partition A_11.
!     A%a11d    i   Location of descriptor for the inverse of
!                  left upper partition A_11 of matrix A, inv(A11).
!     A%a21    i   Location of descriptor for the right upper
!                  partition A_12 of matrix A.
!     x2       i   Solution of reduced equation, x_2.
!     b        io  In:  Right-hand side of the complete equation
!                       A x = b.
!                  Out: --- undefined.
!     x        o   Solution of the complete equation  A x = b.

!#enddoc

!     Local Variables
!     ===============

CHARACTER (LEN=*), PARAMETER :: rounam = 'possred'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Shift x_2 into position in 'x':

x(A%g+1:A%n) = x2

!     {  x(A%G+1:A%N) = x_2  }

!     Compute b(1:A%G) = b(1:A%G) - A(1:A%G,A%G+1:A%N) x(A%G+1:A%N) = b_1 - A_12 x_2:

CALL csrvec (-1.0D0, A%a12, x(A%g+1:A%n), b(1:A%g))

!     Compute  x(1:A%G) = (1/A_11) (b_1 - A_12 x_2)

CALL diavec (.false., 1.0D0, A%a11d, b, x(1:A%g))


END SUBROUTINE possred

END MODULE

