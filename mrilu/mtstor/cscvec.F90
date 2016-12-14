!#begindoc
 
#ifndef WITH_UNION

#define cscmatrix  anymatrix

#endif

MODULE m_cscvec

CONTAINS

SUBROUTINE cscvec (alpha, a, x, y )

USE m_build

DOUBLE PRECISION			, INTENT(IN)            :: alpha
TYPE (cscmatrix)			, POINTER	        :: a
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)            :: x
DOUBLE PRECISION, DIMENSION(:)		, INTENT(OUT)           :: y

!     Calculates the matrix-vector product:
!        y(:) = y(:) + alpha*A(:,1:A%n)*x(1:A%n)

!     The matrix  A is stored in CSC format
!        [A%n, A%beg, A%jco, A%co]

!     Arguments:
!     ==========
!     A%n      i   Number of rows in matrix  A.
!     A%beg    i   A%beg(c): Location in 'A%jco' and 'A%co' of first
!                  non-zero element in each column  c  of matrix  A.
!     A%jco    i   A%jco(A%beg(c):A%beg(c+1)-1): Row numbers of the
!                  non-zero elements in column  c  of matrix  A.
!     A%co     i   A%co(A%beg(c):A%beg(c+1)-1): Values of the non-zero
!                  elements in column  c  of matrix  A.
!     x        i   Vector, containing segment x(1:A%n).
!     y        io  In:  y(:) Input value vector segment.
!                  Out: y(:) Result value vector segment.

!#enddoc

INTEGER 		:: c, nz
DOUBLE PRECISION 	:: fact

!     Calculate  y(:) = y(:) + alpha*A(:,1:A%n)*x(1:A%n)

DO c = 1, A%n
  fact = alpha*x(c)
  FORALL ( nz = A%beg(c) : A%beg(c+1)-1 )   y(A%jco(nz)) = y(A%jco(nz)) + A%co(nz)*fact
END DO

!     End of  cscvec
END SUBROUTINE cscvec

END MODULE