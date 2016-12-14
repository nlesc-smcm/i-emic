!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_csrvec

CONTAINS

SUBROUTINE csrvec (alpha, a, x, y )

USE m_build

DOUBLE PRECISION			, INTENT(IN)             :: alpha
TYPE (csrmatrix)			, POINTER		 :: a
DOUBLE PRECISION, DIMENSION(:)		, INTENT(IN)             :: x
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN OUT)         :: y

!     Calculates the matrix-vector product:
!        y(1:A%n) = y(1:A%n) + alpha*A(1:A%n,:)*x(:)

!     The matrix  A is stored in CSR format
!        [A%n, A%beg, A%jco, A%co]

!     Arguments:
!     ==========
!     A%n      i   Number of rows in matrix  A.
!     A%beg    i   A%beg(r): Location in 'A%jco' and 'A%co' of first
!                  non-zero element in each row  r  of matrix  A.
!     A%jco    i   A%jco(A%beg(r):A%beg(r+1)-1): Column numbers of the
!                  non-zero elements in row  r  of matrix  A.
!     A%co     i   A%co(A%beg(r):A%beg(r+1)-1): Values of the non-zero
!                  elements in row  r  of matrix  A.
!     x        i   Vector.
!     y        io  In:  y(1:A%n) Input value vector segment.
!                  Out: y(1:A%n) Result value vector segment.

!#enddoc

INTEGER 		:: r

!     Calculate  y(1:A%n) := y(1:A%n) + alpha*A(1:A%n,:)*x(:)

FORALL ( r = 1:A%n )
  y(r) = y(r) + alpha*SUM(A%co(A%beg(r):A%beg(r+1)-1)*x(A%jco(A%beg(r):A%beg(r+1)-1)))
END FORALL

END SUBROUTINE csrvec

END MODULE








