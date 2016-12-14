!#begindoc
 
#ifndef WITH_UNION

#define diamatrix  anymatrix

#endif

MODULE m_diavec

CONTAINS
 
SUBROUTINE diavec (add, alpha, A, x, y )

USE m_build

LOGICAL					, INTENT(IN)		:: add
DOUBLE PRECISION			, INTENT(IN)		:: alpha
TYPE (diamatrix)			, POINTER		:: A
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)    	:: x
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN OUT)	:: y

!     Calculates the diagonal-matrix vector product:
!        add = .true.  y(1:A%n) := y(1:A%n) + alpha * A(1:A%n,1:A%n) x(1:A%n)
!        add = .false. y(1:A%n) :=            alpha * A(1:A%n,1:A%n) x(1:A%n)

!     Arguments:
!     ==========

!     A%n     	i   Number of rows/columns in (block-) diagonal matrix A: MOD(A%n,A%blksiz)=0.
!     A%blksiz 	i   Number of rows/columns in the blocks of (block-) diagonal matrix A.
!                  'N' should be an integer multiple of 'A%blksiz'.
!     A%com    	i   Values of the elements of A, stored in column major order per block.
!     x        	i   Input vector.
!     y        	io  Diagonal vector product: 
!			y(1:A%n) := y(1:A%n) + alpha * A(1:A%n,1:A%n) x(1:A%n)

!#enddoc


INTEGER :: rowbeg

IF (add) THEN

!     Calculate:    y(1:A%n) :=  y(1:A%n) + alpha * A(1:A%n,1:A%n) x(1:A%n)

  IF (A%blksiz > 1) THEN
    FORALL(rowbeg = 1:A%n:A%blksiz)  y(rowbeg:rowbeg+A%blksiz-1) = &
    y(rowbeg:rowbeg+A%blksiz-1) + &
    alpha * MATMUL( A%com(:,rowbeg:rowbeg+A%blksiz-1), x(rowbeg:rowbeg+A%blksiz-1) )
  ELSE
    y = y + alpha * A%com(1,:) * x
  END IF

ELSE

!     Calculate:    y(1:A%n) :=  alpha * A(1:A%n,1:A%n) x(1:A%n)

  IF (A%blksiz > 1) THEN
    FORALL (rowbeg = 1:A%n:A%blksiz) y(rowbeg:rowbeg+A%blksiz-1) = &
    alpha * MATMUL( A%com(:,rowbeg:rowbeg+A%blksiz-1), x(rowbeg:rowbeg+A%blksiz-1) )
  ELSE
    y = alpha * A%com(1,:) * x
  END IF

END IF

END SUBROUTINE diavec

END MODULE