!#begindoc
 
MODULE m_dtrsm

CONTAINS

SUBROUTINE dtrsm ( upper, nounit, a, b )
!     .. Scalar Arguments ..

USE m_dump

LOGICAL					, INTENT(IN)	:: upper
LOGICAL					, INTENT(IN)    :: nounit
DOUBLE PRECISION, DIMENSION(:,:)	, INTENT(IN)    :: a
DOUBLE PRECISION, DIMENSION(:,:)	, INTENT(OUT)   :: b

!     .. Array Arguments ..

!     ..

!  Purpose
!  =======

!  DTRSM  solves one of the matrix equations

!     A *X = B

!  where X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  

!  The matrix X is overwritten on B.

!  Parameters
!  ==========

!  UPPER   - LOGICAL.
!           On entry, UPPER specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:

!              .TRUE.   A is an upper triangular matrix.

!              .FALSE.  A is a lower triangular matrix.

!           Unchanged on exit.

!  NOUNIT   - LOGICAL
!           On entry, NOUNIT specifies whether or not A is unit triangular
!           as follows:

!              .FALSE.  A is assumed to be unit triangular.

!              .TRUE.   A is not assumed to be unit triangular.

!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( m, m )
!           The upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPPER = .FALSE., the lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  NOUNIT = .FALSE. the nounitonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.

!  B      - DOUBLE PRECISION array of DIMENSION ( m, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.

!  Level 3 Blas routine.


!#enddoc


!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.

!     .. Local Scalars ..
INTEGER :: j, k, n, m

!     .. Executable Statements ..

!     Test the input parameters.

IF (UBOUND(a, 1) /= UBOUND(a, 2) ) CALL dump(__FILE__,__LINE__,'Matrix A is not a square matrix')
IF (UBOUND(a, 1) /= UBOUND(b, 1) ) CALL dump(__FILE__,__LINE__,'Incompatible dimension of A and B')

n = UBOUND(b, 2)
m = UBOUND(b, 1)

!     Quick return if possible.

IF( n == 0 ) RETURN

!           Form  B := inv( A )*B.
    
IF( upper )THEN
  DO   j = 1, n
    DO   k = 1, m
      IF( b( k, j ) /= 0.0D+0 )THEN
        IF( nounit ) b( k, j ) = b( k, j )/a( k, k )
        b( 1:k-1, j ) = b( 1:k-1, j ) - b( k, j )*a( 1:k-1, k )
      END IF
    END DO
  END DO
ELSE
  DO   j = 1, n
    DO  k = 1, m
      IF( b( k, j ) /= 0.0D+0 )THEN
        IF( nounit ) b( k, j ) = b( k, j )/a( k, k )
        b( k+1:m, j ) = b( k+1:m, j ) - b( k, j )*a( k+1:m, k )
      END IF
    END DO
  END DO
END IF

RETURN

END SUBROUTINE dtrsm

END MODULE