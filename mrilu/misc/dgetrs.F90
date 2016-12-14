!#begindoc
 
MODULE m_dgetrs

CONTAINS

SUBROUTINE dgetrs(n, a, piv, b )

!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993

!     .. Scalar Arguments ..

USE m_dtrsm
USE m_dlaswp

INTEGER					, INTENT(IN)            :: n
DOUBLE PRECISION, DIMENSION(1:n,1:n)	, INTENT(IN)        	:: a
INTEGER, DIMENSION(1:n)			, INTENT(IN)  		:: piv
DOUBLE PRECISION, DIMENSION(1:n,1:1), INTENT(IN OUT)	        :: b

!  Purpose
!  =======

!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  !!!! Currently B is restricted to a vector, otherwise argument conflict 
!  !!!! with calling routine SOLVE
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.

!  Arguments
!  =========

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.

!  A       (input) DOUBLE PRECISION array, dimension (N,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.

!  PIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row PIV(i).

!  B       (input/output) DOUBLE PRECISION array, dimension (N,1)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.

!#enddoc
!  =====================================================================

!     .. Parameters ..

DOUBLE PRECISION, PARAMETER :: one = 1.0D+0
!     ..

!     Test the input parameters.

!     Quick return if possible
!IF( n == 0 ) RETURN
!IF( n == 0 .OR. 1 == 0 ) RETURN
  
    
! Solve A * X = B.
    
! Apply row interchanges to the right hand sides.
    
  CALL dlaswp( b, 1, n, piv )
    
! Solve L*X = B, overwriting B with X.
    
  CALL dtrsm( .FALSE., .FALSE., a, b )
    
! Solve U*X = B, overwriting B with X.
    
  CALL dtrsm( .TRUE., .TRUE., a, b ) 
  RETURN
  
!     End of DGETRS
  
END SUBROUTINE dgetrs

END MODULE