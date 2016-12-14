!#begindoc
 
MODULE m_dgetf2

CONTAINS

SUBROUTINE dgetf2(a, piv, info )

!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992

!     .. Scalar Arguments ..

USE m_dump
#ifdef WITH_ATLAS
EXTERNAL 		:: dswap
#else
USE m_dswap
#endif

DOUBLE PRECISION, DIMENSION(:,:)	, INTENT(IN OUT)        :: a
INTEGER, DIMENSION(:)			, INTENT(OUT)           :: piv
INTEGER					, INTENT(OUT)           :: info

!  Purpose
!  =======

!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.

!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).

!  This is the right-looking Level 2 BLAS version of the algorithm.

!  Arguments
!  =========

!  A       (input/output) DOUBLE PRECISION array, dimension (M,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.

!  PIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row PIV(i).

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.

!#enddoc
!  =====================================================================

INTEGER :: j, jp, m, n, minmn

m = UBOUND(a,1)
n = UBOUND(a,2)
minmn = MIN( m, n )

IF (UBOUND(piv,1)<minmn) CALL dump(__FILE__,__LINE__,'Upperbound array piv too small')

info = 0

!     Quick return if possible

IF( minmn == 0 ) RETURN

DO  j = 1, minmn
  
!        Find pivot and test for singularity.
  
  jp = j - 1 + MAXLOC( a(j:m,j), DIM=1 )
  piv( j ) = jp
  IF( a( jp, j ) /= 0.0D+0 ) THEN
    
!   Apply the interchange to columns 1:N.
    
    IF( jp /= j ) CALL dswap(a(j,:), a(jp,:) )
    
!   Compute elements J+1:M of J-th column.
    
    IF( j < m ) a( j+1:m, j ) = 1.0D+0 / a( j, j ) * a( j+1:m, j )
    
  ELSE IF( info == 0 ) THEN
    
    info = j
  
  END IF
  
  IF( j < minmn ) THEN
    
!    Update trailing submatrix.
    
     a(j+1:m,j+1:n) = a(j+1:m,j+1:n) - a(j+1:m,j:j) * TRANSPOSE(a(j+1:n,j:j))
  
   END IF
END DO
RETURN

!     End of DGETF2

END SUBROUTINE dgetf2

END MODULE