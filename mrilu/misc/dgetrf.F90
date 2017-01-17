!#begindoc
 
MODULE m_dgetrf

CONTAINS

SUBROUTINE dgetrf( a, piv, info )

!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993

!     .. Scalar Arguments ..

USE m_dump
USE m_dtrsm
USE m_dlaswp
USE m_dgetf2

DOUBLE PRECISION, DIMENSION(:,:)	, INTENT(IN OUT)        :: a
INTEGER, DIMENSION(:)			, INTENT(OUT)           :: piv
INTEGER					, INTENT(OUT)           :: info

!     ..
!     .. Array Arguments ..


!     ..

!  Purpose
!  =======

!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.

!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).

!  This is the right-looking Level 3 BLAS version of the algorithm.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (M,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.

!  PIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row PIV(i).

!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.

!#enddoc
!  =====================================================================

!     .. Parameters ..

DOUBLE PRECISION, PARAMETER :: one = 1.0D+0
!     ..
!     .. Local Scalars ..
INTEGER :: i, iinfo, j, jb, jjb, nb, minmn, m, n

m = UBOUND(a,1)
n = UBOUND(a,2)
minmn = MIN(m,n)
info=0

IF (UBOUND(piv,1)<minmn) CALL dump(__FILE__,__LINE__,'Upperbound array piv too small')

!     Test the input parameters.

!     Quick return if possible

IF( m == 0 .OR. n == 0 ) RETURN

!     Determine the block size for this environment.

! nb = ilaenv( 1, 'DGETRF', m, n, -1 )
! ilaenv always produces the result 64

nb = 64
!workaround to avoid blocked code, which contains a bug
nb=1
PRINT *,m,n
IF( nb <= 1 .OR. nb >= minmn ) THEN
  
!        Use unblocked code.
  
  CALL dgetf2( a, piv, info )
  print *,info
ELSE
  
!        Use blocked code.
  
  DO  j = 1, minmn , nb
    jb = MIN( minmn-j+1, nb )
    jjb = j+jb
    
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
    
    CALL dgetf2( a( j:m, j:jjb-1 ), piv( j:jjb-1 ), iinfo )
    
!           Adjust INFO and the pivot indices.
    
    IF( info == 0 .AND. iinfo > 0 ) info = iinfo + j - 1
    PRINT *,j,info,iinfo

    FORALL  (i=j: MIN( m, jjb-1 ))  piv( i ) = j - 1 + piv( i )
    
!           Apply interchanges to columns 1:J-1.
    
    CALL dlaswp( a(:,1:j-1), j, jjb-1, piv(1:minmn) )
    
    IF( jjb <= n ) THEN
      
!              Apply interchanges to columns JJB:N.
      
      CALL dlaswp( a(:, jjb:n ), j, jjb-1, piv(1:minmn) )
      
!              Compute block row of U.
      
!     CALL dtrsm( .FALSE., .FALSE., jb, n-j-jb+1, a( j:j+m-1, j:j+m-1 ), a( j:j+m-1, jjb:jjb+n-1 ) )
!      CALL dtrsm( .FALSE., .FALSE., a( j:n-jb+m, j:j+jb-1 ), a( j:j+jb-1, jjb:jjb+jb-1 ) )
      CALL dtrsm( .FALSE., .FALSE., a( j:jjb-1, j:jjb-1 ), a( j:jjb-1, jjb:n) )
      IF( jjb <= m ) THEN
        
!                 Update trailing submatrix.
        
        a( jjb:m, jjb:n ) = a( jjb:m, jjb:n ) - a( jjb:m, j:jjb-1 ) *  a( j:jjb-1, jjb:n ) 
      END IF
    END IF
  END DO
END IF
RETURN

!     End of DGETRF

END SUBROUTINE dgetrf

END MODULE