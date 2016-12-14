!#begindoc
 
MODULE m_dlaswp

CONTAINS

SUBROUTINE dlaswp( a, k1, k2, ipiv )

#ifdef WITH_ATLAS
EXTERNAL 		:: dswap
#else
USE m_dswap
#endif

!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992

!     .. Scalar Arguments ..

DOUBLE PRECISION, DIMENSION(1:,:), INTENT(IN OUT)         :: a
INTEGER, INTENT(IN)                      :: k1
INTEGER, INTENT(IN)                      :: k2
INTEGER, DIMENSION(1:), INTENT(IN)                      :: ipiv

!     ..
!     .. Array Arguments ..


!     ..

!  Purpose
!  =======

!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.

!  Arguments
!  =========

!  N       (input) INTEGER
!          The number of columns of the matrix A.

!  A       (input/output) DOUBLE PRECISION array, dimension (:,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.

!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.

!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.

!  IPIV    (input) INTEGER array, dimension (M*abs(1))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.


!#enddoc
! =====================================================================

!     .. Local Scalars ..
INTEGER          :: i, ip, j
DOUBLE PRECISION :: temp

!     Interchange row I with row IPIV(I) for each of rows K1 through K2.

!$OMP PARALLEL DO PRIVATE(i,j,temp)
!DEC$ PARALLEL
DO  i = k1, k2
  ip = ipiv( i )
  IF ( ip /= i ) THEN
    DO  j = 1, UBOUND(a, 2)
      temp = a(i,j)
      a(i,j) = a(ip,j)
      a(ip,j) = temp
    END DO
  END IF
END DO
!$OMP END PARALLEL DO
  
RETURN
  
END SUBROUTINE dlaswp

END MODULE