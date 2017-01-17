!#begindoc

MODULE m_dgeco

CONTAINS
 
SUBROUTINE dgeco (a, n, Pvt, rcond)

USE m_dgefa
USE m_dump

INTEGER					, INTENT(IN)            :: n
DOUBLE PRECISION, DIMENSION(1:n,1:n)	, INTENT(IN OUT)	:: a
INTEGER, DIMENSION(1:n)			, INTENT(OUT)           :: Pvt
DOUBLE PRECISION			, INTENT(OUT)           :: rcond

!     Factors a double precision matrix by gaussian elimination
!     and estimates the condition of the matrix.

!     if  rcond  is not needed, dgefa is slightly faster.
!     to solve    A*x = b ,        follow  dgeco  by  dgesl.
!     to compute  inverse(A)*c ,   follow  dgeco  by  dgesl.
!     to compute  determinant(A) , follow  dgeco  by  dgedi.
!     to compute  inverse(A) ,     follow  dgeco  by  dgedi.

!     Arguments:
!     ==========
!     a        	io  Input:  The matrix A to be factored.
!                   Output: An upper triangular matrix and the
!                   multipliers which were used to obtain it.
!                   The factorization can be written  A = L*U  where
!                   L  is a product of permutation and unit lower
!                   triangular matrices and  U  is upper triangular.
!     n        	i   The order of the matrix A.
!     Pvt     	o   An integer vector of pivot indices.
!     rcond    	o   An estimate of the reciprocal condition of  A.
!                   For the system  A*x = b , relative perturbations
!                   in  A  and  b  of size  epsilon  may cause
!                   relative perturbations in  x  of size  epsilon/rcond.
!                   If  rcond  is so small that the logical expression
!                            1.0 + rcond = 1.0
!                   is true, then  A  may be singular to working
!                   precision.  In particular,  rcond  is zero  if
!                   exact singularity is detected or the estimate
!                   underflows.

!#enddoc

!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.


!     subroutines and functions

!     blas     daxpy, dscal
!     fortran  dabs,  dmax1, dsign

!     internal variables

!     z        -   A work vector whose contents are usually unimportant.
!                  If  A  is close to a singular matrix, then  z  is
!                  an approximate null vector in the sense that
!                  norm(A*z) = rcond*norm(A)*norm(z) .

INTEGER						:: ier
DOUBLE PRECISION 				:: ek, t, wk, wkm
DOUBLE PRECISION 				:: anorm, s, sm, ynorm
INTEGER 					:: info, j, k, kp1, l
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: z

!     compute 1-norm of a

anorm = 0.0D0
DO  j = 1, n
  anorm = DMAX1 (anorm, SUM(ABS(a(1:n,j))) )
END DO

!     factor

CALL dgefa (a(1:n,1:n), Pvt, info)

!     rcond = 1/(norm(A)*(estimate of norm(inverse(A)))) .
!     estimate = norm(z)/norm(y) where  A*z = y  and  trans(A)*y = e .
!     trans(a)  is the transpose of A .  The components of  e  are
!     chosen to cause maximum local growth in the elements of w  where
!     trans(u)*w = e .  The vectors are frequently rescaled to avoid
!     overflow.

!     solve trans(u)*w = e

  ALLOCATE( z(1:n), STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

ek = 1.0D0
z = 0.0D0

DO  k = 1, n
  IF (z(k) /= 0.0D0) ek = DSIGN(ek,-z(k))
    IF (ABS(ek-z(k)) > ABS(a(k,k))) THEN
      s = ABS(a(k,k))/ABS(ek-z(k))
      z = s * z
      ek = s*ek
    END IF
    wk = ek - z(k)
    wkm = -ek - z(k)
    s = ABS(wk)
    sm = ABS(wkm)
    IF (a(k,k) /= 0.0D0) THEN
      wk = wk/a(k,k)
      wkm = wkm/a(k,k)
    ELSE
      wk = 1.0D0
      wkm = 1.0D0
    END IF
    kp1 = k + 1
    sm = sm + SUM(ABS(z(kp1:n) + wkm*a(k,kp1:n)))
    z(kp1:n) = z(kp1:n) + wk*a(k,kp1:n)
    s = s + SUM(ABS(z(kp1:n)))
    IF (s < sm) THEN
      z(kp1:n) = z(kp1:n) + (wkm - wk)*a(k,kp1:n)
      wk = wkm
    END IF
    z(k) = wk
  END DO
  s = 1.0D0/SUM(ABS(z))
  z = s * z
  
!     solve trans(l)*y = w
  
  DO  k = n, 1 , -1
    IF (k < n) z(k) = z(k) + DOT_PRODUCT( a(k+1:n,k), z(k+1:n) )
      IF (ABS(z(k)) > 1.0D0) THEN
        s = 1.0D0/ABS(z(k))
        z = s * z
      END IF
      l = Pvt(k)
      t = z(l)
      z(l) = z(k)
      z(k) = t
    END DO
    s = 1.0D0/SUM(ABS(z))
    z = z * s
    
    ynorm = 1.0D0
    
!     solve l*v = y
    
  DO  k = 1, n
    l = Pvt(k)
    t = z(l)
    z(l) = z(k)
    z(k) = t
    IF (k < n) z(k+1:n) = z(k+1:n) + t * a(k+1:n,k)
      IF (ABS(z(k)) > 1.0D0) THEN
        s = 1.0D0/ABS(z(k))
        z = z * s
        ynorm = s*ynorm
      END IF
    END DO
    s = 1.0D0/SUM(ABS(z))
    z = z * s
    ynorm = s*ynorm
    
!     solve  u*z = v
      
  DO  k = n, 1, -1
    IF (ABS(z(k)) > ABS(a(k,k))) THEN
      s = ABS(a(k,k))/ABS(z(k))
      z = z * s
      ynorm = s*ynorm
    END IF
    IF (a(k,k) /= 0.0D0) z(k) = z(k)/a(k,k)
    IF (a(k,k) == 0.0D0) z(k) = 1.0D0
    z(1:k-1) = z(1:k-1) -z(k) * a(1:k-1,k)
  END DO
!     make znorm = 1.0
  s = 1.0D0/SUM(ABS(z))
  z = z * s
  ynorm = s*ynorm
      
  IF (anorm /= 0.0D0) rcond = ynorm/anorm
  IF (anorm == 0.0D0) rcond = 0.0D0

  DEALLOCATE( z, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

END SUBROUTINE dgeco

END MODULE

