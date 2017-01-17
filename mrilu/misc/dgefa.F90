!#begindoc
 
MODULE m_dgefa

CONTAINS

SUBROUTINE dgefa (a, ipvt, info)

DOUBLE PRECISION, DIMENSION(1:,1:)	, INTENT(IN OUT)	:: a
INTEGER, DIMENSION(1:)			, INTENT(OUT)           :: ipvt
INTEGER					, INTENT(OUT)           :: info

!     Factors a double precision matrix, A, by gaussian elimination.

!     dgefa  is usually called by  dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .

!     Arguments:
!     ==========
!     a        io  Input:  The matrix A to be factored.
!                  Output: An upper triangular matrix, U, and the
!                  multipliers which were used to obtain it.
!                  The factorization can be written  A = L*U  where
!                  L  is a product of permutation and unit lower
!                  triangular matrices and  U  is upper triangular.
!     ipvt     o   An integer vector of pivot indices.
!     info     o   = 0  normal value.
!                  = k  if  u(k,k) = 0.0 .  This is not an error
!                  condition for this subroutine, but it does
!                  indicate that  dgesl  or  dgedi  will divide by zero
!                  if called.  Use  rcond  in  dgeco  for a reliable
!                  indication of singularity.

!#enddoc

!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.

!     internal variables

DOUBLE PRECISION 	:: t
INTEGER 		:: j, k, l


INTEGER 	:: n

n = UBOUND(a,1)

IF ( UBOUND(a,1) /= UBOUND(a,2) ) STOP 'Abnormal temination in dgefa: matrix a is not square'
IF ( UBOUND(a,1) /= UBOUND(ipvt,1) ) STOP 'Abnormal temination in dgefa: incompatible a and ipvt'

!     gaussian elimination with partial pivoting

info = 0
DO  k = 1, n-1
  
! find l = pivot index

  l = MAXLOC(ABS(a(k:n,k)),DIM=1) + k - 1
  ipvt(k) = l
  
! zero pivot implies this column already triangularized
  
  IF (a(l,k) /= 0.0D0) THEN
  
!   interchange if necessary
  
    IF (l /= k)  THEN
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    END IF
  
!   compute multipliers
  
    a(k+1:n,k) =  - (1.0D0 / a(k,k)) * a(k+1:n,k)
  
!   row elimination with column indexing
  
    DO  j = k+1, n
      t = a(l,j)
      IF (l /= k) THEN
        a(l,j) = a(k,j)
        a(k,j) = t
      END IF
      a(k+1:n,j) = a(k+1:n,j) + t * a(k+1:n,k)
    END DO
  ELSE
    info = k
  END IF
END DO
ipvt(n) = n
IF (a(n,n) == 0.0D0) info = n

END SUBROUTINE dgefa

END MODULE




