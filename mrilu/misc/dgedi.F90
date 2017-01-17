!#begindoc
 
MODULE m_dgedi

CONTAINS

SUBROUTINE dgedi (a, n, Pvt)

USE m_dump
#ifdef WITH_ATLAS
EXTERNAL 		:: dswap
#else
USE m_dswap
#endif

INTEGER					, INTENT(IN)		:: n
DOUBLE PRECISION, DIMENSION(1:n,1:n)	, INTENT(IN OUT)        :: a
INTEGER, DIMENSION(1:n)			, INTENT(IN)            :: Pvt

!     Computes the inverse of a matrix,
!     using the factors computed by  dgeco  or  dgefa.

!     Arguments:
!     ==========
!     a        io  Input:  the output from  dgeco  or  dgefa.
!                  Output: inverse of original matrix if requested,
!                          otherwise unchanged.
!     n        i   The order of the matrix  A .
!     Pvt     i   The pivot vector from  dgeco  or  dgefa.
!     det      o   Determinant of original matrix if requested,
!                  otherwise not referenced.
!                     determinant = det(1) * 10.0**det(2)
!                  with  1.0 .le. dabs(det(1)) .lt. 10.0
!                  or
!                     determinant = det(1) = 0.0 .
!     work     -   Work vector; contents destroyed.
!                  = 10   Compute determinant only.

!     error condition

!     A division by zero will occur if the input factor contains
!     a zero on the diagonal and the inverse is requested.
!     It will not occur if the subroutines are called correctly
!     and if
!        dgeco  has set  rcond .gt. 0.0
!     or
!        dgefa  has set  info .eq. 0 .

!#enddoc

!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.

!     subroutines and functions

!     internal variables

INTEGER						:: ier
INTEGER 					:: i, j, k, l
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: work
DOUBLE PRECISION				:: fact

DO  k = 1, n
  fact = 1.0D0/a(k,k)
  a(k,k) = fact
  a(1:k-1,k) = -fact * a(1:k-1,k)
  FORALL  (j = k + 1: n)
    a(1:k-1,j) = a(1:k-1,j) + a(k,j) * a(1:k-1,k)
  END FORALL
  a(k,k+1:n) = a(k,k+1:n) * fact
END DO

!        form inverse(u)*inverse(l)

IF (n>1) THEN

  ALLOCATE( work(1:n), STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

  DO  k = n-1, 1, -1
    work(k+1:n) = a(k+1:n,k)
    a(k+1:n,k) = 0.0D0
!    DO  j = k+1, n
!      a(:,k) = a(:,k) + a(:,j)*work(j)
!    END DO
      a(:,k) = a(:,k) + MATMUL( a(:,k+1:n), work(k+1:n))
    l = Pvt(k)
    IF (l /= k) CALL dswap (a(:,k), a(:,l) )
  END DO

  DEALLOCATE( work, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

END IF

END SUBROUTINE dgedi

END MODULE