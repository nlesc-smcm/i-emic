!#begindoc

MODULE m_dgesl

CONTAINS
 
SUBROUTINE dgesl (a, n, ipvt, b, job)

INTEGER					, INTENT(IN)             :: n
DOUBLE PRECISION, DIMENSION(1:n,1:n)	, INTENT(IN)             :: a
INTEGER, DIMENSION(1:n)			, INTENT(IN)             :: ipvt
DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)         :: b
INTEGER					, INTENT(IN)             :: job

!     Solves the double precision linear system
!        A*x = b
!     or
!        trans(A)*x = b
!     using the factors of A computed by  dgeco  or  dgefa.

!     Arguments:
!     ==========
!     a        i   The output from  dgeco  or  dgefa.
!     lda      i   The leading dimension of the array  a .
!     n        i   The order of the matrix  A .
!     ipvt     i   The pivot vector from  dgeco  or  dgefa.
!     b        io  Input:  The right hand side vector.
!                  Output: The solution vector  x.
!     job      i   .EQ. 0   Solve  A*x = b ,
!                  .NE. 0   Solve  trans(A)*x = b  where
!                           trans(A)  is the transpose of  A.

!     error condition

!     A division by zero will occur if the input factor contains a
!     zero on the diagonal.  Technically this indicates singularity
!     but it is often caused by improper arguments or improper
!     setting of lda.  It will not occur if the subroutines are
!     called correctly and if
!        dgeco  has set  rcond .gt. 0.0
!     or
!        dgefa  has set  info .eq. 0 .

!     Example:
!     To compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns, use:

!        CALL dgeco (a, lda, n, ipvt, rcond, z)
!        if (rcond is too small) go to ...
!        do 10 j = 1, p
!           CALL dgesl (a, lda, n, ipvt, c(1,j), 0)
!     10 continue

!#enddoc

!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.

!     subroutines and functions

!     internal variables

DOUBLE PRECISION 	:: t
INTEGER 		:: k, l

IF (job == 0) THEN
  
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
  
    DO  k = 1, n-1
      l = ipvt(k)
      t = b(l)
      IF (l /= k) THEN
        b(l) = b(k)
        b(k) = t
      END IF
     b(k+1:n) = b(k+1:n) + t * a(k+1:n,k)
    END DO
    
!        now solve  u*x = y
    
    DO  k = n, 1, -1
      b(k) = b(k)/a(k,k)
      b(1:k-1) = b(1:k-1) -b(k) * a(1:k-1,k)
    END DO
  ELSE
    
!        job .NE. 0, solve  trans(a) * x = b
    
!        First solve  trans(u)*y = b
    
    DO  k = 1, n
      b(k) = (b(k) - DOT_PRODUCT( a(1:k,k), b(1:k) ))/a(k,k)
    END DO
    
!        now solve trans(l)*x = y
    
    DO  k = n-1, 1, -1
      b(k) = b(k) + DOT_PRODUCT( a(k+1:n,k), b(k+1:n) )
      l    = ipvt(k)
      IF (l /= k) THEN
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      END IF
    END DO
  END IF
  RETURN
  
END SUBROUTINE dgesl

END MODULE