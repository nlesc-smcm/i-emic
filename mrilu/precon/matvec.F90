!#begindoc

#ifndef WITH_UNION

#define scdematrix  anymatrix
#define scbmmatrix  anymatrix

#endif

MODULE m_matvec

CONTAINS
 
SUBROUTINE matvec (n, alpha, S, x, y)

USE m_build
USE m_dump
USE m_csrvec
USE m_cscvec
USE m_diavec

INTEGER					, INTENT(IN)    :: n
DOUBLE PRECISION			, INTENT(IN)	:: alpha
TYPE (anymatrix)			, POINTER	:: S
DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)	:: x
DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT):: y

!     Computes the matrix vector product
!        y := y + alpha S x

!     Arguments:
!     ==========
!     N        	i   Number of rows/columns in matrix  S.
!     S      	i   Location of descriptor of the Schur-complement matrix  S.
!     x        	i   Input vector
!     y       	io  Matrix vector product:  S x

!#enddoc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'matvec'

DOUBLE PRECISION, PARAMETER :: one = 1.0D0
DOUBLE PRECISION, PARAMETER :: zero = 0.0D0

!     Local variables:
!     ================

INTEGER						:: ier
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) 	:: wrk, z
TYPE (scdematrix), POINTER           		:: xscde
TYPE (scbmmatrix), POINTER  			:: xa

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

IF (S%typ /= diatp) THEN
!        Initialize:  y := 0
  y = zero
END IF

IF (S%typ == csctp) THEN
  
!   Calculate:  y := y + alpha S x
    CALL cscvec (alpha, anytocsc (S), x, y )
  
ELSE IF (S%typ == csrtp) THEN
  
!   Calculate:  y := y + alpha S x
    CALL csrvec (alpha, anytocsr(S), x, y )

ELSE IF (S%typ == diatp) THEN

!   Calculate:  y := y + alpha S x
    CALL diavec (.true., alpha, anytodia(S), x, y)
  
ELSE IF (S%typ == scbmtp) THEN
    xa => anytoscbm(S)
#ifdef DEBUG
  IF ( n + xa%g /= xa%n ) CALL dump(__FILE__,__LINE__,'Internal error')
#endif
  
!   Request workspace for DP work array, y:
    ALLOCATE( wrk(1:xa%a12%n), STAT=ier )
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
  
!   Request workspace for DP work array, z:
    ALLOCATE( z(1:xa%a12%n), STAT=ier )
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
 
!   Calculate:  y := y + alpha S x  where S = (A22 - A21 inv(A11) A12)
!   Calculate:  y := y + alpha (A22 - A21 inv(A11) A12) x
  
!        (1) Calculate:  y := y + alpha A22 x
  
    CALL csrvec (alpha, xa%a22, x, y(1:xa%a22%n) )
  
!        (2) Calculate:  wrk := A12 x
  
    wrk = zero
    CALL csrvec (one, xa%a12, x, wrk )
  
!        (3) Calculate:  z := inv(A11) wrk  (= inv(A11) A12 x)
  
#ifdef DEBUG
  IF (n + xa%a11d%n /= xa%n) CALL dump(__FILE__,__LINE__,'Internal error')
#endif
  
    CALL diavec (.false., one, xa%a11d, wrk, z )
  
!        (4) Calculate:  y := y - alpha A21 z = alpha A22 y - alpha A21 inv(A11) A12 x
  
    CALL csrvec (-alpha, xa%a21, z, y(1:xa%a21%n) )
  
!   Free workspace segments for work arrays y, z:
    DEALLOCATE( wrk, STAT=ier )
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
    DEALLOCATE( z, STAT=ier )
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
    
  ELSE IF (S%typ == csrdtp) THEN
      xscde => anytoscde(S)
    
!     Calculate:  y := y + alpha S x
    
!     First the diagonal part:
      CALL diavec (.true., alpha, xscde%dia, x, y)
    
!     Next add the off-diagonal part:
      CALL csrvec (alpha, xscde%offd, x, y)
    
  ELSE
    PRINT '(/, A, 2X, A, A, /, 3X, A, I11, 3X, A)' , 'Internal error in', rounam, '!', 'Storage type', S%typ, 'not implemented!'
  END IF

END SUBROUTINE
  
END MODULE