!#begindoc
 
MODULE m_dperv

CONTAINS

SUBROUTINE dperv (invper, p, x, px)

USE m_dump

LOGICAL					, INTENT(IN)	:: invper
INTEGER, DIMENSION(1:)			, INTENT(IN)    :: p
DOUBLE PRECISION, TARGET, DIMENSION(1:)	, INTENT(IN)	:: x
DOUBLE PRECISION, TARGET, DIMENSION(1:)	, INTENT(OUT)   :: px

!     Permutes an integer vector 'x' and stores the result into
!     'Px', with:
!        Px = P x      , or
!        Px = inv(P) x
!     depending on the value of 'InvPer'.
!     The permutation matrix P is characterised by the permutation
!     vector 'p'.

!     Arguments:
!     ==========
!     InvPer   i   Use the inverse permutation  inv(P) in stead of P.
!     p        i   Permutation vector; p is a permutation of (1:n)
!     x        i   In:  Some input vector  x_in
!     Px       o   The permuted vector 'x':
!                  InvPer    :  Px = P x_in
!                  NOT InvPer:  Px = inv(P) x_in = P' x_in

!     N.B.
!     The actual arguments corresponding with the vectors 'x' and 'Px'
!     should not overlap!

!#enddoc

!     Local variables:
!     ================
INTEGER 					:: i, ier
DOUBLE PRECISION, DIMENSION(:)	, POINTER       :: ptr_x,ptr_px
LOGICAL						:: copbck

ptr_x  => x
ptr_px => px

IF ( ASSOCIATED(ptr_x, ptr_px) ) THEN 
  ALLOCATE( ptr_px(1:UBOUND(x,1)), STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
  copbck = .true.
ELSE
  copbck = .false.
END IF



IF ( UBOUND(p,1) /= UBOUND(x,1) )   STOP 'Abnormal temination in dperv: incompatible upper bounds of p and x'
IF ( UBOUND(p,1) /= UBOUND(ptr_px,1) )  STOP 'Abnormal temination in dperv: incompatible upper bounds of p and ptr_px'

!     Determine the permuted vector to be stored into  Px:

# ifdef DEBUG

IF (invper) THEN
!        Px := inv(P) x
  DO i = 1, UBOUND(p,1)
    IF ( p(i) < 1            )  STOP 'Abnormal temination in dperv: p(i)<1'
    IF ( p(i) > UBOUND(ptr_px,1) )  STOP 'Abnormal temination in dperv: incompatible p(i) and ptr_px'
    ptr_px(p(i)) = x(i)
  END DO
ELSE
!       Px := P x (= x(p)):
  DO i = 1, UBOUND(p,1)
    IF ( p(i) < 1           )  STOP 'Abnormal temination in dperv: p(i)<1'
    IF ( p(i) > UBOUND(x,1) )  STOP 'Abnormal temination in dperv: incompatible p(i) and x'
    ptr_px(i) = x(p(i))
  END DO
  
END IF

# else

IF (invper) THEN
!        Px := inv(P) x
  FORALL (i=1:UBOUND(p,1)) ptr_px(p(i)) = x(i)
ELSE
!       Px := P x (= x(p)):
  FORALL (i=1:UBOUND(p,1)) ptr_px(i) = x(p(i))
END IF

# endif

!     Copy permuted vector back if necessary:
IF (copbck) THEN
  px = ptr_px
  DEALLOCATE( ptr_px, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
END IF

END SUBROUTINE dperv

END MODULE










