!#begindoc
 
MODULE m_filgrd

CONTAINS

SUBROUTINE filgrd (m, strtch, twosid, shiftlb, shiftrb, x, h)

INTEGER					, INTENT(IN)		:: m
DOUBLE PRECISION			, INTENT(IN)         	:: strtch
LOGICAL					, INTENT(IN)            :: twosid
LOGICAL					, INTENT(IN)            :: shiftlb
LOGICAL					, INTENT(IN)            :: shiftrb
DOUBLE PRECISION, DIMENSION(0:m)	, INTENT(OUT)           :: x
DOUBLE PRECISION, DIMENSION(0:m+1)	, INTENT(OUT)           :: h

!     Compute the positions of grid points on the interval [0,1].

!     Arguments:
!     ==========
!     M        i   The number of intervals in the interval.
!     Strtch   i   The stretching factor.
!     TwoSid   i   The type of stretching:
!                  .TRUE.   Two sided stretching with
!                           h((M+1)/2) = h(1) = Strtch,
!                           h(M+1-i) = h(i) for i = 1,2,...,(M+1)/2
!                  .FALSE.: Single side stretching with
!                           h(M) / h(1) = Strtch
!     ShiftLB  i   Shift Left Border between gridpoints
!                  .TRUE.:  Left border, 0, at:  (x(0)+x(1))/2
!                  .FALSE.: Left border, 0, at:  x(0)
!     ShiftRB  i   Shift Right Border between gridpoints
!                  .TRUE.:  Right border, 1, at:  (x(M-1)+x(M))/2
!                  .FALSE.: Right border, 1, at:  x(M)
!     x        o   Positions of the grid points in the interval [0,1]
!                  as indicated by the arguments above.
!     h        o   Sizes of the subintervals, i.e.:
!                     h(1:M) = x(1:M) - x(0:M-1),
!                     h(0) = h(1), h(M+1) = h(M)

!#enddoc

!     Local variables:
!     ================

INTEGER 			:: i, mid
DOUBLE PRECISION 		:: c, fact, lb, rb, nlb, nrb


!     Initialize left- and right boundary:
lb = 0.0D0
rb = 1.0D0

IF (twosid) THEN
!        Two sided stretching with:  h((M+1)/2) = h(1) = Strtch
  
  mid = (m+1) / 2

  c = EXP(LOG(strtch) / DBLE(mid-1))
  x(0) = 0
  DO i = 1, mid
    x(i) = 1 + c * x(i-1)
  END DO
!        {  x(i) = SUM(j :0<=j<i: c**j), i = 0,1,2,...,mid  }
  
  IF (MOD(m,2) == 0) THEN
!           The midpoint lb+(rb-lb)/2 should coincide with x(mid).
    
    fact = (rb - lb) / (2 * x(mid))
!           {  lb + fact * x(mid) = lb + (rb-lb)/2  ,
!              lb + fact * x(0)   = lb              }
    
    x(mid) = lb + fact * x(mid)
  ELSE
!           The midpoint lb+(rb-lb)/2 should coincide with
!           (x(mid)+x(mid-1))/2.
    
    fact = (rb - lb) / (x(mid)+x(mid-1))
!           {  lb + fact * (x(mid)+x(mid-1))/2 = lb + (rb-lb)/2  }
  END IF
  
  DO i = 0, mid-1
    c      = fact * x(i)
    x(i)   = lb + c
    x(m-i) = rb - c
  END DO

ELSE
!        Single side stretching with:  h(M) / h(1) = Strtch
  
  c = EXP(LOG(strtch) / DBLE(m-1))
  x(0) = 0
  DO i = 1, m
    x(i) = 1 + c * x(i-1)
  END DO
!        {  x(i) = SUM(j :0<=j<i: c**j), i = 0,1,2,...,M  }
  
  fact = (rb - lb) / x(m)
!        {  lb + fact * x(M) = rb  ,  lb + fact * x(0) = lb  }
  
  FORALL (i =0:m) x(i) = lb + fact * x(i)
END IF

IF (shiftlb .OR. shiftrb) THEN
  nlb = lb
  IF (shiftlb) THEN
    nlb = (lb + x(1)) / 2.0D0
  END IF
  nrb = rb
  IF (shiftrb) THEN
    nrb = (x(m-1) + rb) / 2.0D0
  END IF
  
!        Transform [nlb,nrb] -> [lb,rb]:
  fact = (rb - lb) / (nrb - nlb)
  x(0) = lb - fact * (nlb - lb)
  
  FORALL (i=1:m) x(i) = x(0) + fact * (x(i) - lb)

END IF


FORALL (i=1:m)  h(i) = x(i) - x(i-1)

h(0)   = h(1)
h(m+1) = h(m)


!     End of  filgrd
END SUBROUTINE filgrd

END MODULE