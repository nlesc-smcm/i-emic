!#begindoc
 
MODULE m_dswap

USE m_dump

CONTAINS

!DEC$ ATTRIBUTES INLINE :: dswap

SUBROUTINE  dswap ( x, y )

DOUBLE PRECISION, DIMENSION(:)	, INTENT(IN OUT)        :: x, y

!     Interchanges the values of the two vectors  x  and  y

!     Arguments:
!     ==========
!     x       io  Input:  Vector  x 
!                 Output: Vector  y
!     y       io  Input:  Vector  y
!                 Output: Vector  x

!#enddoc

DOUBLE PRECISION :: temp
INTEGER		 :: i

IF (UBOUND(x, 1) /= UBOUND(y, 1) ) CALL dump(__FILE__,__LINE__,'Incompatible dimension of x and y')

!$OMP PARALLEL DO PRIVATE(i,temp)
DO  i = 1, UBOUND(x, 1)
  temp = x(i)
  x(i) = y(i)
  y(i) = temp
END DO
!$OMP END PARALLEL DO


END SUBROUTINE  dswap

END MODULE












