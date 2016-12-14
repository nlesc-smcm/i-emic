!#begindoc
 
MODULE m_rdvc

CONTAINS

SUBROUTINE rdvc (filenm, binary, n, v, ier)

USE m_dump
USE m_getunit

CHARACTER (LEN=*)			, INTENT(IN)            :: filenm
LOGICAL					, INTENT(IN)		:: binary
INTEGER					, INTENT(OUT)           :: n
DOUBLE PRECISION, DIMENSION(:)		, POINTER		:: v
INTEGER					, INTENT(OUT)           :: ier

!     Reads a vector V from the ASCII file 'filenm'.
!     The representation of the vector in the file is
!        N, v
!     Numbers in the file should be separated by at least one space or
!     end-of-line!
!     The order of the numbers in the file are
!        1 * integer:                the length of vector  V
!        End-of-line
!        N * double precision:       all values of V(1:N)

!     Arguments:
!     ==========
!     filenm   i   Filename
!     N        o   The actual length of the vector V.
!     v        o   v(i), 1 <= i <= N, i-th element of vector V.
!     ier      o   The error code,
!                  =  0   No error
!                  = -3   Error during closing, opening or rewinding
!                         file 'filenm'.
!                  = -4   Unexpected End Of File encountered.
!                  = -5   Error during reading  'filenm'.
!                  =-12   No logical unit available to read from
!                         'filenm'.



!#enddoc

!     1997-11-12  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!     Local variables:
!     ================
INTEGER :: i, unitnr
DOUBLE PRECISION :: dummy

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'rdvec'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif



!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988

n   = 0

!     Structure of the file:
!     1    * integer                 Actual length of vector  V.
!     N    * double precision

IF (binary) THEN
  OPEN (unitnr, FILE = filenm, STATUS = 'OLD', FORM = 'UNFORMATTED', ERR = 997)
ELSE
  OPEN (unitnr, FILE = filenm, STATUS = 'OLD', ERR = 997)
END IF

REWIND (unitnr, ERR = 997)
!     Count number of elements:
n = 0
10   CONTINUE
IF (binary) THEN
  READ (unitnr, END = 100, ERR = 995)  dummy
ELSE 
  READ (unitnr, *, END = 100, ERR = 995)  dummy
END IF
n = n + 1
GO TO 10
100  CONTINUE
!     End of file reached
PRINT *,n

ALLOCATE( v(1:n), STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

REWIND (unitnr, ERR = 997)
IF (binary) THEN
  READ (unitnr, END = 996, ERR = 995)  (v(i), i = 1, n)
ELSE
  READ (unitnr, *, END = 996, ERR = 995)  (v(i), i = 1, n)
END IF
CLOSE (unitnr, ERR = 997)

!     Normal return:
RETURN

!     Error return:
988  ier = ier - 7
995  ier = ier - 1
996  ier = ier - 1
997  ier = ier - 3


!     End of  rdvc
END SUBROUTINE rdvc

END MODULE