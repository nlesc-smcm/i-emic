!#begindoc
 
MODULE m_rdvec

CONTAINS

SUBROUTINE rdvec (filenm, binary, n, v, ier, errinf)

USE m_getunit

CHARACTER (LEN=*)			, INTENT(IN)            :: filenm
LOGICAL					, INTENT(IN)		:: binary
INTEGER					, INTENT(OUT)           :: n
DOUBLE PRECISION, DIMENSION(1:)		, INTENT(OUT)		:: v
INTEGER					, INTENT(OUT)           :: ier
INTEGER					, INTENT(OUT)           :: errinf

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
!                  = -1   Length of vector >UBOUND(v,1), the actual length
!                         is returned in 'N' and in 'errinf'.
!                  = -3   Error during closing, opening or rewinding
!                         file 'filenm'.
!                  = -4   Unexpected End Of File encountered.
!                  = -5   Error during reading  'filenm'.
!                  =-12   No logical unit available to read from
!                         'filenm'.
!     errinf   o   Extra error information for 'ier' = -1 else
!                  undefined.


!#enddoc

!     1997-11-12  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!     Local variables:
!     ================
INTEGER :: i, unitnr

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'rdvec'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

errinf = 0

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
IF (binary) THEN
  READ (unitnr, END = 996, ERR = 995)  n
ELSE
  READ (unitnr, *, END = 996, ERR = 995)  n
END IF
IF (n > UBOUND(v,1)) THEN
  errinf = n
  GO TO 999
END IF
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
997  ier = ier - 2
999  ier = ier - 1

!     End of  rdvec
END SUBROUTINE rdvec

END MODULE