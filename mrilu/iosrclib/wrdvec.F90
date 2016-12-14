!#begindoc

MODULE m_wrdvec

CONTAINS
 
SUBROUTINE wrdvec (filenm, binary, n, V, ier)

USE m_build
USE m_getunit
USE m_ioerrmsg

CHARACTER (LEN=*)		, INTENT(IN)           		:: filenm
LOGICAL				, INTENT(IN)			:: binary
INTEGER				, INTENT(OUT)                   :: n
DOUBLE PRECISION, DIMENSION(:)	, POINTER			:: V
INTEGER				, INTENT(OUT)                  	:: ier

!     Reads a vector V from the ASCII file 'filenm'. This vector is
!     stored into the newly allocated segment V.

!     The representation of the vector in the file is
!        N, v
!     Numbers in the file should be separated by at least one space or
!     end-of-line!
!     The order of the numbers in the file are
!        1 * integer:                N: the length of vector  V
!        End-of-line
!        N * double precision:       all values of V(1:N)

!     Arguments:
!     ==========
!     filenm   	i   Filename
!     N        	o   The actual length of the vector V.
!     V      	o   Location of the vector V.
!              	    V(i), 1 <= i <= N, i-th element of vector V.
!     ier      	o   Error code:
!              	    =  0  No error occurred, subroutine allocated storage
!              	          as described above.
!              	    <  0  Some error occurred, an error message has been
!              	          written to standard output.

!#enddoc

!     1997-11-12  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!     Local Parameter:
!     ================

CHARACTER (LEN=*), PARAMETER :: rounam = 'wrdvec'

!     Local variables:
!     ================
INTEGER :: errinf, i, unitnr

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988

!     Structure of the file:
!     1    * integer                 Actual length of vector: N
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
IF (n <= 0) THEN
  PRINT '(/, A, 2X, A, /, 3X, A, I11)', 'Fatal error occurred in', rounam,  &
      'Illegal value of number of elements in vector:', n
  GO TO 1000
END IF


!     Request a segment to store the vector:

ALLOCATE( V(1:n), STAT=ier)
IF (ier < 0) GO TO 1000
IF (binary) THEN
  READ (unitnr, END = 996, ERR = 995) (V(i), i = 1, n)
ELSE
  READ (unitnr, *, END = 996, ERR = 995) (V(i), i = 1, n)
END IF
CLOSE (unitnr, ERR = 997)

!     Normal return:
ier = 0
RETURN

!     IO-error:
988  ier = ier - 7
995  ier = ier - 1
996  ier = ier - 1
997  ier = ier - 3

errinf = 0

CALL ioerrmsg (filenm, rounam, 'wrdvec', ier, errinf)


!     Error Return, force value ier < 0:
1000 CONTINUE
ier = -1

!     End of  wrdvec
END SUBROUTINE wrdvec

END MODULE