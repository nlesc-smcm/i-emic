!#begindoc
 
MODULE m_wrtvec

CONTAINS

SUBROUTINE wrtvec (filenm, binary, n, v, wmant, ier)

USE m_getunit

CHARACTER (LEN=*)			, INTENT(IN)		:: filenm
LOGICAl					, INTENT(IN)		:: binary
INTEGER					, INTENT(IN)            :: n
DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: v
INTEGER					, INTENT(IN)            :: wmant
INTEGER					, INTENT(OUT)           :: ier

!     Writes a vector V  to the ASCII file 'filenm'.
!     The representation of the vector in the file is
!        N, v
!     Each number is written on a separate line.
!     The order of the numbers in the file are
!        1 * integer:          the length of the vector v
!        N * double precision: all values of v

!     Arguments:
!     ==========
!     filenm   i   Filename
!     N        i   The length of the vector v
!     v        i   v(i), 1 <= i <= N, i-th element of the vector v
!     wmant    i   Number of digits in mantisse to represent the value
!                  of v(i).  (2 <= wmant <= 16)
!     ier      o   The error code,
!                  =  0   No error
!                  = -3   Error during closing, opening or rewinding
!                         file 'filenm'
!                  = -5   Error during writing 'filenm'
!                  =-12   No logical unit available to write on
!                         'filenm'.

!#enddoc

!     1996-07-12  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)


INTEGER :: i, unitnr, wm
CHARACTER (LEN=16) :: vfmt

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988

!     Structure of the file:
!     1    * integer                 Actual length of vector N
!     N    * double precision

!     Compute the formats

wm = wmant
IF ( wm < 2 ) wm =  2
IF ( wm > 16 ) wm = 16
WRITE ( vfmt, 903 )  wm+6, wm-1
!      PRINT '(A, A, A)', 'Variable format: >', vfmt, '<'

IF (binary) THEN
  OPEN (unitnr, FILE = filenm, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', ERR = 997)
ELSE
  OPEN (unitnr, FILE = filenm, STATUS = 'UNKNOWN', ERR = 997)
END IF

REWIND (unitnr, ERR = 997)

IF (binary) THEN

  WRITE (unitnr, ERR = 995)  n

  WRITE (unitnr, ERR = 995)  (v(i), i = 1, n)

ELSE

  WRITE (unitnr, '(I10)', ERR = 995)  n

  WRITE (unitnr, vfmt, ERR = 995)  (v(i), i = 1, n)

END IF

CLOSE (unitnr, ERR = 997)

!     Normal return:
RETURN

903  FORMAT ( '(1P,SP,(E', i2, '.', i2, '))' )

!     Error return:
988  ier = ier - 7
995  ier = ier - 2
997  ier = ier - 3

!     End of  wrtvec
END SUBROUTINE wrtvec

END MODULE