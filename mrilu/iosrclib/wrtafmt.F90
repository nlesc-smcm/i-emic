!#begindoc

MODULE m_wrtafmt

CONTAINS

SUBROUTINE wrtafmt (filenm, n, a, wmant, ier)

USE m_getunit

CHARACTER (LEN=*)				, INTENT(IN)		:: filenm
INTEGER						, INTENT(IN)            :: n
DOUBLE PRECISION, DIMENSION(1:n,1:n)		, INTENT(IN)         	:: a
INTEGER						, INTENT(IN)            :: wmant
INTEGER						, INTENT(OUT)           :: ier


!     Writes a Full matrix  A  to the ASCII file 'filenm'.
!     The representation of the vector in the file is
!        N,
!        A(i,1:N)
!     Each number is written on a separate line.
!     The order of the numbers in the file are
!        1 * integer:          the number of rows/columns in the
!                              matrix  A
!        N * N * double precision: all values of  A  in row major order

!     Arguments:
!     ==========
!     filenm   i   Filename
!     N        i   The number of rows and columns of the matrix  A.
!     A        i   A(i,j), 1 <= i <= N, 1 <= j <= N, the j-th element
!                  in the i-th row of the matrix  A.
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

!     1997-09-19  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)


INTEGER 		:: i, j, unitnr, wm
CHARACTER (LEN=16) 	:: vfmt

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'wrtamft'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988

!     Structure of the file:
!     1 * integer                 Number of rows/columns of A.
!     N *  a row of the  N  double precision numbers  A(,1:N)

!     Compute the formats

wm = wmant
IF ( wm < 2 ) wm =  2
IF ( wm > 16 ) wm = 16
WRITE ( vfmt, 903 )  wm+6, wm-1

OPEN (unitnr, FILE = filenm, STATUS = 'UNKNOWN', ERR = 997)

REWIND (unitnr, ERR = 997)

WRITE (unitnr, '(I10)', ERR = 995)  n

DO i = 1, n
  WRITE (unitnr, vfmt, ERR = 995)  (a(i,j), j = 1, n)
END DO

CLOSE (unitnr, ERR = 997)

!     Normal return:
RETURN

903  FORMAT ( '(1P,SP,(E', i2, '.', i2, '))' )

!     Error return:
988  ier = ier - 7
995  ier = ier - 2
997  ier = ier - 3

!     End of  wrtafmt
END SUBROUTINE wrtafmt

END MODULE