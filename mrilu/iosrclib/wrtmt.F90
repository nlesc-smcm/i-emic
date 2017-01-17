!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_wrtmt

CONTAINS

SUBROUTINE wrtmt (filenm, binary, a, wmant, ier)

USE m_getunit
USE m_build

CHARACTER (LEN=*), 			INTENT(IN)	:: filenm
LOGICAL,				INTENT(IN)	:: binary
TYPE (csrmatrix),			POINTER     	:: a
INTEGER, 				INTENT(IN)      :: wmant
INTEGER, 				INTENT(OUT)     :: ier

!     Writes a sparse matrix A in Compressed Sparse Row (CSR) format
!     to the ASCII file 'filenm'.  The representation of the matrix
!     is stored, rowwise, in the 4 parameters
!        A%N, a%beg, a%jco and a%co
!     Each number is written on a separate line.
!     The order of the numbers in the file are
!        integer:        A%N, the order of the matrix A.
!        A%N+1 * integer:  Indices indicating the beginning of a row
!        A%NNZ * integer:  column numbers, with A%NNZ
!                        the number of nonzero elements in the
!                        matrix A.
!        A%NNZ * double precision: the values of the nonzero elements
!                        in the matrix A.

!     Arguments:
!     ==========
!     filenm   i   Filename
!     A%N        i   The order of the matrix A
!     a%beg     i   a%beg(r), 1<=r<=A%N: Index in 'a%jco' and 'a%co' of the
!                  first nonzero element A(r,j) in row  r in the
!                  matrix A.
!                  a%beg(A%N+1): Index in 'a%jco' and 'a%co' of last nonzero
!                  element + 1.
!     a%jco     i   a%jco(nz), a%beg(r)<=nz<a%beg(r+1), 1<=r<=A%N: Column
!                  number of the nonzero element a%co(nz) in row  r  of
!                  the matrix  A.
!     a%co      i   a%co(a%beg(r):a%beg(r+1)-1), 1<=r<=A%N: Values of the
!                  nonzero elements in row  r  of the matrix  A.
!     wmant    i   Number of digits in mantisse to represent the value
!                  of a%co(v).  (2 <= wmant <= 16)
!     ier      o   The error code,
!                  =  0   No error
!                  = -3   Error during closing, opening or rewinding
!                         file 'filenm'
!                  = -5   Error during writing 'filenm'
!                  =-12   No logical unit available to write on
!                         'filenm'.

!#enddoc

!     1996-07-10  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)


INTEGER 		:: i, nrd, unitnr, v, wm
CHARACTER (LEN=8) 	:: jcofmt
CHARACTER (LEN=8) 	:: begfmt
CHARACTER (LEN=16) 	:: cofmt

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'wrtmt'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988

!     Structure of the file:
!     1    * integer                 Actual number of rows A%N
!     A%N+1  * integer
!     A%NNZ  * integer
!     A%NNZ  * double precision

!     Compute the formats

nrd     = 1 + INT(LOG10(REAL(a%nnz+1)))
IF ( nrd <= 9 ) THEN
  WRITE ( begfmt, 901 )  nrd
ELSE
  WRITE ( begfmt, 902 )  nrd
END IF

nrd     = 1 + INT(LOG10(REAL(a%n+2)))
IF ( nrd <= 9 ) THEN
  WRITE ( jcofmt, 901 )  nrd
ELSE
  WRITE ( jcofmt, 902 )  nrd
END IF

wm = wmant
IF ( wm < 2 ) wm =  2
IF ( wm > 16 ) wm = 16
WRITE ( cofmt, 903 )  wm+6, wm-1

IF (binary) THEN
  OPEN (unitnr, FILE = filenm, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', ERR = 997)
ELSE
  OPEN (unitnr, FILE = filenm, STATUS = 'UNKNOWN', ERR = 997)
END IF

REWIND (unitnr, ERR = 997)


IF (binary) THEN

  WRITE (unitnr, ERR = 995)  a%n

  WRITE (unitnr, ERR = 995)  (a%beg(i), i = 1, a%n+1)

  WRITE (unitnr, ERR = 995)  (a%jco(v), v = 1, a%nnz)

  WRITE (unitnr, ERR = 995)  (a%co(v), v = 1, a%nnz)

ELSE

  WRITE (unitnr, '(I10)', ERR = 995)  a%n

  WRITE (unitnr, begfmt, ERR = 995)  (a%beg(i), i = 1, a%n+1)

  WRITE (unitnr, jcofmt, ERR = 995)  (a%jco(v), v = 1, a%nnz)

  WRITE (unitnr, cofmt, ERR = 995)  (a%co(v), v = 1, a%nnz)

END IF

CLOSE (unitnr, ERR = 997)

!     Normal return:
RETURN

901  FORMAT ( '((I', i1, '))' )
902  FORMAT ( '((I', i2, '))' )
903  FORMAT ( '(1P,SP,(E', i2, '.', i2, '))' )

!     Error return:
988  ier = ier - 7
995  ier = ier - 2
997  ier = ier - 3

!     End of  wrtmt
END SUBROUTINE wrtmt

END MODULE