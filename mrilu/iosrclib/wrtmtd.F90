!#begindoc

#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_wrtmtd

CONTAINS

SUBROUTINE wrtmtd (filenm, binary, a, ad, wmant)

USE m_build
USE m_getunit
USE m_dump
 
CHARACTER (LEN=*)				, INTENT(IN)		:: filenm
LOGICAL						, INTENT(IN)		:: binary
TYPE (csrmatrix)				, POINTER		:: a
TYPE (diamatrix)				, POINTER		:: ad
INTEGER						, INTENT(IN)            :: wmant

!     Writes a sparse matrix  A + D  in Compressed Sparse Row (CSR)
!     format to the ASCII file 'filenm'.
!     The (block-) diagonal of A is zero.  The non-zero entries of the
!     matrix A are stored, in CSR format, in:
!        [N, a%beg, a%jco, a%co]
!     The values of the (block-)diagonal matrix  D  are stored in  ad%com,
!     in column major order.

!     Each number is written on a separate line.
!     The order of the numbers in the file are
!        1   * integer:  the order of the matrix A (= N).
!        N+1 * integer:  indices indicating the beginning of a row
!        NNZ * integer:  column numbers of the entries in  D  and the
!                        non-zeros of  A, with
!                        NNZ := a%beg(A%N+1)-1 + A%N*Ad%blksiz
!        NNZ * double precision:the corresponding entries in  D and
!                               the non-zeros of  A.

!     Arguments:
!     ==========
!     filenm   i   Filename
!     a%n        i   The order of the matrix A
!     Ad%blksiz   i   Number of rows/columns in a diagonal block.
!                  'N' should be an integer multiple of 'Ad%blksiz'.
!     a%beg     i   a%beg(i): index in 'a%jco' and 'a%co' of the first
!                  off-diagonal element in row i of matrix A.
!                  a%beg(N+1) = index of last nonzero element + 1
!     a%jco     i   a%jco(nz): column number of off-diagonal element
!                  a%co(nz).
!     a%co      i   a%co(nz): value of off-diagonal element.
!     ad%com     i   Contains the elements of the diagonal blocks,
!                  stored in column major order per block.
!     wmant    i   Number of digits in mantisse to represent the value
!                  of a%co(v).  (2 <= wmant <= 16)

!#enddoc

!     1998-03-06  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'wrtmtd'
#endif

!     Local Variables:
!     ================
INTEGER :: i, j, nrd, nz, unitnr, wm, ier, r
CHARACTER (LEN=8) :: jcofmt
CHARACTER (LEN=8) :: begfmt
CHARACTER (LEN=16) :: cofmt

!     Statement functions:
!     ====================

INTEGER :: firrow, lasrow

!     FirRow(r)  returns the first row number of the block containing
!                row number 'r'.

firrow(r)  = ((r-1)/ad%blksiz)*ad%blksiz + 1

!     LasRow(r)  returns the last row number of the block containing
!                row number 'r'.

lasrow(r)  = ((r-1)/ad%blksiz)*ad%blksiz + ad%blksiz

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam
#endif


!     Get an free logical unit number for writing:
CALL getunit (unitnr, ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error during reading filenm')

!     Structure of the file:
!     1    * integer                 Actual number of rows N
!     N+1  * integer
!     NNZ  * integer
!     NNZ  * double precision

!     Compute the formats

nrd     = 1 + INT(LOG10(REAL(csrnnz(a) + a%n*ad%blksiz+1)))
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

!     Write data to file:

IF (binary) THEN

  WRITE (unitnr, ERR = 995)  a%n

  WRITE (unitnr, ERR = 995) ( a%beg(i) + (i-1)*ad%blksiz, i = 1, a%n+1 )

  WRITE (unitnr, ERR = 995) ( ( j                             , j = firrow(i), lasrow(i) ),  &
  ( a%jco(nz), nz = a%beg(i), a%beg(i+1)-1 ), i = 1, a%n )

  WRITE (unitnr, ERR = 995) ( ( ad%com(MOD(i-1,ad%blksiz)+1,j), j = firrow(i), lasrow(i) ),  &
  ( a%co(nz) , nz = a%beg(i), a%beg(i+1)-1 ), i = 1, a%n )

ELSE

  WRITE (unitnr, '(I10)', ERR = 995)  a%n

  WRITE (unitnr, begfmt, ERR = 995) ( a%beg(i) + (i-1)*ad%blksiz, i = 1, a%n+1 )

  WRITE (unitnr, jcofmt, ERR = 995) ( ( j, j = firrow(i), lasrow(i) ), &
  ( a%jco(nz), nz = a%beg(i), a%beg(i+1)-1 ), i = 1, a%n )

  WRITE (unitnr, cofmt, ERR = 995) ( ( ad%com(MOD(i-1,ad%blksiz)+1,j), &
  j = firrow(i), lasrow(i) ), ( a%co(nz), nz = a%beg(i), a%beg(i+1)-1 ), &
  i = 1, a%n )


END IF

CLOSE (unitnr, ERR = 997)

!     Normal return:
RETURN

901  FORMAT ( '((I', i1, '))' )
902  FORMAT ( '((I', i2, '))' )
903  FORMAT ( '(1P,SP,(E', i2, '.', i2, '))' )

995  CALL dump(__FILE__,__LINE__,'No logical unit available to write on filenm')
997  CALL dump(__FILE__,__LINE__,'Error during closing, opening or rewinding file filenm')

END SUBROUTINE wrtmtd

END MODULE
