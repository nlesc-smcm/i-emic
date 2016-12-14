!#begindoc
 
#ifndef WITH_UNION

#define partmatrix  anymatrix

#endif

MODULE m_wrtbldu

CONTAINS

SUBROUTINE wrtbldu (basenm, Part, ier)

USE m_getunit
USE m_fstrlen
USE m_build

CHARACTER (LEN=*)			, INTENT(IN)	:: basenm
TYPE (partmatrix)			, POINTER	:: Part
INTEGER					, INTENT(OUT)	:: ier

!     Writes the sparse matrices  L, D and U  in CSR format into the
!     3 binary files
!        '<basenm>.L', '<basename>.D'  and  '<basenm>.U'

!     The off-diagonal elements of the matrices L and U are stored in
!     CSR format in [Part%offd%beg,Part%offd%jco,Part%offd%co]; all diagonal elements of L and
!     U have the value 1.
!     The inverse of the diagonal matrix D is stored in the vector
!     [part%dia%com].
!     The elements of L are stored in
!        [PART%N, Part%offd%beg, Part%offd%jco, Part%offd%co, Part%lnzl]
!     and are written into  '<basenm>.L'.
!     The elements of U are stored in
!        [PART%N, Part%offd%beg, Part%offd%jco, Part%offd%co, Part%lnzl]
!     and are written into  '<basenm>.U'.

!     File structure
!     --------------
!     The '.L', '.D' and '.U' files will each contain 4 records, with:
!       record 1   1   * integer:          order of matrix [ = PART%N ]
!       record 2   PART%N+1 * integer:          indices indicating beginning
!                                          of a row of  L, D  or  U
!       record 3   NNZ * integer:          column numbers of nonzero
!                                          elements of  L, D  or  U.
!                                          N.B. Not necessarily in
!                                          increasing order!!!
!       record 4   NNZ * double precision: nonzero values of the
!                                          matrices  L, D  or  U

!     Arguments:
!     ==========
!     basenm   i   Base of the output filenames
!     PART%N        i   The order of the matrices  L, D and U
!     Part%offd%beg    i   Part%offd%beg(i) = index in array Part%offd%jco or Part%offd%co for the first
!                  nonzero off-diagonal element (L+U)(i,j) in row i
!                  (1 <= i <= PART%N+1)
!                  Part%offd%beg(PART%N+1) = 1 + index of last nonzero off-diagonal
!                  element of (L+U)
!     Part%offd%jco    i   Part%offd%jco(nz) = column number of nz-th nonzero
!                  off-diagonal element of L+U
!     Part%offd%co     i   Part%offd%co(nz)  = value of nz-th nonzero off-diagonal
!                  element of L+U  (1 <= nz <= Part%offd%beg(PART%N+1) - 1
!     Part%lnzl     i   Part%lnzl(i) = index in array Part%offd%jco or Part%offd%co of the last
!                  nonzero off-diagonal element of L in row i
!                  (1 <= i <= PART%N)
!     part%dia%com    i   part%dia%com(1,i)  = value of i-th diagonal element  inv(D)
!                  (1 <= i <= PART%N)
!     ier      o   The error code,
!                  =  0   No error
!                  = -3   Error during closing, opening or rewinding
!                         file <basenm>.D , <basenm>.L  or  <basenm>.U
!                  = -5   Error during writing <basenm>.D , <basenm>.L
!                         or  <basenm>.U
!                  =-12   No logical unit available to write on any of
!                         the files  <basenm>.D , <basenm>.L  or
!                         <basenm>.U

!#enddoc

!     1998-08-26  ,  Doeke D. de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!-----------------------------------------------------------------------

INTEGER, DIMENSION(1:Part%N+1)		:: beg

!     Local parameter:
!     ================

INTEGER, PARAMETER :: maxfnmlen = 80

!     Local variables:
!     ================
INTEGER :: baslen, i, iunit, nz

CHARACTER (LEN=maxfnmlen) :: filenm

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'wrtbldu'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

baslen =  MIN (fstrlen(basenm), maxfnmlen-4)

!     Get an free logical unit number for writing

CALL getunit (iunit, ier)
IF ( ier /= 0 )  GO TO 988


!     Write the data to the '.D' file:

filenm(1:baslen+2) = basenm(1:baslen) // '.D'

OPEN  (iunit, FILE = filenm(1:baslen+2), STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', ERR = 997)
REWIND (iunit, ERR = 997)

WRITE (iunit, ERR = 995) Part%N

WRITE (iunit, ERR = 995)  (i, i = 1, Part%N+1)

WRITE (iunit, ERR = 995)  (i, i = 1, Part%N)

WRITE (iunit, ERR = 995) ( 1.0D0 / part%dia%com(1,i), i = 1, Part%N )

CLOSE (iunit)


!     Write the data to the '.L' file:

beg(1) = 1
DO i = 1, Part%N
  beg(i+1) = beg(i) + (Part%lnzl(i) - Part%offd%beg(i) + 2)
END DO

filenm(1:baslen+2) = basenm(1:baslen) // '.L'

OPEN  (iunit, FILE = filenm(1:baslen+2), STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', ERR = 997)
REWIND (iunit, ERR = 997)

WRITE (iunit, ERR = 995)  Part%N

WRITE (iunit, ERR = 995)  (beg(i), i = 1, Part%N+1)

WRITE (iunit, ERR = 995)  ( ( Part%offd%jco(nz), nz = Part%offd%beg(i), Part%lnzl(i) ), i, i = 1, Part%N )

WRITE (iunit, ERR = 995)  ( ( Part%offd%co(nz), nz = Part%offd%beg(i), Part%lnzl(i) ), 1.0D0, i = 1, Part%N )

CLOSE (iunit, ERR = 997)


!     Write the data to the '.U' file:

beg(1) = 1
DO i = 1, Part%N
  beg(i+1) = beg(i) + (Part%offd%beg(i+1) - Part%lnzl(i))
END DO

filenm(1:baslen+2) = basenm(1:baslen) // '.U'

OPEN  (iunit, FILE = filenm(1:baslen+2), STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', ERR = 997)
REWIND (iunit, ERR = 997)

WRITE (iunit, ERR = 995)  Part%N

WRITE (iunit, ERR = 995)  (beg(i), i = 1, Part%N+1)

WRITE (iunit, ERR = 995) ( i, ( Part%offd%jco(nz), nz = Part%lnzl(i)+1, Part%offd%beg(i+1)-1 ), i = 1, Part%N )

WRITE (iunit, ERR = 995) ( 1.0D0, ( Part%offd%co(nz), nz = Part%lnzl(i)+1, Part%offd%beg(i+1)-1 ), i = 1, Part%N )

CLOSE (iunit, ERR = 997)


!     Normal return:
RETURN


!     Error return:
988  ier = ier - 7
995  ier = ier - 2
997  ier = ier - 3

!     End of  wrtbldu
END SUBROUTINE wrtbldu

END MODULE