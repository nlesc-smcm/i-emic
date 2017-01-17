!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

PROGRAM bmt2amt

USE m_build
USE m_rdfilenm
USE m_fstrlen
USE m_rdmt
USE m_wrtmt
USE m_ioerrmsg
USE m_wacsr

!     Conversion from a binary, i.e. unformatted, matrix to an
!     ASCII matrix.

!     Reads a coefficient matrix from an input file;  the matrix is
!     stored in binary, i.e. Fortan unformatted, CSR format on file.
!     The matrix is written to an ASCII file, also in CSR format.
!     The names of the 2 files are read from standard input, i.e. the
!     terminal.  Blanks in file names are not allowed!
!     The maximum number of equations in the linear system allowed is
!     1 800 000.
!     The maximum number of non-zeros in the sparse matrix is
!     18 000 000

!#enddoc

!     1997-09-18  ,  Doeke de Vries
!     2003-03-06  ,  Doeke de Vries  Increase max. number of equations
!                                    and increase number of non-zeros.

!     Parameter constants:
!     ====================
!     fnmlgth   Maximum length of the file names.

INTEGER :: fnmlgth
PARAMETER      ( fnmlgth = 128)

CHARACTER (LEN=1), PARAMETER      :: rounam = 'bmt2amt'

!     Local variables:
!     ================

CHARACTER (LEN=fnmlgth) 		:: infnm, outfnm
INTEGER 				:: errinf, ier, fnmlen, neq
TYPE (csrmatrix), POINTER		:: a

! SAVE             a%beg, a%jco, a%co


WRITE ( 6, '(/, 5X, A )' ) 'Conversion from Binary CSR matrix to ASCII CSR matrix'

!     Read the 2 filenames, with trailing blanks, stored in variables
!     from standard input.

CALL rdfilenm ('Enter name Binary input file:', infnm)

CALL rdfilenm ('Enter name ASCII output file:', outfnm)

!     Read matrix from file INFNM:

fnmlen = fstrlen(infnm)
CALL rdmt (infnm(1:fnmlen), .true., a)

!     Write matrix to file OUTFNM:

fnmlen = fstrlen(outfnm)
CALL wrtmt (outfnm(1:fnmlen), .false., a, 16, ier)

IF ( ier /= 0 ) THEN
  CALL ioerrmsg (outfnm(1:fnmlen), 'wrtmt', rounam, ier, 0)
  STOP
END IF

!     End of  bmt2amt
END PROGRAM bmt2amt
