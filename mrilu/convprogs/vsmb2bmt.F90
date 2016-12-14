!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

PROGRAM vsm2bmt

USE m_build
USE m_dump
USE m_rdfilenm
USE m_ioerrmsg
USE m_wrtmt
USE m_rdvsmb
USE m_fstrlen
USE m_wacsr

!     Conversion from an ASCII CSR format matrix to a binary,
!     unformatted, CSR format matrix

!     Reads a coefficient matrix from an input file;  the matrix is
!     stored in ASCII CSR format.  The matrix is written to a
!     binary, unformatted file, also in CSR format.
!     The names of the 2 files are read from standard input, i.e. the
!     terminal.  Blanks in file names are not allowed!
!     The maximum number of equations in the linear system allowed is
!     1 800 000.
!     The maximum number of non-zeros in the sparse matrix is
!     18 000 000

!#enddoc

!     1997-09-17  ,  Doeke de Vries
!     2003-03-06  ,  Doeke de Vries  Increase max. number of equations
!                                    and increase number of non-zeros.

!     Parameter constants:
!     ====================

!     fnmlgth   Maximum length of the file names.

INTEGER :: fnmlgth
PARAMETER      ( fnmlgth = 128)

CHARACTER (LEN=*), PARAMETER      :: rounam = 'vsm2bmt'

!     Local variables:
!     ================

CHARACTER (LEN=fnmlgth) 		:: infnm, outfnm
INTEGER 				:: errinf, ier, fnmlen
TYPE (csrmatrix), POINTER		:: a

! SAVE             a%beg, a%jco, a%co

WRITE ( 6, '(/, 5X, A )' ) 'Conversion from ASCII CSR matrix to Binary CSR matrix'

!     Read the 2 filenames, with trailing blanks, stored in variables
!     from standard input.

CALL rdfilenm ('Enter name ASCII input file:', infnm)

CALL rdfilenm ('Enter name binary output file:', outfnm)

!     Read matrix from file INFNM:

fnmlen = fstrlen(infnm)
CALL rdvsmb (infnm(1:fnmlen), a, ier, errinf)

IF (ier /= 0) THEN
  CALL ioerrmsg (infnm(1:fnmlen), 'rdvsm', rounam, ier, errinf)
  STOP
END IF

!     Write matrix to file OUTFNM:

fnmlen = fstrlen(outfnm)
CALL wrtmt (outfnm(1:fnmlen), .true., a, 16, ier)

IF ( ier /= 0 ) THEN
  CALL ioerrmsg (outfnm(1:fnmlen), 'wrtbmt', rounam, ier, 0)
  STOP
END IF

!     End of  vsm2bmt
END PROGRAM vsm2bmt
