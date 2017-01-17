!#begindoc
 
PROGRAM avczn2avc

USE m_rdfilenm
USE m_fstrlen
USE m_ioerrmsg
USE m_rdvc
USE m_wrtvec

!     Adds the length of an ASCII vector in front of it.

!     Reads a vector from an input file.  The vector is written to a
!     binary, unformatted file.
!     The names of the 2 files are read from standard input, i.e. the
!     terminal.  Blanks in file names are not allowed!
!     The maximum number of elements allowed is  10 000 000.

!#enddoc

!     1997-09-17  ,  Doeke de Vries
!     2003-03-06  ,  Doeke de Vries     Increase max. vector length.

!     Parameter constants:
!     ====================
!     MaxN      Maximum number of elements.
!     fnmlgth   Maximum length of the file names.



INTEGER :: fnmlgth
PARAMETER      ( fnmlgth = 128)

CHARACTER (LEN=*), PARAMETER      :: rounam = 'avcnz2avc' 

!     Local variables:
!     ================

CHARACTER (LEN=fnmlgth) 		:: infnm, outfnm
INTEGER 				:: errinf, ier, fnmlen, nelm
DOUBLE PRECISION, DIMENSION(:), POINTER	:: v

! SAVE             v

WRITE ( 6, '(/, 5X, A )' ) 'Puts the number of elements in front of ASCII vector'
!     Read the 2 filenames, with trailing blanks, stored in variables
!     from standard input.

CALL rdfilenm ('Enter name ASCII input file:', infnm)

CALL rdfilenm ('Enter name ASCII output file:', outfnm)

!     Read vector from file INFNM:

fnmlen = fstrlen(infnm)
CALL rdvc (infnm(1:fnmlen), .false., nelm, v, ier)

IF (ier /= 0) THEN
  CALL ioerrmsg (infnm(1:fnmlen), 'rdvc', rounam, ier, errinf)
  STOP
END IF

!     Write vector to file OUTFNM:

fnmlen = fstrlen(outfnm)
CALL wrtvec (outfnm(1:fnmlen), .false. , nelm, v, 16, ier)

IF ( ier /= 0 ) THEN
  CALL ioerrmsg (outfnm(1:fnmlen), 'wrtvec', rounam, ier, 0)
  STOP
END IF

!     End of  avczn2avc
END PROGRAM avczn2avc
