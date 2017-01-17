!#begindoc
 
PROGRAM bvc2avc

USE m_rdfilenm
USE m_fstrlen
USE m_rdvec
USE m_ioerrmsg
USE m_wrtvec
USE m_getunit

!     Conversion from binary, unformatted, vector to an ASCII vector.

!     Reads a vector from a binary, i.e. unformatted, input file.
!     The vector is written to an ASCII file.
!     The names of the 2 files are read from standard input, i.e. the
!     terminal.  Blanks in file names are not allowed!
!     The maximum number of elements allowed is  10 000 000.

!#enddoc

!     1997-09-18  ,  Doeke de Vries
!     2003-03-06  ,  Doeke de Vries     Increase max. vector length.

!     External function:
!     ==================
!INTEGER :: fstrlen
!EXTERNAL fstrlen

!     Parameter constants:
!     ====================
!     MaxN      Maximum number of elements.
!     fnmlgth   Maximum length of the file names.

INTEGER :: maxn
PARAMETER      ( maxn = 10000000 )

INTEGER :: fnmlgth
PARAMETER      ( fnmlgth = 128)

CHARACTER (LEN=1), PARAMETER      :: rounam = 'bvc2avc'

!     Local variables:
!     ================

CHARACTER (LEN=1) :: infnm*(fnmlgth), outfnm*(fnmlgth)

INTEGER 				:: errinf, ier, fnmlen, nelm
DOUBLE PRECISION, DIMENSION(1:maxn) 	:: v

! SAVE             v


WRITE ( 6, '(/, 5X, A )' ) 'Conversion from Binary vector to ASCII vector'

!     Read the 2 filenames, with trailing blanks, stored in variables
!     from standard input.

CALL rdfilenm ('Enter name Binary input file:', infnm)

CALL rdfilenm ('Enter name ASCII output file:', outfnm)

!     Read vector from file INFNM:

fnmlen = fstrlen(infnm)
CALL rdvec (infnm(1:fnmlen), .true., nelm, v, ier, errinf)

IF (ier /= 0) THEN
  CALL ioerrmsg (infnm(1:fnmlen), 'rdvec', rounam, ier, errinf)
  STOP
END IF

!     Write vector to file OUTFNM:

fnmlen = fstrlen(outfnm)
CALL wrtvec (outfnm(1:fnmlen), .false.,  nelm, v, 16, ier)

IF ( ier /= 0 ) THEN
  CALL ioerrmsg (outfnm(1:fnmlen), 'wrtvec', rounam, ier, 0)
  STOP
END IF

!     End of  bvc2avc
END PROGRAM bvc2avc

