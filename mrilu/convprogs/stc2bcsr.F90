!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

PROGRAM stc2bcsr

USE m_build
USE m_rdfilenm
USE m_fstrlen
USE m_ioerrmsg
USE m_wrtvec
USE m_wrtmt
USE m_rdstc
USE m_wacsr

!     Conversion from stencil format to binary CSR format.

!     Reads a coefficient matrix and right hand side of a linear system
!     from an input file;  the matrix and right hand side are stored in
!     stencil format (ASCII).  The matrix is written to a binary, i.e.
!     unformatted, file in CSR format and the right hand side is written
!     to a separate binary file.
!     The names of the 3 files are read from standard input, i.e. the
!     terminal.  Blanks in file names are not allowed!
!     The maximum number of equations in the linear system allowed is
!     2 400 000.

!#enddoc

!     1997-09-17  ,  Doeke de Vries
!     2003-03-06  ,  Doeke de Vries     Increase number of equations.

!     Parameter constants:
!     ====================
!     MaxNeq    Maximum number of equations.
!     MaxNNZ    Maximum number of nonzeros.
!     fnmlgth   Maximum length of the file names.

INTEGER :: maxneq, maxnnz
PARAMETER      ( maxneq = 2400000, maxnnz = 5*maxneq )

INTEGER :: fnmlgth
PARAMETER      ( fnmlgth = 128 )

CHARACTER (LEN=1), PARAMETER      :: rounam = 'stc2bcsr' 

!     Local variables:
!     ================

CHARACTER (LEN=fnmlgth) 		:: infnm, matfnm, rhsfnm
INTEGER 				:: errinf, fnmlen, ier
TYPE (csrmatrix), POINTER		:: a
DOUBLE PRECISION, DIMENSION(1:maxneq) 	:: rhs

! SAVE             a%beg, a%jco, a%co, rhs

CALL wacsr (maxneq, maxnnz, a)

WRITE ( 6, '(/, 5X, A )' )  &
    'Conversion from stencil format to binary CSR format'

!     Read the 3 filenames, with trailing blanks stored in variables
!     from standard input.

CALL rdfilenm ('Enter name input file with stencil:', infnm)

CALL rdfilenm ('Enter name output file for matrix:', matfnm)

CALL rdfilenm ('Enter name input file for right hand size:', rhsfnm)

!     Read stencil from file INFNM

fnmlen = fstrlen(infnm)
CALL rdstc (infnm(1:fnmlen), .false., a, rhs, ier, errinf)

IF (ier /= 0) THEN
  CALL ioerrmsg (infnm(1:fnmlen), 'rdstc', rounam, ier, errinf)
  STOP
END IF

!     Write right hand side vector to file RHSFNM

fnmlen = fstrlen(rhsfnm)
CALL wrtvec (rhsfnm(1:fnmlen), .true., a%n, rhs, 2, ier)

IF (ier /= 0) THEN
  CALL ioerrmsg (rhsfnm(1:fnmlen), 'wrtvec', rounam, ier, 0)
  STOP
END IF

!     Write matrix to file MATFNM

fnmlen = fstrlen(matfnm)
CALL wrtmt (matfnm(1:fnmlen), .true., a, 2, ier)

IF (ier /= 0) THEN
  CALL ioerrmsg (matfnm(1:fnmlen), 'wrtbmt', rounam, ier, errinf)
  STOP
END IF

!     End of  stc2bcsr
END PROGRAM stc2bcsr
