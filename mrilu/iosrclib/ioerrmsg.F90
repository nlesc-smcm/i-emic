!#begindoc
 
MODULE m_ioerrmsg

CONTAINS

SUBROUTINE ioerrmsg (filenm, subrnm, callnm, ier, errinf)

USE m_fstrlen

CHARACTER (LEN=*)	, INTENT(IN)		:: filenm
CHARACTER (LEN=*)	, INTENT(IN)            :: subrnm
CHARACTER (LEN=*)	, INTENT(IN)            :: callnm
INTEGER			, INTENT(IN)            :: ier
INTEGER			, INTENT(IN)            :: errinf

!     Prints an error message corresponding with the value of the error
!     code in IER (!= 0), returned from a call of one of the subroutines
!     from the module "iosrclib".
!     See also the documentation files for the subroutines rd....txt and
!     wrt....txt.

!     Arguments:
!     ==========
!     filenm   i   The name of the file.
!     subrnm   i   The name of the subroutine that returned the error
!                  code in IER.
!     callnm   i   The name of the calling program unit, i.e. the unit
!                  that containes the call of the subroutine SUBRNM.
!     ier      i   The error code:
!                  =  0   No error, no error message printed.
!                  <  0   An fatal error code.
!                  >  0   A nonfatal error code.
!     errinf   i   Extra error information, interpretation dependent on
!                  the value of 'ier'.

!#enddoc

!     1997-07-19  ,  Doeke de Vries.
!     2003-03-05  ,  Last update (Doeke de Vries)

LOGICAL :: isvsm
INTEGER :: LEN

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'ioerrmsg'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

LEN = fstrlen(subrnm)
isvsm = .false.
IF (LEN > 3) THEN
  IF (subrnm(LEN-2:LEN) == 'vsm') THEN
    isvsm = .true.
  END IF
END IF

LEN = fstrlen(filenm)

IF ( ier /= 0 ) THEN
  PRINT '( /, A, 2X, A, 2X, A, 2X, A )' , 'Fatal error occurred in', subrnm, 'called from', callnm
  IF ( ier == -1 ) THEN
    PRINT '(A, /, 3X, A, A, A, X, I11, /)' , 'Too many rows/columns in matrix/vector on file!',  &
        'Number of rows/columns in ', filenm(1:LEN), ' is:', errinf
  ELSE IF ( ier == - 2 ) THEN
    PRINT '(A, /, 3X, A, A, A, X, I11, /)' , 'Too many nonzero elements in matrix on file!',  &
        'Number of nonzeros in ', filenm(1:LEN), ' is:', errinf
  ELSE IF ( ier == - 3 ) THEN
    IF (isvsm) THEN
      PRINT '(A, 2X, A, A, /)' , 'Error during opening, rewinding or closing file', filenm(1:LEN), '.beg!'
    ELSE
      PRINT '(A, 2X, A, /)' , 'Error during opening, rewinding or closing file', filenm
    END IF
  ELSE IF ( ier == - 4 ) THEN
    IF (isvsm) THEN
      PRINT '(A, 2X, A, A, /)' , 'Unexpected End Of File encountered when reading', filenm(1:LEN), '.beg!'
    ELSE
      PRINT '(A, 2X, A, /)' , 'Unexpected End Of File encountered when reading', filenm
    END IF
  ELSE IF ( ier == - 5 ) THEN
    IF (isvsm) THEN
      PRINT '(A, 2X, A, A, /)' , 'Error during IO-operation on file', filenm(1:LEN), '.beg!'
    ELSE
      PRINT '(A, 2X, A, /)' , 'Error during IO-operation on file', filenm
    END IF
  ELSE IF ( ier == - 6 ) THEN
    PRINT '(A, 2X, A, A, /)' , 'Error during opening, rewinding or closing file', filenm(1:LEN), '.jco!'
  ELSE IF ( ier == - 7 ) THEN
    PRINT '(A, 2X, A, A, /)' , 'Unexpected End Of File encountered when reading', filenm(1:LEN), '.jco!'
  ELSE IF ( ier == - 8 ) THEN
    PRINT '(A, 2X, A, A, /)' , 'Error during IO-operation on file', filenm(1:LEN), '.jco!'
  ELSE IF ( ier == - 9 ) THEN
    PRINT '(A, 2X, A, A, /)' , 'Error during opening, rewinding or closing file', filenm(1:LEN), '.co!'
  ELSE IF ( ier == -10 ) THEN
    PRINT '(A, 2X, A, A, /)' , 'Unexpected End Of File encountered when reading', filenm(1:LEN), '.co!'
  ELSE IF ( ier == -11 ) THEN
    PRINT '(A, 2X, A, A, /)' , 'Error during IO-operation on file', filenm(1:LEN), '.co!'
  ELSE IF ( ier == -12 ) THEN
    PRINT '(A, /)' , 'No logical unit available!'
  ELSE IF ( ier == -13 ) THEN
    PRINT '(A, 2X, A, 2X, A, /)' , 'Contents of file', filenm, 'is not consistent!'
  ELSE
    PRINT '(A, /, 3X, A, X, I11, /)' , 'Internal error!', 'Illegal error code:', ier
  END IF
ELSE IF ( ier > 0 ) THEN
  PRINT '(A, /, 3X, A, X, I11, /)' , 'Internal error!', 'Illegal error code:', ier
  STOP
END IF

!     End of  ioerrmsg
END SUBROUTINE ioerrmsg

END MODULE