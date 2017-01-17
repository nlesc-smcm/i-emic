!#begindoc
 
MODULE m_rdfilenm

CONTAINS

SUBROUTINE rdfilenm (prompt, filenm)


CHARACTER (LEN=*)	, INTENT(IN)		:: prompt
CHARACTER (LEN=*)	, INTENT(OUT)           :: filenm

!     Writes the prompt PROMPT to standard output and reads the file
!     name from standard input in argument FILENM.  Blanks in the
!     file name are not allowed!

!     Arguments:
!     ==========
!     prompt   i   Character string written to the terminal to prompt
!                  the user to type the name of a file.
!     filenm   o   Character string containing the name of the file.
!                  The file name is left justified in FILENM and the
!                  remaining characters are blanks.

!#enddoc

!     1997-01-22  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)


INTEGER :: begfnm, endfnm, lgth

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'rdfilenm'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

lgth = LEN (filenm)

100  CONTINUE
WRITE (6, '(/, A )') prompt

READ (5, '( A )') filenm

!       Search for the start of first nonblank in filenm

begfnm = 1
110    IF ( begfnm <= lgth  .ANd. filenm(begfnm:begfnm) == ' ' ) THEN
  begfnm = begfnm + 1
  GO TO 110
END IF
IF ( begfnm > lgth )  GO TO 100

!     Determine the position of the last nonblank character of the
!     word starting in position BEGFNM in FILENM.

endfnm = begfnm
120  IF ( endfnm < lgth  .ANd. filenm(endfnm+1:endfnm+1) /= ' ' ) THEN
  endfnm = endfnm + 1
  GO TO 120
END IF

!     Shift word to beginning of FILENM and replace remaining
!     characters by (trailing) blanks.

filenm = filenm(begfnm:endfnm)

!     End of  rdfilenm
END SUBROUTINE rdfilenm

END MODULE
