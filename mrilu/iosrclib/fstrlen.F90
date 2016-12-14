!#begindoc
 
MODULE m_fstrlen

CONTAINS

INTEGER FUNCTION fstrlen (str)

CHARACTER (LEN=*)	, INTENT(IN)	:: str

!     Returns the length of the Fortran string, i.e. the character
!     value, in STR without trailing blanks; this length equals the
!     position of the last nonblank character in STR.
!     If STR contains only blanks the value 1 is returned.

!     Argument:
!     =========
!     str     i   the Fortran string (charater value) which "length" is
!                 returned.

!#enddoc

INTEGER :: lgth

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'fstrlen'

PRINT '(A, X, A)' , 'Entry:', rounam
PRINT '(A, X, A)' , 'str = ', str

#endif

lgth = LEN(str)
100  IF ( lgth > 1  .ANd. str(lgth:lgth) == ' ' ) THEN
  lgth = lgth - 1
  GO TO 100
END IF

fstrlen = lgth

!     End of  fstrlen
END FUNCTION fstrlen

END MODULE






