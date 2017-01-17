!#begindoc
 
MODULE m_xfnminr

CONTAINS

SUBROUTINE xfnminr (fnm, ndd, inr, resfnm)

USE m_fstrlen

CHARACTER (LEN=*)	, INTENT(IN)		:: fnm
INTEGER			, INTENT(IN)            :: ndd
INTEGER			, INTENT(IN)        	:: inr
CHARACTER (LEN=*)	, INTENT(OUT)           :: resfnm

!     Extend the file name 'fnm' by appending the last 'ndd' decimal
!     digits of the non-negative integer number 'inr', and place the
!     result in 'resfnm'.
!     The program stops with an error message when the actual argument
!     'resfnm' is too short.

!     Arguments:
!     ==========
!     fnm      i   File name string
!     ndd      i   Number of decimal digits with leading zeros,
!                  1 <= ndd <= 4.
!     inr      i   Integer number, inr >= 0
!     resfnm   o   File name string =
!                  fnm // decimal representation of MOD(ABS(inr),10**ndd)

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'xfnminr'

!     Local variables:
!     ================
CHARACTER (LEN=8) 	:: vfmt
INTEGER 		:: fnmlen, nd, nr

nd = MIN(4,MAX(1,ndd))
nr = ABS(inr)
WRITE ( vfmt, '( ''(A,I'', I1, ''.'', I1, '')'' )' ) nd, nd

fnmlen = fstrlen(fnm)
IF (LEN(resfnm) < fnmlen+nd) THEN
  PRINT '(/, A, 2X, A, 2(/, 3X, A, I4))', 'Fatal error occurred in', rounam,  &
      'Length CHARACTER argument of result is:', LEN(resfnm), 'Length should be at least:             ', fnmlen+nd
  STOP 'in xfnminr'
END IF

WRITE (resfnm, vfmt)  fnm(1:fnmlen), MOD(nr,10**nd)

!     End of  xfnminr
END SUBROUTINE xfnminr

END MODULE