!#begindoc

MODULE m_guerrmsg

CONTAINS

SUBROUTINE guerrmsg (callnm, ier)

CHARACTER (LEN=*)	, INTENT(IN)            :: callnm
INTEGER			, INTENT(IN)            :: ier


!     Prints an error message corresponding with the value of the error
!     code in IER (!= 0), returned from a call to getunit.
!     See also the documentation file getunit.txt!

!     Arguments:
!     ==========
!     callnm   i   The name of the calling program unit.
!     ier      i   The error code (see also getunit):
!                  =  0   No error
!                  = -1   No free logical unit number available.
!                  >  0   A Fortran error code.

!#enddoc

!     1997-05-12  ,  Doeke de Vries.
!     2003-03-05  ,  Last update (Doeke de Vries).

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'guerrmsg'

!     TRACE INFORMATION
WRITE (6,'(A)') rounam
#endif

IF ( ier /= 0 ) THEN
  PRINT '(/, A, A)' , 'Fatal error occurred in  getunit  called from  ', callnm
  IF ( ier == -1 ) THEN
    PRINT '(A, /)' , 'Too many open units in use!'
  ELSE IF ( ier > 0 ) THEN
    PRINT '(A, /, 3X, A, I11, /)' , 'Fortran I/O error!',  'Error number: ', ier
  ELSE
    PRINT '(A, /, A, I11, /)' , 'Internal error!',  'Illegal error code: ', ier
    STOP
  END IF
END IF

!     End of  guerrmsg
END SUBROUTINE guerrmsg

END MODULE