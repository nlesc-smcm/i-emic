!#begindoc

MODULE m_mterrmsg 

CONTAINS

SUBROUTINE mterrmsg (callnm, atp, etp)

USE m_build

CHARACTER (LEN=*), INTENT(IN)            :: callnm
INTEGER, INTENT(IN)                      :: atp
INTEGER, INTENT(IN)                      :: etp

!     Matrix Type error message.
!     Prints an error message on the standard output file, if the actual
!     matrix/partition type differs from the expected type or that the
!     type used is not implemented or unknown.

!     Arguments:
!     ==========
!     callnm   i   The name of the program unit where the mismatch of
!                  storage type was detected.
!     atp      i   The actual matrix storage type.
!     etp      i   The expected matrix storage type.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'mterrmsg'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


PRINT '(/, A, 2X, A, /, A)' , 'Fatal error detected in', callnm,  &
    'Internal error!'

IF ( atp /= etp ) THEN
  PRINT '( A, /, 2(3X, A, I4, /) )' , 'Mismatch of matrix/partition type:',  &
      'The expected type is:', etp, 'the actual type is:  ', atp
ELSE IF ( etp == bdiartp  .Or. etp == jdstp  .Or.  &
      etp == symtp  ) THEN
  PRINT '(A, /, 3X, A, I4)' , 'Matrix type not implemented!',  &
      'Matrix type is:', etp
ELSE
  PRINT '(A, /, 3X, A, I11)' , 'Illegal matrix/partition type!',  &
      'Matrix/partition type is:', etp
END IF

!     End of  mterrmsg
END SUBROUTINE mterrmsg

END