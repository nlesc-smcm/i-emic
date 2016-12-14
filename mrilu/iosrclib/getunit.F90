!#begindoc
 
MODULE m_getunit

CONTAINS

SUBROUTINE getunit (unitnr, ier)

INTEGER		, INTENT(OUT)		:: unitnr
INTEGER		, INTENT(OUT)           :: ier

!     Returns a free, i.e. unopened, logical unit number in UNITNR,
!     which can be used as logical unit number in a subsequent OPEN
!     statement.

!     Arguments:
!     ==========
!     unitnr   o   A free logical unit number if IER = 0.
!                  Undefined if IER != 0.
!     ier      o   The error code:
!                  =  0   No error
!                  = -1   No free logical unit number available.
!                  >  0   A Fortran error code.

!#enddoc

!     1997-01-26  ,  Latest version of Mark Aves
!     1997-03-07  ,  Return error code and change FUNCTION to SUBROUTINE
!                    Change non-ANSI statements.


!     MACHINE DEPENDENT FUNCTION.
!     The search begins from unit number 99, the highest possible value,
!     down to  1, the lowest possible value.

! --- Local Parameters

INTEGER, PARAMETER :: start = 99

! --- Local Variables
LOGICAL :: lopen

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'getunit'

!     TRACE INFORMATION
PRINT '(A, 1X, A)' , 'Entry:', rounam
#endif

! +++ Start the search from known position
unitnr = start + 1

! +++ While the channel number is being used, search for an new channel
lopen = .true.
DO WHILE ( lopen  .ANd. unitnr >= 1)
  unitnr = unitnr - 1
  INQUIRE (UNIT=unitnr, ERR=1000, IOSTAT=ier, OPENED=lopen)
  1000 CONTINUE
END DO

IF ( unitnr < 1 ) ier = - 1

#ifdef DEBUG

PRINT '(A, I0, A)' , 'unitnr = ', unitnr
PRINT '(A, X, A)' , 'Leave:', rounam

#endif


END SUBROUTINE getunit

END MODULE





