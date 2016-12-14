!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define prcmatrix  anymatrix

#endif

MODULE m_cmpsolprc

CONTAINS

SUBROUTINE cmpsolprc (BlkSiz, A, b, x, NrSteps, TestX, xRef)

USE m_dump
USE m_build
USE m_wfree
USE m_cmpsol
USE m_cmpprc
USE m_prepro
USE m_solprc
USE m_applprc
USE m_prcpars

INTEGER					, INTENT(IN)            :: BlkSiz
TYPE (csrmatrix)			, POINTER		:: A
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)       	:: b
DOUBLE PRECISION, DIMENSION(:)		, POINTER               :: x
INTEGER					, INTENT(IN)            :: NrSteps
LOGICAL					, INTENT(IN)            :: TestX
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)         	:: xRef

!     Test the computation of the MRILU preconditioner for the left hand
!     side of a linear system, with matrix A, together with the test of
!     the computation of the solution of this linear system,  A x = b,
!     using this MRILU preconditioner.

!     Arguments:
!     ==========
!     BlkSiz   	i   Number of rows/columns in a diagonal block.
!                   The value of 'A%N' should be an integer multiple of 'BlkSiz'.
!     A      	io  In:  Location of the descriptor of the
!                       matrix  A which is stored in CSR format.
!                   Out: -- Undefined!
!                       (The storage for the matrix A has been released)
!     b       	i   In:  Right-hand side of the matrix equation Aorig x = b
!     X        	io  In:  -- Undefined! Not allocated.
!                   Out: The solution of the linear system  Aorig X = b
!     NrSteps  	i   Number of times the same computed preconditioner is
!                   used to solve the linear system.
!     TestX    	i   To compare the computed solution with the exact solution in  xRef.
!     xRef      	i   If TestX: The exact solution of the linear system A xRef = b

!#enddoc

!     Local Parameters:
!     =================
CHARACTER (LEN=*), PARAMETER :: rounam = 'cmpsolprc'

!     Local Variables:
!     ================

INTEGER						:: ier, step
TYPE (prcmatrix), POINTER			:: Prc
DOUBLE PRECISION 				:: abserr

!-----------------------------------------------------------------------

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

! Allocation and zero the initial solution:

ALLOCATE ( x(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
x(1:A%n) = 0.0D0

! Construct complete preconditioner and store location in 'Prc':

CALL cmpprc (BlkSiz, A, Prc)

! Loop to test the REUSE of the PRECONDITIONER

DO step = 1, NrSteps

! Solve the linear system
!   A x = b
! using the complete preconditioner [Prc]:
  
  IF (TestPrec) THEN
    CALL applprc (Prc, x, b)
  ELSE
    CALL solprc (Prc, x, b)
  END IF
  
  IF (TestX) THEN

!   Test problem with exact solution in xRef.
    
!   Calculate the absolute error.
  
    abserr = MAXVAL( ABS((x(1:Prc%n)-xRef(1:Prc%n)) - (x(1)-xRef(1)) ) )

    PRINT '(/, A, 1P, E22.16)', 'Absolute error:  ', abserr

  END IF

END DO

! Free the storage of the preconditioner:

CALL prcfree(Prc)

END SUBROUTINE cmpsolprc

END MODULE