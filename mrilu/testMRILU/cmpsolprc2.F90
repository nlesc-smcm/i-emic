!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define prcmatrix  anymatrix

#endif

MODULE m_cmpsolprc2

CONTAINS

SUBROUTINE cmpsolprc2 (BlkSiz, A, b, x, NrSteps, TestX, xRef)

USE m_dump
USE m_build
USE m_wfree
USE m_cmpsol
USE m_cmpprc
USE m_prepro
USE m_solprc
USE m_applprc
USE m_prcpars
USE m_wacsr
USE m_chmat
USE m_matvec
USE m_rdmt

INTEGER					, INTENT(IN)            :: BlkSiz
TYPE (csrmatrix)			, POINTER		:: A
DOUBLE PRECISION, DIMENSION(1:A%n)	       	:: b
DOUBLE PRECISION, DIMENSION(:)		, POINTER               :: x
INTEGER					, INTENT(IN)            :: NrSteps
LOGICAL					, INTENT(IN)            :: TestX
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)         	:: xRef

!     Test for the case that we have a changing matrix while we keep the
!     preconditioning the same. 

!     Arguments:
!     ==========
!     BlkSiz   	i   Number of rows/columns in a diagonal block.
!                   The value of 'A%N' should be an integer multiple of 'BlkSiz'.
!     A      	io  In:  Location of the descriptor of the
!                       matrix  A which is stored in CSR format.
!                   Out: -- Undefined!
!                       (The storage for the matrix A has been released)
!     b       	i   In:  Right-hand side of the matrix equation Aorig x = b
!     x        	io  In:  -- Undefined! Not allocated.
!                   Out: The solution of the linear system  Aorig X = b
!     NrSteps  	i   Number of times the same computed preconditioner is
!                   used to solve the linear system.
!     TestX    	i   True if we can compare the computed solution with the 
!                   exact solution in  xRef.
!     xRef      i   If TestX: The exact solution of the linear system A xRef = b

!#enddoc

!     Local Parameters:
!     =================
CHARACTER (LEN=*), PARAMETER :: rounam = 'cmpsolprc2'

!     Local Variables:
!     ================

INTEGER						:: ier, step, i, j
TYPE (csrmatrix), POINTER                       :: Ac
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
! make a copy of A in Ac for later use
CALL wacsr (A%n, A%nnz, Ac)

    Ac%beg = A%beg
    Ac%jco = A%jco
    Ac%co  = A%co
    Ac%nnz = A%nnz
    Ac%n = A%n
    
! Compute first rhs
 b(1:A%n) = 0.0D0
 CALL matvec (A%n, 1.0D0, csrtoany(A), xRef, b)

! compute preconditioner
CALL cmpprc (BlkSiz, A, Prc)
  
!solve problem Ax=b  
  IF (TestPrec) THEN
    CALL applprc (Prc, x, b)
  ELSE
    CALL solprc (Prc, x, b)
  END IF
  ! check the accuracy of the solution
  IF (TestX) THEN

!   Test problem with exact solution in xRef.
    
!   Calculate the absolute error.
  
    abserr = MAXVAL( ABS((x(1:Prc%n)-xRef(1:Prc%n)) - (x(1)-xRef(1)) ) )

    PRINT '(/, A, 1P, E22.16)', 'Absolute error:  ', abserr

  END IF
 
 !Possibility to read new matrix
 ! READ new matrix
!  CALL rdmt ('Apert', .false.,  Ac)

! Or perturb the copy of matrix A
    DO i=1,Ac%n
      DO j=Ac%beg(i), Ac%beg(i+1)-1
        Ac%co(j)=(1+i/10)*Ac%co(j)
      ENDDO
    ENDDO    

    b(1:Ac%n) = 0.0D0
! Construct rhs for this matrix
  CALL matvec (Ac%n, 1.0D0, csrtoany(Ac), xRef, b)
! incorporate the new matrix in the existing preconditioner. 
  CALL chmat(BlkSiz,Ac,Prc)
 
! Solve the new system  
  x(1:Prc%n) = 0.0D0
  IF (TestPrec) THEN
    CALL applprc (Prc, x, b)
  ELSE
    CALL solprc (Prc, x, b)
  END IF
  
! check the accuracy of the solution
  IF (TestX) THEN

!   Test problem with exact solution in xRef.
    
!   Calculate the absolute error.
  
    abserr = MAXVAL( ABS((x(1:Prc%n)-xRef(1:Prc%n)) - (x(1)-xRef(1)) ) )

    PRINT '(/, A, 1P, E22.16)', 'Absolute error:  ', abserr

  END IF



! Free the storage of the preconditioner:

CALL prcfree(Prc)

END SUBROUTINE cmpsolprc2

END MODULE
