!#begindoc
 
MODULE m_rdsolpars

CONTAINS
 
SUBROUTINE rdsolpars (iunit)

USE m_glbpars
USE m_inisol
USE m_dump

INTEGER, INTENT(IN)                  :: iunit

!     Read the parameters to control the solution the linear system of
!     equations with CG-like method, and initialise common block
!     /solpars/.

!     The expected format of the input file consists of the lines with:
!        A dummy line
!        Integer Number        CGType, the Type number of a Congruent
!                              Gradient type solver:
!                              CGType = 0  Use CG
!                              CGType = 1  Use Bi-CGSTAB using a Left
!                                          preconditioner
!                              CGType = 2  Use Bi-CGSTAB using a Right
!                                          preconditioner
!                              CGType = 3  Use GMRES
!        A dummy line
!        Integer Number        M, the dimension of the subspace used in
!                              GMRES
!        Integer Number        MaxNIts: Maximum number of iterations in
!                              the CG-like solver
!        Double Precision      RedTol: Factor indicating the required
!                              reduction of the residual of the
!                              preconditioned system.
!        Double Precision      AbsTol: Upper bound residual of the
!                              preconditioned system.
!     The iteration process in the CG-type solver stops as soon as one
!     the 3 conditions is satisfied:
!        Number of iterations = MaxNIts           .OR.
!        Reduction factor of residuals <= RedTol  .OR.
!        Residual <= AbsTol

!     Argument:
!     =========
!     iunit    i   Unit number of the open input file containing the
!                  input parameters.

!#enddoc


!     Global Parameters:
!     ==================

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'rdsolpars'
#endif

CHARACTER (LEN=80) :: line

INTEGER :: cgtype, maxnits, mgmres
DOUBLE PRECISION :: abstol, redtol

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Skip the first line:
READ (iunit, '(A80)') line

!     Solution method and convergence criteria
READ (iunit, *) cgtype
READ (iunit, '(A80)') line

IF (cgtype < 0  .OR. cgtype > 3) CALL dump(__FILE__,__LINE__, 'Value CG type NOT in [0,1,2,3]!')

READ (iunit, *) mgmres
IF (cgtype == 3 .AND. mgmres <= 0.0D0) CALL dump(__FILE__,__LINE__, 'Value Mgmres should be > 0!')

READ (iunit, *) maxnits
IF (maxnits <= 0.0D0) CALL dump(__FILE__,__LINE__, 'Value MaxNIts should be > 0!')

READ (iunit, *) redtol
IF (redtol < 0.0D0) CALL dump(__FILE__,__LINE__,  'Negative value RedTol NOT allowed!')

READ (iunit, *) abstol
IF (abstol < 0.0D0) CALL dump(__FILE__,__LINE__, 'Negative value VisNSC NOT allowed!')


IF (outlev >= 1) THEN
  PRINT '(A)', ' '
  
!        Print CG type:
  IF (cgtype == 0) THEN
    PRINT '(A)', 'Method: CG'
  ELSE IF (cgtype == 1) THEN
    PRINT '(A)', 'Method: Bi-CGSTAB'
  ELSE IF (cgtype == 2) THEN
    PRINT '(A)', 'Method: Bi-CGSTABr'
  ELSE IF (cgtype == 3) THEN
    PRINT '(A, I3)', 'Method: GMRES(m), with m =', mgmres
  END IF
  
  PRINT '(1P, A, 1X, I4, 2(A, 2X, A, 1X, E8.2))',  &
      'Max Nr Iters =', maxnits, ',', 'Residual <=', abstol, ',',  &
      'Decrease residual <=', redtol
END IF

CALL inisol (cgtype, mgmres, maxnits, redtol, abstol)

END SUBROUTINE rdsolpars

END MODULE