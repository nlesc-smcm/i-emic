!#begindoc
 
MODULE m_rdpblpars

CONTAINS

SUBROUTINE rdpblpars (iunit, readeq, blksiz, nrintv, probty, strtch, cvx, cvy, testx)

USE m_glbpars
USE m_dump


INTEGER			, INTENT(IN)		:: iunit
INTEGER			, INTENT(OUT)           :: readeq
INTEGER			, INTENT(OUT)           :: blksiz
INTEGER			, INTENT(OUT)           :: nrintv
INTEGER			, INTENT(OUT)           :: probty
DOUBLE PRECISION	, INTENT(OUT)         	:: strtch
DOUBLE PRECISION	, INTENT(OUT)         	:: cvx
DOUBLE PRECISION	, INTENT(OUT)         	:: cvy
LOGICAL			, INTENT(OUT)           :: testx

!     Read the parameters which define the problem from the open file
!     indicated by unit number 'iunit'.

!     The expected format of the input file consists of the lines with:
!        A dummy line
!        Integer            Type number ReadEq: Read or Generate the
!                           linear system.
!                           ReadEq = 1  Convection-Diffusion
!                           ReadEq = 2  Convection-Diffusion-3D
!                           ReadEq = 3  Read ASCII file (Stencil format)
!                           ReadEq = 4  Read binary file (CSR format)
!                           ReadEq = 5  Read ASCII file (CSR format)
!        A dummy line
!        Integer            The size of the diagonal blocks in the
!                           matrix of the linear system.
!        Integer (1,2,...)  NrIntv, the number of intervals in  X- and
!                           Y-direction.
!                           Only used if ReadEq = 1,2
!        Integer            ProbTy, the problem type:
!                           ProbTy = 1,4 DDDD
!                           ProbTy = 2,5 DDNN
!                           ProbTy = 3,6 NNNN
!                           Only used if ReadEq = 1,2
!        Double Precision   Strtch, grid stretching:
!                           0.0001 <= Strtch <= 10000
!                           Only used if ReadEq = 1,2
!        Double Precision   CVX, Convection term: coefficient of  u_x
!                           Only used if ReadEq = 1,2
!        Double Precision   CVY, Convection term: coefficient of  u_y
!                           Only used if ReadEq = 1,2
!        Letter (T/F)       Construct right-hand side for the solution
!                           x = (1,2,...)

!     Argument:
!     =========
!     iunit    i   Unit number of the open input file containing the
!                  input parameters.
!     ReadEq   o   Read or generate equation
!     BlkSiz   o   Size of diagonal blocks in matrix linear system.
!     NrIntv   o   Number of intervals in each of  X- and Y-direction.
!     ProbTy   o   Problem type
!     Strtch   o   Grid stretching
!     CVX      o   Coefficient of  u_x
!     CVY      o   Coefficient of  u_y
!     TestX    o   Construct  rhs  for a given solution.

!#enddoc


!     Global Parameters:
!     ==================

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'rdpblpars'
#endif

!     Local Variables:
!     ================
CHARACTER (LEN=80) :: line

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


!     Skip dummy line:
READ (iunit, '(A80)') line

!     Read or set up system of equations
READ (iunit, *) readeq

IF (readeq < 1  .Or. readeq > 5) CALL dump(__FILE__,__LINE__, 'Code to read or generate is NOT known!')

!     Skip dummy line:
READ (iunit, '(A80)') line

!     Read size of diagonal blocks in matrix linear system:
READ (iunit, *) blksiz
IF (blksiz < 0) CALL dump(__FILE__,__LINE__, 'Size of diagonal blocks is NOT correct!')

!     Number intervals, if equations generated
READ (iunit, *) nrintv
IF ((readeq == 1  .Or. readeq == 2)  .AND. (nrintv < 1))  CALL dump(__FILE__,__LINE__, 'Illegal number of intervals')

READ (iunit, *) probty
IF ((readeq == 1  .Or. readeq == 2) .AND. (probty < 1  .Or. probty > 6)) & 
 CALL dump(__FILE__,__LINE__, 'Problem Type: is NOT known!')

READ (iunit, *) strtch
IF ((readeq == 1  .Or. readeq == 2) .AND. (strtch < 1.0D-4 .OR. 1.0D4 < strtch)) &
  CALL dump(__FILE__,__LINE__, 'Value grid stretching NOT in [1.0D-4,1.0D+4]!')

READ (iunit, *) cvx
READ (iunit, *) cvy

!     Construct  rhs  to give solution X=(1,2,...)
READ (iunit,*) testx


IF (outlev >= 1) THEN
  PRINT '(A)', ' '
  
!        Print problem type:
  IF (readeq == 1) THEN
    PRINT '(A)', 'Construct Convec-Diff Eqn'
  ELSE IF (readeq == 2) THEN
    PRINT '(A)', 'Construct 3D Convec-Diff Eqn'
  ELSE IF (readeq == 3) THEN
    PRINT '(A)', 'Eqn from ASCII files (Stencil)'
  ELSE IF (readeq == 4) THEN
    PRINT '(A)', 'Eqn from binary files (CSR)'
  ELSE IF (readeq == 5) THEN
    PRINT '(A)', 'Eqn from ASCII files (CSR)'
  END IF
  
  IF (readeq == 1  .Or. readeq == 2) THEN
    PRINT '(A, 1X, I5)',  'Nr Intvals =', nrintv
    PRINT '(1P, A, 1X, E8.2, 2(A, 2X, A, 1X, E8.2))',  &
        'Stretching =', strtch, ',', 'Coeff U_x =', cvx, ',', 'Coeff U_y =', cvy
  END IF
  
  IF (testx) THEN
    PRINT '(A)', 'Construct RHS for given solution.'
  END IF
END IF


END SUBROUTINE rdpblpars

END MODULE