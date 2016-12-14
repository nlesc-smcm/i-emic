!#begindoc
 
MODULE m_rdprcpars

CONTAINS
 
SUBROUTINE rdprcpars (iunit)

USE m_glbpars
USE m_iniprc
USE m_dump

INTEGER, INTENT(IN)                  :: iunit

!     Read the parameters to
!     control the preprocessing of a linear system, and to
!     control the construction of a preconditioner of the matrix of
!     this system.
!     Initialise the common blocks /prcpars/.

!     The expected format of the input file consists of the lines with:
!        Dummy line
!        Letter (T/F)       Scale rows of input matrix, so that for
!                           the scaled matrix A:
!                              A(i:: SUM(j:: ABS(A(i,j))) = 1)
!        Letter (T/F)       Apply Reverse Cuthill-McKee reordering
!                           of the matrix A.
!        Letter (T/F)       XactElm.  Use exact elimination in the
!                           original linear system to reduce the size of
!                           the system to be solved with an iterative
!                           method.

!        Double Precision   GusFctr. Apply Gustafsson's Modification's
!                           type of lumping on the block-diagonal
!                           sub-matrices, using the factor 'GusFctr'.
!                           0.0 <= GusFctr <= 1.0

!        Letter (T/F)       CLSOnce.  Compute the Lump Space once in
!                           stead of at each new Schur-complement.
!        Dummy line
!        Double Precision   NLSFctr.  New Lump Space Factor:
!                           The new lump space used
!                           >= NLSFctr * Lump space from new Schur-compl.
!                           Only used when  .NOT. CLSOnce.
!                           0 <= NLSFctr <= 1
!        Double Precision   EpsW.  Lump space tolerance:
!                           Lumpspace = Epsw / MAX(ABS(inv(diag(SC))))
!                           0 <= EpsW
!        Double Precision   ElmFctr, element factor times the free Lump
!                           Space gives maximum value of an element that
!                           may be lumped.
!        Double Precision   RedFctr, reduction factor for Lump Space
!                           from the new Schur-complement.
!                           Free Lump Space
!                           = RedFctr * Lump Space Schur-complement.
!                           0 <= RedFctr <= 1.0D0

!        Double Precision   SchTol, Schur Tolerance: the minimum
!                           value of elements in the newly computed
!                           Schur-Complement.  0 <= SchTol

!        Double Precision   DensLim, density limit for the last
!                           partitioned block.
!        Double Precision   GlobFrac, the global fraction limit for the
!                           last partitioned block.
!        Double Precision   LocFrac, local fraction limit for the last
!                           partitioned block.
!        Double Precision   SparsLim, limit for a Sparse last block.
!                           If density last block
!                           <= SparsLim, use some ILU-factorisation,
!                           else, use LU-factorisation of last block.

!        Integer            ILUType.  Type of ILU factorisation of sparse
!                           last block of Multi Level Preconditioner:
!                           0 = (M)ILU0
!                           1 = (M)ILU1  Level of fill = 1
!                           2 = (M)ILU2  Level of fill = 2
!                           ...
!                           9 = (M)ILUT
!        Double Precision   DropTol.  Drop Tolerance in ILUT
!                           factorization of last block.  0<=DropTol.
!        Double Precision   CompFctr.  Compensation Factor for the
!                       f    diagonal in the ILUT factorization of the
!                           last block if  DropTol>0.
!                           0<=CompFctr<=1.
!        Dummy line
!        Double Precision   CPivTol.  Change Pivot Tolerance for a
!                           non-diagonal pivot element in the
!                           Incomplete LDU factorization of the last
!                           block.
!                           ABS(U(i,j)) * CPivTol > ABS(U(i,i))
!        Double Precision   LUTol, LU Tolerance; last block is
!                           singular if
!                              MIN(|diag(U)|)/MAX(|diag(U)|)
!                              <
!                              LUTol * MAX(i::ABS(U(i,i)))
!        Dummy line
!        Letter (T/F)       SingLU.  Allow singular U in LU
!                           factorization of last block.
!                           Only in last diagonal element of U!

!     Argument:
!     =========
!     iunit    i   Unit number of the open input file containing the
!                  input parameters.

!#enddoc


!     Global Parameters:
!     ==================

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'rdprcpars'
#endif

!     Local Variables:
!     ================
CHARACTER (LEN=80) :: line

!     Preprocessing parameters:
LOGICAL :: scarow, usercm, xactelm, TestPrec

!     Preconditioner parameters:
DOUBLE PRECISION :: compfctr, cpivtol, denslim, gusfctr, globfrac, locfrac,  &
    sparslim, droptol, elmfctr, epsw, lutol, nlsfctr, redfctr, schtol
INTEGER :: ilutype
LOGICAL :: clsonce, singlu

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


!     Skip the first line:
READ (iunit, '(A80)') line


!     Preprocessing parameters:

!     Scale A such that 1-norm of A is unity
READ (iunit, *) scarow

!     Use Reverse Cuthill-McKee reordering
READ (iunit, *) usercm

!     Use exact elimination to reduce system size:
READ (iunit, *) xactelm


!     Preconditioner parameters:

!     Gustafsson's modification, Lump space tolerance, factors and
!     Last block density

READ (iunit, *) gusfctr

READ (iunit, *) clsonce
READ (iunit, '(A80)') line
READ (iunit, *) nlsfctr

IF (.NOT. clsonce .AND. (nlsfctr < 0.0D0 .OR. nlsfctr > 1.0D0)) & 
   CALL dump(__FILE__,__LINE__, 'Value NLSFctr NOT in [0.0D0,1.0D0]!')

READ (iunit, *) epsw
IF (epsw < 0.0D0 .OR. epsw > 10.0D0) CALL dump(__FILE__,__LINE__, 'Value EpsW NOT in [0.0D0,1.0D1]!')

READ (iunit, *) elmfctr
IF (elmfctr < 0.0D0 .OR. elmfctr > 1.0D0) CALL dump(__FILE__,__LINE__, 'Value ElmFctr NOT in [0.0D0,1.0D0]!')

READ (iunit, *) redfctr
IF (redfctr < 0.0D0 .OR. redfctr > 1.0D2) CALL dump(__FILE__,__LINE__, 'Value RedFctr NOT in [0.0D0,1.0D2]!')

READ (iunit, *) schtol
IF (schtol < 0.0D0 .OR. schtol > 1.0D-1) CALL dump(__FILE__,__LINE__, 'Value SchTol NOT in [0.0D0,1.0D-1]!')

READ (iunit, *) denslim
IF (denslim < 0.0D0 .OR. denslim > 1.0D0) CALL dump(__FILE__,__LINE__, 'Value DensLim NOT in [0.0D0,1.0D0]!')

READ (iunit, *) globfrac
IF (globfrac < 0.0D0 .OR. globfrac > 1.0D0) CALL dump(__FILE__,__LINE__, 'Value GlobFrac NOT in [0.0D0,1.0D0]!')

READ (iunit, *) locfrac
IF (locfrac < 0.0D0 .OR. locfrac > 1.0D0) CALL dump(__FILE__,__LINE__, 'Value LocFrac NOT in [0.0D0,1.0D0]!')

READ (iunit, *) sparslim
IF (sparslim < 0.0D0 .OR. sparslim > 1.0D0) CALL dump(__FILE__,__LINE__, 'Value SparsLim NOT in [0.0D0,1.0D0]!')

READ (iunit, *) ilutype
IF (ilutype < 0  .Or. ilutype > 9) CALL dump(__FILE__,__LINE__, 'Value ILUType NOT in [0:1:9]!')

READ (iunit, '(A80)') line

READ (iunit, *) droptol
IF (ilutype == 9 .AND. droptol < 0.0D0) CALL dump(__FILE__,__LINE__, 'Negative value DropTol NOT allowed!')

READ (iunit, *) compfctr
IF (compfctr < 0.0D0 .OR. compfctr > 1.0D0) CALL dump(__FILE__,__LINE__, 'Value CompFctr NOT in [0.0D0,1.0D0]!')

READ (iunit, '(A80)') line

READ (iunit, *) cpivtol
IF (ilutype == 9 .ANd. (cpivtol < 0.0D0 .OR. cpivtol > 1.0D0)) CALL dump(__FILE__,__LINE__, 'Value CPivTol NOT in [0.0D0,1.0D0]!')

READ (iunit, *) lutol
READ (iunit, '(A80)') line
READ (iunit, *) singlu
READ (iunit, *) TestPrec

IF (singlu .AND. lutol < 0.0D0) CALL dump(__FILE__,__LINE__, 'Negative value LUTol NOT allowed!')

IF (outlev >= 1) THEN
  PRINT '(A)', ' '
  IF (scarow) THEN
    PRINT '(A)', 'Scale rowsums of Matrix |A| to 1.'
  END IF
  IF (usercm) THEN
    PRINT '(A)', 'Use Reverse Cuthill-McKee reordering.'
  END IF
  IF (xactelm) THEN
    PRINT '(A)', 'Use Exact elimination.'
  END IF
  IF ( gusfctr /= 1.0D0 ) THEN
    PRINT '(A)', 'Use Gustafsson in Multi Level Preconditioner.'
  END IF
  

  PRINT '(1P, A, 1X, E8.2)', 'GusFctr =', gusfctr
  
  PRINT '(A, 1X, L2, A, 2X, A, 1X, F5.3, A)',  &
      'CLSOnce =', clsonce, ',', 'NLSFctr =', nlsfctr, ','
  
  PRINT '(3(A, 1X, F5.3, A, 2X), 1P, A, 1X, E8.2, A)', 'EpsW =', epsw, ',',  &
      'ElmFctr =', elmfctr, ',', 'RedFctr =', redfctr, ',',  &
      'SchTol =', schtol, ','
  
  PRINT '(A, 1X, F5.3, A, 3(2X, A, 1X, F5.3, A))',  &
      'DensLim =', denslim, ',', 'GlobFrac =', globfrac, ',',  &
      'LocFrac =', locfrac, ',', 'SparsLim =', sparslim, ','
  
  IF (ilutype == 9) THEN
    PRINT '(/, A)', 'ILUT factorisation of last sparse block.'
    PRINT '(1P, A, 1X, E8.2, 2(2X, A, 2X, A, 1X, E8.2))',  &
        'CompFctr =', compfctr, ',', 'DropTol =', droptol, ',',  &
        'CPivTol =', cpivtol
  ELSE
    PRINT '(/, A, I1, A)', 'ILU(', ilutype,  &
        ') factorisation of last sparse block.'
    PRINT '(1P, A, 1X, E8.2)', 'CompFctr =', compfctr
  END IF
  
  PRINT '(1P, A, 1X, E8.2, A, 2X, A, L2)', 'LUTol =', lutol, ',',  &
      'SingLU =', singlu

!  IF (TestPrec) THEN
!    PRINT '(A)', 'Calls applprc instead of solprc.'
!  END IF

END IF


CALL iniprc (usercm, scarow, xactelm, clsonce, nlsfctr, epsw, elmfctr, gusfctr, redfctr, schtol,  &
    denslim, globfrac, locfrac, sparslim, ilutype, droptol, compfctr, cpivtol, lutol, singlu, testprec)

END SUBROUTINE rdprcpars

END MODULE
