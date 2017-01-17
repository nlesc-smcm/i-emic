!#begindoc
 
MODULE m_prpars

CONTAINS

SUBROUTINE prpars

USE m_prcpars
USE m_solpars
USE m_vispars

!     Print the value of each parameter on the standard output file.

!#enddoc

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'prpars'

!     TRACE INFORMATION
      PRINT '(A, X, A)' , 'Entry:', RouNam

#endif


!     Parameters preprocessor construction:
!     -------------------------------------
 
PRINT '(/, A, /, A)'    , 'Precontioner Parameters:'          , '------------------------'
PRINT '(A, X, L8)'      , 'Apply Reverse Cuthill-McKee:      ', CutMck
PRINT '(A, X, L8)'      , 'Scale rowsums matrix to 1:        ', ScaRow
PRINT '(A, X, L8)'      , 'Use exact elimination:            ', XactElm
PRINT '(A, X, L8)'      , 'Compute Lump space once:          ', CLSOnce
PRINT '(A, X, 1P, E8.2)', 'New Lump space factor:            ', NLSFctr
PRINT '(A, X, 1P, E8.2)', 'Lump space tolerance:             ', EpsW
PRINT '(A, X, 1P, E8.2)', 'Element factor:                   ', ElmFctr
PRINT '(A, X, 1P, E8.2)', 'Gustafsson''s factor:             ', GusFctr
PRINT '(A, X, 1P, E8.2)', 'Reduction factor:                 ', RedFctr
PRINT '(A, X, 1P, E8.2)', 'Schur tolerance:                  ', SchTol
PRINT '(A, X, 1P, E8.2)', 'Density limit last block:         ', DensLim
PRINT '(A, X, 1P, E8.2)', 'Global fraction limit:            ', GlobFrac
PRINT '(A, X, 1P, E8.2)', 'Local fraction limit:             ', LocFrac
PRINT '(A, X, 1P, E8.2)', 'Sparse limit last block:          ', SparsLim
PRINT '(A, X, I8)'      , 'Type ILU factorisation last block:', ILUType
IF (ILUType .EQ. 9) THEN
  PRINT '(A, X, 1P, E8.2)','Drop tolerance ILUT factorisation:', DropTol
  PRINT '(A, X, 1P, E8.2)','Pivot tolerance in ILUT fact.:    ', CPivTol
ENDIF


!     Parameters solver:
!     ------------------

PRINT '(/, A, /, A)'      ,'Solver Parameters:'                , '------------------'
PRINT '(A, X, I8)'        ,'Type Congruent Gradient solver:   ', CGType
IF (CGType .EQ. 3) THEN
  PRINT '(A, X, I8)'      ,'Dimension subspace in GMRES:      ', Mgmres
ENDIF
PRINT '(A, X, I8)'        ,'Maximum number iterations:        ', MaxNIts
PRINT '(A, X, 1P, E8.2)'  ,'Reduction tolerance:              ', RedTol
PRINT '(A, X, 1P, E8.2)'  ,'Absolute tolerance:               ', AbsTol


!     Visualisation parameters:
!     -------------------------

PRINT '(/, A, /, A)'      ,'Visualisation Parameters:'         , '-------------------------'
PRINT '(A, X, L8)'        ,'Original (scaled) matrix:         ', VisAsc
PRINT '(A, X, I8)'        ,'Number reoordered matrices:       ', VisNRO
PRINT '(A, X, I8)'        ,'Number Schur-complement matrices: ', VisNSC
PRINT '(A, X, L8)'        ,'Last Schur-complement matrix:     ', VisLSC
PRINT '(A, X, L8)'        ,'Incomplete LDU factorisation:     ', VisILDU

PRINT '(A)', ' '

END SUBROUTINE

END MODULE
