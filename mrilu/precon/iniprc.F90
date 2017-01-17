!#begindoc
 
MODULE m_iniprc

CONTAINS

SUBROUTINE iniprc (acutmck, ascarow, axactelm,  &
    aclsonce, anlsfctr, aepsw, aelmfctr, agusfctr, aredfctr, aschtol,  &
    adenlim, aglobfrac, alocfrac, asparlim,  &
    ailutyp, adroptol, acompfct, acpivtol, alutol, asinglu, aTestPrec)

USE m_prcpars

LOGICAL, INTENT(IN)                      :: acutmck
LOGICAL, INTENT(IN)                      :: ascarow
LOGICAL, INTENT(IN)                      :: axactelm
LOGICAL, INTENT(IN)                      :: aclsonce
DOUBLE PRECISION, INTENT(IN)             :: anlsfctr
DOUBLE PRECISION, INTENT(IN)             :: aepsw
DOUBLE PRECISION, INTENT(IN)             :: aelmfctr
DOUBLE PRECISION, INTENT(IN)             :: agusfctr
DOUBLE PRECISION, INTENT(IN)             :: aredfctr
DOUBLE PRECISION, INTENT(IN)             :: aschtol
DOUBLE PRECISION, INTENT(IN)             :: adenlim
DOUBLE PRECISION, INTENT(IN)             :: aglobfrac
DOUBLE PRECISION, INTENT(IN)             :: alocfrac
DOUBLE PRECISION, INTENT(IN)             :: asparlim
INTEGER, INTENT(IN)                      :: ailutyp
DOUBLE PRECISION, INTENT(IN)             :: adroptol
DOUBLE PRECISION, INTENT(IN)             :: acompfct
DOUBLE PRECISION, INTENT(IN)             :: acpivtol
DOUBLE PRECISION, INTENT(IN)             :: alutol
LOGICAL, INTENT(IN)                      :: asinglu
LOGICAL, INTENT(IN)                      :: aTestPrec








!     INItialise parameters for the construction of the
!     PReConditioner.

!     Initialises the variables in prcpars.
!     These variables are the parameters which control the computation of
!     the preconditioner.

!     See also the description in the file  'prcpars.F90'

!     Arguments:
!     ==========
!     aCutMcK  i   Apply Reverse Cuthill-McKee ordering of the original
!                  matrix.
!     aScaRow  i   Scale rows of input matrix, so that for the scaled
!                  matrix A, A_sc:
!                  .TRUE.:  Row sums of absolute values of elements are
!                           equal 1:
!                           A(i :: SUM(j :: ABS(A_sc(i,j))) = 1)
!                  .FALSE.: Maximum of abolute values of the elements
!                           of  A_sc  is in the range [1,2)  and
!                           all diagonal elements are non-negative:
!                           1 <= MAX(i,j :: ABS(A_sc(i,j))) < 2  AND
!                           A(i :: A_sc(i,i) >= 0)
!     aXactElm i   Use exact elimination in the original linear system
!                  to reduce the size of the system to be solved with an
!                  iterative method.
!     aCLSOnce i   Compute Lump Sum Once
!                  .TRUE.:  Compute Lump Space for the 1st Schur-
!                           complement only.
!                  .FALSE.: Compute Lump Space for each new computed
!                           Schur-complement.
!     aNLSFctr i   New Lump Space Factor:
!                  The new lump space used
!                   >= NLSFctr * Lump space from new Schur-complement.
!                  Only used when  .NOT. CLSOnce.
!                  0 <= NLSFctr <= 1
!     aEpsW    i   Drop tolerance for the lumping strategy:
!                  Lumpspace = Epsw / MAX(ABS(inv(diag(Schur-compl.))))
!                  0 <= EpsW
!     aElmFctr i   Element Factor times the free Lump Space gives
!                  maximum value of an element that may be lumped.
!                  0 <= aElmFctr <= 1.0D0
!     aGusfctr  i   Apply Gustafsson modification's type of lumping on
!                  the (block-)diagonal submatrices.
!     aRedFctr i   Reduction Factor for the Lump Space from the Schur-
!                  complement:
!                  Free Lump Space = RedFctr * Lump Space Schur-compl.
!                  0 <= aRedFctr <= 1.0D0
!     aSchTol  i   Schur Tolerance: the non-zero off-diagonal elements
!                  of the newly computed Schur-complement to be stored
!                  should be greater than 'SchTol'.
!                  0 <= aSchTol
!     aDenLim  i   Block density limit for the last partitioned block.
!                  The density of the last constructed Schur-complement
!                  should be less than 'DensLim' in order to become the
!                  last block of the Multi Level Preconditioner.
!     GlobFrac  i  Global fraction limit for the last partitioned block.
!                  The ratio of the order of the last constructed
!                  Schur-complement and the order of the original matrix
!                  should be less than 'GlobFrac' in order to become the
!                  last block of the Multi Level Preconditioner.
!     LocFrac   i  Local fraction limit for the last partitioned block.
!                  The ratio of the order of the last eliminated (left
!                  upper) block and the order of the last constructed
!                  Schur-complement should be less than 'LocFrac' in
!                  order to become the last block of the Multi Level
!                  Preconditioner.
!     aSparLim i   Limit for a Sparse last block.
!                  If the density of the last block is at most
!                  'SparsLim' the representation of the last block
!                  and an Incomplete LDU factorization is made of this
!                  block.  If the density of the last block is greater
!                  than 'SparsLim' the representation of the last block
!                  is changed to a full matrix and a LU factorization is
!                  made.
!     aILUTyp  i   Type of ILU factorisation of sparse last block
!                  of Multi Level Preconditioner:
!                  0 = (M)ILU0
!                  1 = (M)ILU1  Level of fill = 1
!                  2 = (M)ILU2  Level of fill = 2
!                  ...
!                  9 = (M)ILUT
!     aDropTol i   Drop Tolerance for the Incomplete LDU factorization
!                  of the last block.  Setting  aDropTol = 0  produces
!                  the complete LDU factorization.
!                  0<=DropTol.
!     aCompFct i   Compensation Factor for the diagonal in the Incomplete
!                  LDU factorization of the last block.  The sum of the
!                  discarded (dropped) elements in a row of L and DU,
!                  multiplied by 'CompFctr', are added to the diagonal
!                  element of D.
!                  0<=CompFctr<=1.
!                  Used only if  DropTol > 0
!     aCPivTol i   Change Pivot Tolerance for a non-diagonal pivot
!                  element in the Incomplete LDU factorization of the
!                  last block.
!                  Two columns 'i' and 'j' in row 'i' of the matrix DU
!                  are permuted when:
!                     ABS(DU(i,j)) * aCPivTol > ABS(DU(i,i))
!     aLUTol   i   Tolerance to determine the singularity of factor U,
!                  in LU-factorisation of last block.
!                     "U is singular" <==>
!                     MIN(i::ABS(U(i,i))) <= LUTol * MAX(i::ABS(U(i,i)))
!     aSingLU  i   Singular U factor allowed in LU-factorisation of last
!                  block. Only in last diagonal element of U!
!     aTestPrec i  Calls applprc instead of solprc
!#enddoc

!     Global variables:
!     =================

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'iniprc'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

cutmck  = acutmck
scarow  = ascarow
xactelm = axactelm

clsonce = aclsonce
nlsfctr = anlsfctr
epsw    = aepsw
elmfctr = aelmfctr
gusfctr = agusfctr
redfctr = aredfctr

schtol  = aschtol

denslim  = adenlim
globfrac = aglobfrac
locfrac  = alocfrac
sparslim = asparlim

ilutype  = ailutyp
droptol  = adroptol

compfctr = acompfct
cpivtol  = acpivtol

lutol    = alutol
singlu   = asinglu
TestPrec = aTestPrec

!     End of  iniprc
END SUBROUTINE iniprc


END MODULE