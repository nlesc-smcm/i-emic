MODULE m_prcpars

!-----------------------------------------------------------------------

!     Preconditioner parameters.

!    

!     The parameters in the COMMON BLOCK /prcpars/ get a default value
!     during declaration. The values of these
!     parameters can be changed by the SUBROUTINE 'iniprc'.


!     Description of the parameters in /prcpars/:
!     ===========================================
!     CutMcK    Apply Reverse Cuthill-McKee ordering of the original
!               matrix.
!     ScaRow    Scale rows of input matrix, so that for the scaled
!               matrix A, A_sc:
!               .TRUE.:  Row sums of absolute values of elements are
!                        equal 1:
!                        A(i :: SUM(j :: ABS(A_sc(i,j))) = 1)
!               .FALSE.: Maximum of abolute values of the elements
!                        of  A_sc  is in the range [1,2)  and
!                        all diagonal elements are non-negative:
!                        1 <= MAX(i,j :: ABS(A_sc(i,j))) < 2  AND
!                        A(i :: A_sc(i,i) >= 0)
!     XactElm   Use exact elimination in the original linear system to
!               reduce the size of the system to be solved with an
!               iterative method.

!     CLSOnce   Compute Lump Space Once:
!               .TRUE.:  Compute Lump Space for the 1st Schur-
!                        complement only.
!               .FALSE.: Compute Lump Space for each newly computed
!                        Schur-complement.
!     NLSFctr   New Lump Space Factor:
!               The new lump space used
!               >= NLSFctr * Lump space from new Schur-complement.
!               Only used when  .NOT. CLSOnce.
!               0 <= NLSFctr <= 1
!     EpsW      Drop tolerance for the lumping strategy;
!               Lumpspace = Epsw / MAX(ABS(inv(diag(Schur-complement))))
!               EpsW >= 0.
!     ElmFctr   Element Factor times the free Lump Space gives maximum
!               value of an element that may be lumped.
!               0 <= ElmFctr <= 1.0D0
!     GusFctr   Apply Gustafsson modification's type of lumping on the
!               block-diagonal sub-matrices, using the factor 'GusFctr'.
!     RedFctr   Reduction factor for the Lump Space from the Schur-
!               complement:
!               Free Lump Space = RedFctr * Lump Space Schur-complement.
!               0 <= RedFctr <= 1.0D0

!     SchTol    Schur Tolerance: the non-zero off-diagonal elements of
!               the newly computed Schur-complement to be stored should
!               be greater than 'SchTol'.

!               The decision whether or not the last partitioned block
!               has been reached is based on the values of the variables
!               DensLim, GlobFrac and LocFrac. The elimination of blocks
!               is stopped when one of the 3 conditions is satisfied.
!     DensLim   Density limit for the last partitioned block.
!               The density of the last constructed Schur-complement
!               should be at least 'DensLim' in order to become the
!               last block of the Multi Level Preconditioner.
!     GlobFrac  Global fraction limit for the last partitioned block.
!               The ratio of the order of the last constructed
!               Schur-complement and the order of the original matrix
!               should be less than 'GlobFrac' in order to become the
!               last block of the Multi Level Preconditioner.
!     LocFrac   Local fraction limit for the last partitioned block.
!               The ratio of the order of the last eliminated (left
!               upper) block and the order of the last constructed
!               Schur-complement should be less than 'LocFrac' in order
!               to become the last block of the Multi Level
!               Preconditioner.


!     SparsLim  Limit for a Sparse last block.
!               If the density of the last block is at most 'SparsLim'
!               the representation of the last block remains sparse and
!               an Incomplete LDU factorization is made of this block.
!               If the density of the last block is greater than
!               'SparsLim' the representation of the last block is
!               changed to a full matrix and a LU factorization is made.

!     ILUType   ILUType.  Type of ILU factorisation of sparse last block
!               of Multi Level Preconditioner:
!               0 = (M)ILU0
!               1 = (M)ILU1  Level of fill = 1
!               2 = (M)ILU2  Level of fill = 2
!               ...
!               9 = (M)ILUT
!     DropTol   Drop Tolerance for the ILUT factorization of the
!               sparse last block.  Setting  DropTol = 0  produces the
!               complete LDU factorization.
!               0<=DropTol.

!     CompFctr  Compensation Factor for the diagonal in the Incomplete
!               LDU factorization of the last block.  The sum of the
!               discarded (dropped) elements in a row of L and DU,
!               multiplied by 'CompFctr', are added to the diagonal
!               element of D.
!               0<=CompFctr<=1.
!               Used only if  DropTol > 0
!     CPivTol   Change Pivot Tolerance for a non-diagonal pivot element
!               in the Incomplete LDU factorization of the last block.
!               Two columns 'i' and 'j' in row 'i' of the matrix DU are
!               permuted when:
!                  ABS(DU(i,j)) * CPivTol > ABS(DU(i,i))

!     LUTol     Tolerance to determine the singularity of factor U,
!               in LU-factorisation of last block.
!                  "U is singular" <==>
!                  MIN(i::ABS(U(i,i))) <= LUTol * MAX(i::ABS(U(i,i)))
!     SingLU    Singular DU factor allowed in the complete LDU
!               factorization of last block. Only in last diagonal
!               element of  DU!
!     TestPrec  Calls applprc instead of solprc


DOUBLE PRECISION ::		&
	compfctr=1.00D0		,&
	cpivtol	=8.75D-1	,&
	denslim	=1.00D-2	,& 
	droptol	=1.00D-3	,&
	epsw	=1.00D-1	,&
	elmfctr	=2.00D-1	,& 
	globfrac=0.0		,&
	gusfctr	=1.00D0		,&
	locfrac	=0.0		,&
	lutol	=1.0D-10	,&
	nlsfctr	=1.00D-1	,&
	redfctr	=8.00D-1	,&
	schtol	=0.0D0		,&
	sparslim=6.66D-1

INTEGER :: 			&
	ilutype	=9

LOGICAL ::			&
	clsonce	=.false.	,&
	cutmck	=.false.	,& 
	scarow	=.false.	,& 
	singlu	=.false.	,& 
	xactelm	=.true.         ,&
        testprec=.false.


!-----------------------------------------------------------------------

END MODULE