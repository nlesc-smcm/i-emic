!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define csrdmatrix anymatrix
#define prcmatrix  anymatrix
#define scdematrix anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_prepro

CONTAINS
 
SUBROUTINE prepro (blksiz, a, Prc, S)

USE m_build
USE m_wacsr
USE m_wacsrd
USE m_wamlp
USE m_wfree
USE m_copymt
USE m_crbper
USE m_glbpars
USE m_prcpars
USE m_eblkdia
USE m_invbdia
USE m_reordrb
USE m_scalmat
USE m_schurcmpl
USE m_stomat
USE m_wadia

INTEGER			, INTENT(IN)            :: blksiz
TYPE (csrmatrix)	, POINTER		:: a
TYPE (prcmatrix)	, POINTER		:: Prc
TYPE (scdematrix)	, POINTER		:: S

!     Preprocess the matrix A:
!     . Update the fields and segments in [Prc]:
!       .. Compute the scale factors in 'ScaFctr' and scale the
!          matrix A,  Asc := diag(ScaFctr) A(:,:).
!       .. Compute the permutationvector in 'PerRB' of the Red-Black
!          permutation, Prb, and reorder the matrix Asc accordingly.
!          Aro := Prb Asc Prb'.
!       .. Compute the order in 'G' of the left upper sub_matrix Aro_11
!          of this reordered matrix  Aro.
!          Aro_11 is a block-diagonal matrix.
!          (G = 0 if NOT XactElm, in /prcpars/.)
!       .. Compute the Schur-complement  SC  of  Aro_11  in  Aro.
!     . Release the storage of matrix A and
!       store the reordered matrix Aro.

!     Arguments:
!     ==========
!     N        	i   Number of rows/columns in the matrix  A.
!     Prc%G           Number of rows/columns in the sub_matrix  Aro_11  in
!                   the reordered matrix  Aro.
!     BlkSiz   	i   Number of rows/columns in the blocks of the block-
!                   diagonal of matrix.
!                   'N' is an integral multiple of 'BlkSiz'.
!     a      	io  In:  Location of the descriptor of the
!                       matrix of the linear system stored in CSR format.
!                   Out: -- Undefined!
!     Prc    	io  Location of descriptor for the
!                   preconditioner of the original matrix A.
!     S      	o   Schur-complement of A11 of the partitioned Red-Black
!                   reordered matrix  A, possibly with the small elements.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'prepro'

!     Local Variables:
!     ================


TYPE (csrdmatrix), POINTER			:: b
TYPE (diamatrix), POINTER  			:: Ad

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Calculate scaling factors and scale matrix  A  accordingly,
!        A := diag([Prc%scale]) A:
CALL scalmat (a, Prc%scale)

! Request the workspace for a copy of the Block Diagonal matrix:
CALL wadia (A%n, blksiz, Ad)

!     Extract Block Diagonal from A, updating [begA,jcoA,coA], and store
!     the block diagonal into [Ad%co]:
CALL eblkdia (a, Ad)

!     Compute the permutationvector for a Red-Black reordering in
!     [Prc%perrb], and
!     compute the size of the left upper block diagonal sub-matrix of
!     the reordered matrix  B  in 'G':
CALL crbper (a, Ad, Prc)

! begin added 1-7-2005 according to old fortran files

!IF (Prc%g .EQ. Prc%N) THEN
!
!  Reordered matrix is a non-singular block-diagonal matrix.
!  Do NOT use exact elimination.

!   Prc%g = 0
!   Prc%nschur = Prc%N

! ENDIF

! end added

!     Allocate workspace for the matrix B  and store the reordered
!     matrix into  B:

!     Request workspace for the SCD matrix  B:

CALL wacsrd (A%n, blksiz, csrnnz(A), b)
!PRINT *,Prc%g
!     Store the RB-reordered matrix  Prb A inv(Prb) into matrix 'B':
CALL reordrb (Prc%g, Prc%perrb, a, Ad, b)
!PRINT *,Prc%g
CALL diafree(Ad)

!     Free the workspace of the matrix A:
CALL csrfree(a)

  IF (Prc%g > 0) THEN
!        Invert left upper partition of block diagonal matrix  Bd:
!           Bd(1:Prc%G,1:Prc%G) := inv(Bd(1:Prc%G,1:Prc%G))
    CALL invbdia (Prc%g, B%dia)
  END IF
  
!     Place reordered matrix  B  into a new matrix  Aro, [Prc%aro]:
  CALL stomat (Prc%n, Prc%g, b, Prc%aro)

!     Allocate a descriptor for a new MLP type preconditioner and a segment 
!     for the permutation vector:
  CALL wamlp (Prc%nschur, Prc%mlp)
  
  IF (Prc%g > 0) THEN
!        Free the workspace for  B:
    CALL csrdfree(b)
    
!        Construct the Schur-complement of Aro_11 in the reordered matrix
!        'Aro', into 'S', dropping the small elements.
!        The matrix  S  will have (Prc%NSCHUR) rows/columns:
    CALL schurcmpl (Prc%g, Prc%aro, schtol, S)
  ELSE
!        Copy Schur-complement [s] into [s],
!        dropping the small elements:
    CALL copymt (csrdtoscde(b), schtol, S)
    
!        Free the workspace for  B:
    CALL csrdfree(b)
  END IF

IF (outlev >= 4) THEN
  PRINT '(A, I3, A, 2(I7, 1X, A))' , 'Step', 0, ':', Prc%g, 'of the', Prc%n, 'unknowns eliminated exactly.'
END IF


END SUBROUTINE prepro

END MODULE



