!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define csrdmatrix anymatrix
#define prcmatrix  anymatrix
#define scdematrix anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_chmat

CONTAINS
 
SUBROUTINE chmat (blksiz, a, Prc)

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
USE m_wfree

INTEGER			, INTENT(IN)            :: blksiz
TYPE (csrmatrix)	, POINTER		:: a
TYPE (prcmatrix)	, POINTER		:: Prc

! replace the matrix contained in Prc structure by a
!

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




!     Allocate workspace for the matrix B  and store the reordered
!     matrix into  B:

!     Request workspace for the SCD matrix  B:

CALL wacsrd (A%n, blksiz, A%nnz, b)

!     Store the RB-reordered matrix  Prb A inv(Prb) into matrix 'B':
CALL reordrb (Prc%g, Prc%perrb, a, Ad, b)

CALL diafree(Ad)

!     Free the workspace of the matrix A:
CALL csrfree(a)

  IF (Prc%g > 0) THEN
!        Invert left upper partition of block diagonal matrix  Bd:
!           Bd(1:Prc%G,1:Prc%G) := inv(Bd(1:Prc%G,1:Prc%G))
    CALL invbdia (Prc%g, B%dia)
  END IF

! free the storage taken by the old matrix in the preconditioner
  call anyfree(Prc%aro)

!     Place reordered matrix  B  into a new matrix  Aro, [Prc%aro]:
  CALL stomat (Prc%n, Prc%g, b, Prc%aro)
! free storage of b
  CALL csrdfree(b)
 

IF (outlev >= 4) THEN
  PRINT '(A, I3, A, 2(I7, 1X, A))' , 'Step', 0, ':', Prc%g, 'of the', Prc%n, 'unknowns eliminated exactly.'
END IF


END SUBROUTINE chmat

END MODULE



