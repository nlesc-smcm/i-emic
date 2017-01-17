!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define prcmatrix  anymatrix
#define scdematrix anymatrix

#endif

MODULE m_cmpprc

CONTAINS

SUBROUTINE cmpprc (blksiz, A, Prc)

USE m_build
USE m_waprc
USE m_bepnum
USE m_glbpars
USE m_prepro
USE m_prpars
USE m_dump

INTEGER			, INTENT(IN)    	:: blksiz
TYPE (csrmatrix)	, POINTER		:: A
TYPE (prcmatrix)	, POINTER		:: Prc

!     Computes the Preconditioner for a A%NxA%N matrix A, in MLP format,
!     and stores it.
!     This preconditioner can be referenced through the descriptor
!     at location 'Prc'.
!     The storage of the original matrix A, referenced through 'A',
!     is released.
!     The scaled and partitioned Red-Black reordered matrix A, is stored
!     in a new area and can be reached via the
!     descriptor of the preconditioner [Prc].

!     Arguments:
!     ==========
!     A%N        i   Number of rows/columns in the matrix  A.
!     BlkSiz   i   Number of rows/columns in a diagonal block.
!                  The value of 'A%N' should be an integer multiple of
!                  'BlkSiz'.
!     A      io  In:  Location of the descriptor of the
!                       matrix of the linear system stored in CSR format.
!                  Out: -- Undefined!
!     Prc    o   Location of descriptor for the
!                  preconditioner of matrix A.

!     It is assumed that the global variables in the common blocks have
!     been initialised by calling the subroutines 'iniprc' and 'inivis'.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'cmpprc'

!     Local Variables:
!     ================

TYPE (scdematrix), POINTER		:: S
DOUBLE PRECISION 			:: begtim, cumtim, endtim

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


!     Check the input parameters:
!     ---------------------------

IF (blksiz <= 0) THEN
  PRINT '(/, A, 2X, A, /, 3X, I9, 2X, A, /)' ,  &
      'Fatal error occurred in', rounam,  &
      blksiz, 'rows/columns in diagonal blocks not allowed!'
  CALL dump(__FILE__,__LINE__,'Fatal error occurred, rows/columns in diagonal blocks not allowed!')
END IF

IF (MOD(A%n,blksiz) /= 0) THEN
  PRINT '(/, A, 2X, A, /, 3X, A, /)' , 'Fatal error occurred in', rounam,  &
      'Block size not compatible with order of matrix!'
  CALL dump(__FILE__,__LINE__,'Fatal error occurred, block size not compatible with order of matrix!')
END IF

IF (OutLev .GE. 5) THEN
   CALL prpars
ENDIF

!     Compute the preconditioner:
!     ---------------------------

!     Begin timing:
CALL CPU_TIME(begtim)

!     Allocate storage for a descriptor and the segments directly
!     referenced from this descriptor:
CALL waprc (A%n, Prc)

!     Preprocess matrix  A  by row and column permutations to give the
!     matrix  Aro  and construct in case of Exact Elimination the block
!     partitioning of  Aro. The storage of the matrix  A  is released.
!     The (partitioned) reordered matrix  Aro is stored
!     and can be referenced  via the descriptor of the preconditioner
!     [Prc].
!     The Schur-complement of the block  Aro_11  in  Aro  is stored 
!     and can be referenced through 'S'.
!     If no Exact Elimination is requested  S  equals  A_ro, possibly
!     with the small elements dropped.
!     The fields and segments in [Prc] are updated.

CALL prepro (blksiz, A, Prc, S)

!     Compute and store multi-level preconditioner at [Prc%mlp]:

CALL bepnum (blksiz, S, Prc)

CALL CPU_TIME(endtim)
  
  IF (outlev >= 3) THEN
    cumtim = endtim - begtim
    PRINT '(A, X, F7.2, /)' , 'Total time Preconditioner:', cumtim
  END IF
  
  
END SUBROUTINE cmpprc

END MODULE