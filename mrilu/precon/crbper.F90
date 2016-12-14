!#begindoc

#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix
#define prcmatrix  anymatrix

#endif

MODULE m_crbper

USE m_prcpars

CONTAINS

SUBROUTINE crbper (a, Ad, Prc)

USE m_dump
USE m_invblkcutmck 
USE m_schaak
USE m_build

TYPE (csrmatrix)				, POINTER		:: a
TYPE (diamatrix)				, POINTER		:: Ad
TYPE (prcmatrix)				, POINTER		:: Prc

!     Computes some Red-Black permutation, Prb, for the matrix  A,
!     and stores the permutationvector of Prb  into 'x%perrb'.

!     The  reordered matrix  B := Prb A (1/Prb)  will be partitioned as:
!            ( B_11 | B_12 )
!        B = (------+------)
!            ( B_21 | B_22 )
!     where the sub-matrix  B_11  is a non-singular block_diagonal
!     sub-matrix of  B.

!     The order of this sub-matrix  B_11  is stored in 'Prc%G'.

!     Arguments:
!     ==========
!     A%N	i   Number of rows/columns in the matrix  A.
!     Prc%G       o   Number of rows/columns in the sub_matrix  B_11.
!                   'Prc%G' is an integral multiple of 'Ad%blksiz'.
!     Ad%blksiz i   Number of rows/columns in the blocks of the block-
!                   diagonal matrix B_11.
!                   Both 'A%N' and 'Prc%G' are an integral multiple of 'Ad%blksiz'.
!     A%beg     i   A%beg(r), 1<=r<=A%N: Index in 'A%jco' of the first non-
!                   zero off-diagonal element in row  r  of matrix  A.
!                   A%beg(A%N+1): Index in 'A%jco' of the last nonzero
!                   element + 1.
!     A%jco     i   A%jco(nz), A%beg(r)<=nz<A%beg(r+1), 1<=r<=A%N: Column
!                   number of the nz-th non-zero off-diagonal element
!                   in row  r  of matrix  A.
!     Ad%co     i   The main (block-)diagonal of the matrix  A.
!                   Each block is stored in column major order.
!     Prc%perrb   o   Permutationvector of the Red-Black permutation, Prb.
!                   Prc%perrb(newrow), 1<=newrow<=A%N, is the old row/column
!                   number corresponding with row/column  newrow  in the
!                   transformed matrix   Prb A inv(Prb).

!#enddoc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'crbper'

!     Local Variables:
!     ================

INTEGER					:: ier
INTEGER 				:: i
INTEGER, ALLOCATABLE, DIMENSION(:)	:: perm

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( perm(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Calculate permutation, P1 in perm:
  
  IF (cutmck) THEN
!        Preprocessing with Reverse Block Cuthill-McKee ordering
    
    CALL invblkcutmck (Ad%blksiz, a, perm )
    
  ELSE
!        Construct identity permutation, P1 = I
    FORALL (i= 1: A%n) perm(i) = i
  END IF
  
  
  IF (XactElm) THEN

!   Calculate Red-Black permutation  Prb:

    CALL schaak (a, Ad, perm, Prc )

  ELSE

!   Copy RCM permutation into Red-Black permutation  Prb:

    Prc%perrb = perm
    
!   Indicate empty sub-matrix  B_11:

    Prc%g = 0
    Prc%nschur = A%n

  END IF

DEALLOCATE( perm, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  
!     Normal Return:
ier = 0
RETURN

!     Error Return, force value ier < 0:
1000 CONTINUE
ier = -1

!     End of  crbper
END SUBROUTINE crbper

END MODULE
