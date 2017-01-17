!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_redblack

CONTAINS
 
SUBROUTINE redblack (g, Nlim, blksiz, a, perRB )

USE m_dump
USE m_build

INTEGER				, INTENT(OUT)           :: g
INTEGER				, INTENT(IN)            :: Nlim
INTEGER				, INTENT(IN)            :: blksiz
TYPE (csrmatrix)		, POINTER		:: a
INTEGER, DIMENSION(:)		, INTENT(IN OUT)        :: perRB

!     Computes the permutation vector of the permutation  Prb, the
!     Red-Black order, so that the reordered matrix  A  is partitioned
!     as:
!                         ( B_11 | B_12 )
!        Prb A inv(Prb) = (------+------)
!                         ( B_21 | B_22 )
!        and
!        "the sub_matrix  B_11  is a non-singular block-diagonal matrix".

!     The returned permutationvector 'perRB' corresponds with the inverse
!     permutation  Prb'  of  Prb.

!     Arguments:
!     ==========
!     A%N       i   Number of rows/columns in the matrix  A.
!     G        	o   Number of rows/columns in the sub_matrix  B_11.
!                  'G' is an integral multiple of 'BlkSiz'.
!     BlkSiz   	i   Number of rows/columns in a diagonal block.
!                  'A%N' is an integral multiple of 'BlkSiz'.
!     A%beg     i   A%beg(r), 1<=r<=A%N: Index in 'A%jco' of the first non-
!                   zero off-diagonal element in row  r  of matrix  A.
!                   A%beg(A%N+1): Index in 'A%jco' of the last nonzero
!                   element + 1.
!     A%jco     i   A%jco(nz), A%beg(r)<=nz<A%beg(r+1), 1<=r<=A%N: Column
!                   number of the nz-th non-zero off-diagonal element
!                   in row  r  of matrix  A.
!     perRB    	io  In:  Permutationvector of some initial permutation
!                        Pini.  perRB(newrow), 1<=newrow<=UBOUND(perRB,1), where
!                        UBOUND(perRB,1)<=A%N, is the old
!                        row/column number corresponding with row/column
!                        newrow  in the transformed matrix Pini A inv(Pini).
!                   Out: Permutationvector of the Red-Black permutation
!                        Prb.  perRB(newrow), 1<=newrow<=UBOUND(perRB,1), is the old
!                        row/column number corresponding with row/column
!                        newrow  in the transformed matrix Prb A inv(Prb).

!#enddoc

!     Array arguments used as local variables:
!     ========================================
!     Cup2Up    In the first part of the subroutine:
!               Cup2Up(OldRow) = .TRUE.:  row 'OldRow' is directly
!               coupled to one of the rows in the upper partition.
!               In the last part:
!               Cup2Up contains the input value of the permutation
!               vector PerRB.
!     Inuppr    Inuppr(OldRow) = .TRUE.: Row 'OldRow' will be moved to
!               the upper partition.

!     Local Variables:
!     ================
!     Nlim	Number of rows in the upper partition of the matrix
!		to consider. Nlim = UBOUND(perRB,1) <= A%N. Operation is limited to the Nlim x Nlim sub-matrix.
!     FirRow    First Row in current diagonal block.
!     Nlow      Number of rows to be moved to the lower partition.
!     NewRow    Row number in the matrix Pini A inv(Pini).
!     OldRow    Row number in original matrix  A, corresponding with
!               row number 'NewRow' in the matrix  Pini A inv(Pini).

INTEGER				       	:: i, ier, firrow, newrow, oldrow, nlow, nupp
LOGICAL 				:: coupled, newcoupled
LOGICAL, ALLOCATABLE, DIMENSION(:)	:: cup2up, inuppr
INTEGER, ALLOCATABLE, DIMENSION(:)      :: hulp

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'redblack'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

IF (Nlim>A%N) CALL dump(__FILE__,__LINE__,'Wrong Nlim > A%N')

!     Initialise number of equations in upper partition and
!     number of equations at start of lower partition:

g = 0

ALLOCATE( cup2up(1:A%N), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( inuppr(1:A%N), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Data Initialisation:
cup2up = .false.
inuppr = .false.

!  Determine for each equation whether it is to be placed in
!  the upper or in the lower partition.

DO  firrow = 1, Nlim, blksiz
  !#begindoc

#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix

#endif

!        Check each row of the block
  DO newrow =firrow, firrow + blksiz - 1

    oldrow = perRB(newrow)
    
!           Is equation coupled directly to one of equations already in upper partition?
!           Is equation coupled via one of non-zeros in this equation to upper partition?

    coupled = cup2up(oldrow) .OR. ANY( inuppr( A%jco( A%beg(oldrow):A%beg(oldrow+1)-1 ) ) )
    IF (coupled) EXIT

  END DO
  
  IF (.NOT. coupled) THEN

!   Rows in current block will reside in upper partition:

    DO newrow = firrow, firrow + blksiz - 1

      oldrow = perRB(newrow)
      
!     Increment number of equations in upper partition:

      g        = g + 1
      
!     OldRow becomes part of new upper partition:

      inuppr(oldrow) = .true.
      
!     Indicate that corresponding equations are coupled to upper partition:

      cup2up(A%jco(A%beg(oldrow):A%beg(oldrow+1)-1)) = .true.

    END DO
  END IF
END DO


! Update the permutation vector 'perRB' to reflect the new permutation:
DEALLOCATE( cup2up, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

ALLOCATE( hulp(1:A%N), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')


hulp = perRB
nupp = 0
nlow = g

DO newrow = 1, a%n

! 'hulp' contains the initial permutation in 'perRB'!

  oldrow = hulp(newrow)
  IF (inuppr(oldrow) ) THEN
    nupp = nupp + 1
    perRB(nupp) = oldrow
  ELSE
    nlow = nlow + 1
    perRB(nlow) = oldrow
  END IF

END DO

DEALLOCATE( hulp, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( inuppr, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

!     End of  redblack
END SUBROUTINE redblack

END MODULE