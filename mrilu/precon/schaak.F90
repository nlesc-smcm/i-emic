!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix
#define prcmatrix  anymatrix

#endif

MODULE m_schaak

CONTAINS
 
SUBROUTINE schaak (a, Ad, numRCM, x)

USE m_dump
USE m_redblack
USE m_dgeco
USE m_build

TYPE (csrmatrix)		, POINTER		:: a
TYPE (diamatrix)		, POINTER		:: Ad
INTEGER, DIMENSION(1:a%n)	, INTENT(IN)            :: numRCM
TYPE (prcmatrix)		, POINTER		:: x

!     Computes the permutation vector of the permutation  Prb, the
!     Red-Black order, so that the reordered matrix  A  is partitioned
!     as:
!                         ( B_11 | B_12 )
!        Prb A inv(Prb) = (------+------)
!                         ( B_21 | B_22 )
!        and
!        "the sub_matrix  B_11  is a non-singular block-diagonal matrix".
!     The returned permutationvector 'x%perrb' corresponds with the
!     permutation  Prb.

!     Arguments:
!     ==========
!     A%N      	i   Number of rows/columns in the matrix  A.
!     X%G	o   Number of rows/columns in the sub_matrix  B_11.
!                  'X%G' is an integral multiple of 'Ad%blksiz'.
!     Ad%blksiz i   Number of rows/columns in a diagonal block.
!                  'A%N' is an integral multiple of 'Ad%blksiz'.
!     a%beg     i   a%beg(r), 1<=r<=A%N: Index in 'a%jco' of the first non-
!                   zero off-diagonal element in row  r  of matrix  A.
!                   a%beg(A%N+1): Index in 'a%jco' of the last nonzero
!                   element + 1.
!     a%jco     i   a%jco(nz), a%beg(r)<=nz<a%beg(r+1), 1<=r<=A%N: Column
!                   number of the nz-th non-zero off-diagonal element
!                   in row  r  of matrix  A.
!     Ad%com     i   The main (block-)diagonal of the matrix  A.
!                   Each block is stored in column major order.
!     numRCM   	i   Permutationvector of the Reverse Cuthill-McKee
!                   permutation, Prcm.
!                   numRCM(NewRow), 1<=NewRow<=A%N, is the old row/column
!                   number corresponding with row/column  NewRow  in the
!                   transformed matrix   Prcm A inv(Prcm).
!     x%perrb  	o   Permutationvector of the Red-Black permutation, Prb.
!                   x%perrb(NewRow), 1<=NewRow<=A%N, is the old row/column
!                   number corresponding with row/column  NewRow  in the
!                   transformed matrix   Prb A inv(Prb).

!#enddoc

!     Array arguments used as local variables:
!     ========================================
!     nnzic     1. nnzic(OldCol) = Number of non-zeros in column
!                  'OldCol'.
!               2. Work array in subroutine 'redblack'.
!     Moved     1. Moved(NewRow) = .TRUE.:  Row 'NewRow' Moved to lower
!                  partition.
!               2. Work array in subroutine 'redblack'.

!     Local Variables:
!     ================
!     Nlim      Number of rows in the upper partition of the matrix
!               to consider.
!     Nlow      Number of rows to be Moved to the lower partition.
!     NewRow    Row number in the matrix Prcm A inv(Prcm).
!     OldRow    Row number in original matrix  A, corresponding with
!               row number 'NewRow' in the matrix Prcm A inv(Prcm).

INTEGER 						:: FirRow, LasRow, ier, Nreord, i, annz
INTEGER 						:: NewRow, OldRow, nz, largnz, Nlim
DOUBLE PRECISION 					:: rcond
LOGICAL 						:: MustMove
LOGICAL, ALLOCATABLE, DIMENSION(:)            		:: Moved
INTEGER, ALLOCATABLE, DIMENSION(:)             		:: nnzic
INTEGER, DIMENSION(1:Ad%blksiz)             		:: iwrk
DOUBLE PRECISION, DIMENSION(1:Ad%blksiz,1:Ad%blksiz), TARGET	:: AdCoCopy
DOUBLE PRECISION, DIMENSION(:,:), POINTER		:: AdCo

CHARACTER (LEN=*), PARAMETER :: rounam = 'schaak'


IF ( a%n /= Ad%n ) STOP 'Abnormal temination in schaak: incompatible dimensions of a and Ad'
#ifdef DEBUG
! TRACE INFORMATION
  PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( Moved(1:a%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( nnzic(1:a%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Data Initialisation
Moved = .false.

Nlim   = a%n

!     Put equations corresponding with singular diagonal block to the
!     end of the lower partition:
IF (Ad%blksiz == 1) THEN
  
  DO NewRow = 1, a%n
    OldRow = numRCM(NewRow)
    IF (ABS(Ad%com(1,OldRow)) <= 1.0D-8) THEN
!     Store permutation index at end of lower partition.
#ifdef DEBUG
      PRINT '(A, 2X, A, /, 3X, A)' , 'Warning in', rounam, 'Zero on main diagonal, row put to end.'
#endif
      x%perrb(Nlim)   = OldRow
      Nlim          = Nlim - 1
      Moved(NewRow) = .true.
    END IF
  END DO

ELSE
!        {  Ad%blksiz > 1  }

!        For each block
  DO FirRow = 1, a%n, Ad%blksiz
    LasRow = FirRow + Ad%blksiz - 1
  
    AdCoCopy = Ad%com(1:Ad%blksiz,FirRow:LasRow)
  
!           Factor the matrix and estimate condition number:
    CALL dgeco (AdCoCopy, Ad%blksiz, iwrk, rcond )    
  
    IF (rcond <= 1.0D-8) THEN
!              Condition number of diagonal block >= 1.0D+3.
!              Store permutation indices at end of lower partition.
#ifdef DEBUG
      PRINT '(A, 2X, A, /, 3X, A)' , 'Warning in', rounam, 'Singular block on main diagonal, rows put to end.'
#endif

      Nlim                           = Nlim - Ad%blksiz
      x%perrb(Nlim+1:Nlim+Ad%blksiz) = numRCM(FirRow:LasRow)
      Moved(FirRow:LasRow)           = .true.

   END IF
  END DO
END IF

!     Compute number of non-zeros in columns (is known for rows):

nnzic = 0
annz = csrnnz(a)
nnzic(a%jco(1:annz)) = nnzic(a%jco(1:annz)) + 1
!FORALL (i=1:a%n) nnzic(a%jco(a%beg(i):a%beg(i+1)-1)) = nnzic(a%jco(a%beg(i):a%beg(i+1)-1)) + 1

!     Put relatively full rows and columns to end of lower partition.
!     Relatively full means that total number of non-zeros is more than
!     7 times the average number of non-zeros in each row.

largnz = 7.0D0*(annz)/DBLE(a%n)

DO FirRow = 1, a%n, Ad%blksiz
  LasRow = FirRow + Ad%blksiz - 1
  
!        Only if not already permuted to end of lower partition.
  IF (.not. Moved(FirRow)) THEN
    
!           Check each row of the block

    MustMove = .false.
    DO NewRow = FirRow, LasRow
      OldRow = numRCM(NewRow)
      IF ( (a%beg(OldRow+1)-a%beg(OldRow) > largnz) .OR. (nnzic(OldRow) > largnz) ) THEN
        MustMove = .true.
        EXIT
      END IF
    END DO
      
    IF (MustMove) THEN
!              Block must be Moved:
#ifdef DEBUG
      PRINT '(A, 2X, A, /, 3X, A)' , 'Warning in', rounam, 'Row(s) too full, put to end.'
#endif

!     Store permutation directly in right order:

      Nlim                           = Nlim - Ad%blksiz
      x%perrb(Nlim+1:Nlim+Ad%blksiz) = numRCM(FirRow:LasRow)
      Moved(FirRow:LasRow)           = .true.

    END IF
  END IF
END DO

Nreord = Nlim

!     Update the permutation vector  x%perrb(1:Nreord):
DO NewRow = a%n, 1, -1
  IF (.not. Moved(NewRow)) THEN
    x%perrb(Nlim) = numRCM(NewRow)
    Nlim        = Nlim - 1
  END IF
END DO
#ifdef DEBUG
IF (Nlim /= 0) THEN
  PRINT '(A, 2X, A, /, 3X, A, I11)' , 'Internal error in', rounam, 'Illegal value Nlim =', Nlim 
  STOP 'In schaak'
END IF
#endif

DEALLOCATE( Moved, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( nnzic, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

!     Determine Permutation vector with red-black order:

CALL redblack (x%g, Nreord, Ad%blksiz, a, x%perrb)

x%nschur = x%n - x%g


END SUBROUTINE schaak

END MODULE
