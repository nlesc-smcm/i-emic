!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix
#define partmatrix  anymatrix
#define mlpmatrix  anymatrix

#endif

MODULE m_permpr

CONTAINS
 
SUBROUTINE permpr (P)

USE m_dump
USE m_build
USE m_mterrmsg
USE m_permrc

TYPE (mlpmatrix)	, POINTER		:: P

!     Permute row/column numbers in off-diagonal partitions.
!     Permute the column numbers in the right upper blocks and the row
!     numbers in the left lower blocks of all PLDU type partitions;
!     i.e. all except the last PFFP/PSFP partition, in the MLP type
!     matrix P.

!     Arguments:
!     ==========
!     P      i   Location of the descriptor of the MLP type
!                  matrix P.
!     P%perm     i   The inverse of permutation 'Perm' is applied to the
!                  row and column numbers in the above mentioned
!                  partitions.

!#enddoc

!     Array argument used as local variable:
!     ======================================
!     PerInv       Permutation to be applied to the row and column
!                  numbers in the above mentioned blocks.

!     Local Variables:
!     ================

INTEGER					:: ier
INTEGER 				:: i, matsiz
TYPE (partmatrix), POINTER		:: Par
TYPE (csrmatrix), POINTER		:: Upp
TYPE (cscmatrix), POINTER		:: Low
INTEGER, ALLOCATABLE, DIMENSION(:)	:: perinv

CHARACTER (LEN=*), PARAMETER :: rounam = 'permpr'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( perinv(1:P%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

! Construct Inverse of global Permutation vector 'Perm':

FORALL (i = 1:P%n) perinv(P%perm(i)) = i

! Loop over linked-list to consider each partition
! The last partition is not added to the list yet

! Locations of descriptors of first partition:

Par  => P%first

DO WHILE ( ASSOCIATED(Par))
  
! No permutation for Diagonal blocks!
! Location descriptor Upper triangular block:
  Upp => Par%utr
  
  IF (Upp%typ == csrtp) THEN

!   Upper Triangular block stored in CSR format
!           Get row/column indices and apply permutation PerInv:

    CALL permrc (Upp, perinv)

  ELSE

!   Stored in SYMMETRIC format, NOT IMPLEMENTED YET

    CALL mterrmsg (rounam, Upp%typ, symtp)
    CALL dump(__FILE__,__LINE__,'Not implemented yet')

  END IF
  
! Location descriptor Lower triangular block:

  Low => Par%ltr
  
  IF (Low%typ == csctp) THEN

!   Lower Triangular block stored in CSC format
!   Get column/row indices and apply permutation PerInv:

    CALL permrc (csctocsr(Low), perinv)

  ELSE

!   Stored in SYMMETRIC format, NOT IMPLEMENTED YET

    CALL mterrmsg (rounam, Upp%typ, symtp)
    CALL dump(__FILE__,__LINE__,'Not implemented yet')

  END IF
  
! Next partition

  Par => Par%next

END DO

DEALLOCATE( perinv, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')


END SUBROUTINE permpr

END MODULE