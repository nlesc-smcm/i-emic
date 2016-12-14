!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_invblkcutmck 

CONTAINS

SUBROUTINE invblkcutmck (BlkSiz, a, numRCM)

USE m_dump
USE m_build

INTEGER					, INTENT(IN)            :: BlkSiz
TYPE (csrmatrix)			, POINTER		:: a
INTEGER, DIMENSION(1:A%n)		, INTENT(OUT)           :: numRCM

!     Calculates the Reverse Cuthill-Mckee permutation, Prcm, which
!     gives the minimum band-width when applied to the matrix  A:
!        Prcm A inv(Prcm)
!     Signals if the matrix is reducible and prints the number of
!     independent systems, except if there is only one independent
!     system.
!     The reverse ordering, in general gives less fill-in than the
!     ordinary Cuthill-Mckee permutation.

!     Arguments:
!     ==========
!     A%N	i   Number of rows/colums in matrix A.
!     BlkSiz   	i   Number of rows/columns in a diagonal block.
!                   'A%N' is an integral multiple of 'BlkSiz'.
!     A%beg     i   A%beg(r), 1<=r<=A%N: Index in 'A%jco' of the first non-
!                   zero off-diagonal element in row  r  of matrix  A.
!                   A%beg(A%N+1): Index in 'A%jco' of the last nonzero
!                   element + 1.
!     A%jco     i   A%jco(nz), A%beg(r)<=nz<A%beg(r+1), 1<=r<=A%N: Column
!                   number of the nz-th non-zero off-diagonal element
!                   in row  r  of matrix  A.
!     numRCM   	o   Permutationvector of the Reverse Cuthill-McKee
!                   permutation, Prcm.
!                   numRCM(newrow), 1<=newrow<=A%N, is the old row/column
!                   number corresponding with row/column  newrow  in the
!                   transformed matrix   Prcm A inv(Prcm).

!#enddoc

!     Local Variables:
!     ================

INTEGER 		:: i, temp
INTEGER 		:: newblk, firrow, lasrow

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'invblkCutMck'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


!     Calculate the Block Cuthill-McKee ordering:
CALL blkcutmck (BlkSiz, a, numRCM)


!     Reverse the ordering:

!     First, Reverse the ordering on a block level
firrow = 1
lasrow = BlkSiz
DO newblk = 1, A%n/BlkSiz
  
!        Swap the pairs in a block:
  DO i = 1, BlkSiz/2
    temp               = numRCM(firrow-1+i)
    numRCM(firrow-1+i) = numRCM(lasrow+1-i)
    numRCM(lasrow+1-i) = temp
  END DO
  firrow = lasrow+1
  lasrow = lasrow + BlkSiz
END DO

!     Finally, Reverse the ordering. This combination maintains
!     the original equation ordering on a block level.

!     Swap the pairs:
DO i = 1, A%n/2
  temp          = numRCM(i)
  numRCM(i)     = numRCM(A%n+1-i)
  numRCM(A%n+1-i) = temp
END DO

!     End of  invblkcutmck
END SUBROUTINE invblkcutmck
!#begindoc

!=======================================================================

SUBROUTINE blkcutmck (BlkSiz, a, numcm)

USE m_dump
USE m_build

INTEGER					, INTENT(IN)	:: BlkSiz
TYPE (csrmatrix)			, POINTER	:: a
INTEGER, DIMENSION(1:A%n)		, INTENT(OUT)   :: numcm

!     Calculates the Cuthill-Mckee permutation, Pcm, of the matrix  A:
!        Prcm A inv(Prcm)
!     Signals if the matrix is reducible and prints the number of
!     independent systems, except if there is only one independent
!     system.

!     Arguments:
!     ==========
!     A%N      i   Number of rows/colums in matrix A.
!     BlkSiz   i   Number of rows/columns in a diagonal block.
!                  'A%N' is an integral multiple of 'BlkSiz'.
!     A%beg     i   A%beg(r), 1<=r<=A%N: Index in 'A%jco' of the first non-
!                  zero off-diagonal element in row  r  of matrix  A.
!                  A%beg(A%N+1): Index in 'A%jco' of the last nonzero
!                  element + 1.
!     A%jco     i   A%jco(nz), A%beg(r)<=nz<A%beg(r+1), 1<=r<=A%N: Column
!                  number of the nz-th non-zero off-diagonal element
!                  in row  r  of matrix  A.
!     numCM    o   Permutationvector of the the Cuthill-McKee
!                  permutation, Pcm.
!                  numCM(newrow), 1<=newrow<=A%N, is the old row/column
!                  number corresponding with row/column  newrow  in the
!                  transformed matrix   Pcm A inv(Pcm).

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'BlkCutMcK'

!     Local Variables:
!     ================
!     NSystm    Number of independent systems in matrix.
!     Seed      Seed for search process.
!     NxtSee    Next seed to restart process when matrix reducible.
!     NCheck    Number of equations checked.
!     TotLS     Number of equations in the first 'NSystem' level sets.
!     TotLS1    Number of equations in the first (reduced) 'NSystem'-1
!               level sets.
!     Chcked       Chcked(r) = .TRUE.:  Indication that row/column number
!                  'r' has been added to the level sets.

INTEGER 				:: ier, i
INTEGER 				:: nsystm
INTEGER 				:: seed, nxtsee, ncheck, totls, totls1
INTEGER 				:: fircol, firrow, nz
LOGICAL, ALLOCATABLE, DIMENSION(:)	:: chcked

!     Statement functions:
!     ====================

INTEGER :: corblk, corrow, blknr, rownr

!     CorBlk(rownr) returns the block number corresponding with row
!     'rownr'.

corblk(rownr) = ((rownr - 1) / BlkSiz + 1)

!     CorRow(blknr) returns the number of the first row/column in
!     diagonal block 'blknr'.

corrow(blknr) = ((blknr - 1) * BlkSiz + 1)

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( chcked(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

! Initialise:  Chcked(1:A%N) = .FALSE.

! Add rows/columns of the first block into the first level set:

  FORALL (i = 1:BlkSiz) numcm(i)   = i
  chcked(1:BlkSiz) = .true.
  chcked(BlkSiz+1:) = .false.
  totls  = BlkSiz
  totls1 = 0

!     Node to start search for another level set if system is reducible
nsystm = 1
nxtsee = 1

ncheck = 0
DO WHILE (totls < A%n)
  seed   = numcm(ncheck+1)
  firrow = seed
  
!        Find nodes connected to seed
  DO nz = A%beg(firrow), A%beg(firrow+BlkSiz) - 1
    fircol = corrow(corblk(A%jco(nz)))
    
    IF (.not. chcked(fircol)) THEN

!     Node not in a previous level set, add to current
!     collection of nodes in level sets:
      
      FORALL (i = 0:BlkSiz-1) numcm(totls+i+1)   = fircol+i
      chcked(fircol:fircol+BlkSiz-1) = .true.
      totls          = totls + BlkSiz

    END IF
  END DO
  
!        Increment number of equations checked
  ncheck = ncheck + BlkSiz
  
  IF (ncheck == totls) THEN
!           No connections (system is reducible)
#ifdef DEBUG
    IF (nsystm == 1) THEN
      PRINT '(/, A, 2X, A, A)', 'Warning from', rounam, ':  Matrix is reducible!'
      PRINT '(3X, I9, X, A, I9, X, A, /)', totls, 'of the', A%n, 'rows/columns in first set.'
    ELSE
      PRINT '(3X, I9, X, A, I9, X, A, I9, X, A, /)', totls1, 'of the', A%n, 'rows/columns in the first', nsystm, 'level sets.'
    END IF
#endif
  totls1 = totls
  nsystm = nsystm + 1
  
!           Find start block for next level set
  DO WHILE (chcked(nxtsee))
    nxtsee = nxtsee + 1
  END DO

  
! Add to current collection of nodes in level sets

  FORALL (i = 0:BlkSiz-1) numcm(totls+i+1)   = corrow(corblk(nxtsee))+i
  chcked(corrow(corblk(nxtsee)):corrow(corblk(nxtsee))+Blksiz-1) = .true.
  totls = totls + BlkSiz

END IF
END DO

DEALLOCATE( chcked, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

#ifndef DEBUG
IF (nsystm > 1) THEN
  PRINT '(/, A, 2X, A, A)', 'Warning from', rounam, ':  Matrix is reducible!'
  IF (nsystm-1 == 1) THEN
    PRINT '(3X, I9, X, A, I9, X, A, /)', totls1, 'of the', A%n, 'rows/columns in the first level set.'
  ELSE
    PRINT '(3X, I9, X, A, I9, X, A, I9, X, A, /)', totls1, 'of the', A%n, 'rows/columns in the first', nsystm-1, 'level sets.'
  END IF
END IF
#endif

!     End of  BlkCutMcK
END SUBROUTINE blkcutmck

END MODULE