!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_order

CONTAINS
 
SUBROUTINE order (blksiz, nupp, a, newper, cspace, rspace )

USE m_dump
USE m_prcpars
USE m_redblack
USE m_build
USE m_wacsr
USE m_wfree

INTEGER					, INTENT(IN)            :: blksiz
INTEGER					, INTENT(OUT)        	:: nupp
TYPE (csrmatrix)			, POINTER		:: a
INTEGER, DIMENSION(1:A%n)		, INTENT(OUT)           :: newper
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)        	:: cspace
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)        	:: rspace

!     Computes a new order of rows/columns in the  N x N  sub-
!     matrix  A,  and puts the reordering permutation in the vector
!     NewPer.
!     Note: The relative order of the rows/columns within each
!           diagonal block of A remains unchanged!

!     Let  P  be a permutation matrix and let

!        S := P' A P            ( or    A = P S P' )

!     with

!              ( S11 | S12 )
!        S  =  ( ----+---- )
!              ( S21 | S22 )

!     and the reordered sub-matrix, S, is partitioned in such a way
!     that the block  S11  in the left upper part is a strong diagonal
!     dominant sub-matrix.
!     The permutation matrix, P, is characterised by the computed vector
!     'NewPer', i.e. NewPer(1:A%N)' = P (1:A%N)'.

!     Arguments:
!     ==========
!     A%N        i   Number of rows/columns in the sub-matrix A.
!     BlkSiz   i   Number of rows/columns in a diagonal block.
!                  'A%N' is an integral multiple of 'BlkSiz'.
!     Nupp     o   Number of rows in the upper partitions S11 and S12
!                  of matrix S.
!                  'Nupp' is an integral multiple of 'BlkSiz'.
!     A%beg     i   A%beg(i): the index in A%jco and A of the first
!                  off-diagonal element in row i.
!     A%jco     i   A%jco(nz): the column number of off-diagonal element
!                  A%co(nz).
!     A%co      i   A%co(nz):  the value of an off_diagonal element.
!     NewPer   o   Permutation of  1:A%N, with
!                  NewPer(i): row/column number in matrix A corresponding
!                  with row/column i in partitioned matrix S.
!                  The row numbers  NewPer(1:Nupp)  of A will belong to
!                  the upper part of S.
!     CSpace   i   CSpace(1:A%N) Column lump space in A.
!     RSpace   i   RSpace(1:A%N) Row lump space in A.

!#enddoc

!     Global Parameters:
!     ==================

!     Array arguments used as local variables:
!     ========================================
!     InUppr       1. A work array for subroutine 'redblack'.
!                  2. InUppr(r):  Indication whether row r is added to
!                     the upper partition; value .TRUE. or .FALSE.
!     Perm         1. A work array for subroutine 'redblack'.
!                  2. The permutation vector of the inverse permutation
!                     of the permutation in 'NewPer'.
!                  3. A copy of the initial value of 'NewPer'.
!     ar%beg       ar%beg(i): the index in ar%jco of the first large
!                  off-diagonal element of A in row i.
!     ar%jco        ar%jco(lnz): the column number of a large off-diagonal
!                  element of A.
!     sumic        Sum in Column of off-diagonal elements in rows in new
!                  upper border.
!                  sumic(j): Sum(i: 1<=new(i)<=Upp , off-diagonal(i,j) :
!                                   abs(A(i,j)) )
!     sumir        Sum in Row of off-diagonal elements in columns of new
!                  left border.
!                  sumir(i): Sum(j: 1<=new(j)<=Upp , off-diagonal(i,j) :
!                                   abs(A(i,j)) )

!     Local Variables:
!     ================
!     CanAdd       Indicates if rows in actual block can be added to
!                  the upper part.
!     ncr          New Column/Row number in (reordered) matrix  A.
!     ocr          Old Column/Row number in (original) matrix  A.
!     FirCol       First column, or row, number in a block.
!     LasCol       Last column, or row, number in a block.
!     Nlow         Number of rows in lower partition.
!     Nz           Index of nonzero entry in A%jco, A%co and later Row.

LOGICAL 					:: canadd
INTEGER 					:: fircol, lascol, newupp, nlow, nz, nzle, ier
INTEGER 					:: oldcol, ncr, ocr
INTEGER 					:: newrow, oldrow, i
TYPE (csrmatrix), POINTER			:: ar
LOGICAL, ALLOCATABLE, DIMENSION(:)		:: InUppr
INTEGER, ALLOCATABLE, DIMENSION(:)		:: perm
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: sumic, sumir

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'order'
#endif


#ifdef DEBUG

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Allocate the storage for the CSR representation of the matrix:
CALL wacsr (A%n, csrnnz(A), ar)

!     Determine the intial block partition of the matrix A:

!     Store the positions of the "large" elements of A in 'ar%beg'
!     and 'ar%jco':
nzle = A%beg(1)
ar%beg(1) = nzle
DO oldrow = 1, A%n
  DO nz = A%beg(oldrow), A%beg(oldrow+1)-1
    oldcol = A%jco(nz)
    IF ( ABS(A%co(nz)) >= elmfctr*redfctr*MIN(cspace(oldcol),rspace(oldrow)) ) THEN
    ar%jco(nzle) = oldcol
    nzle = nzle + 1
  END IF
END DO
ar%beg(oldrow+1) = nzle
END DO

ALLOCATE( InUppr(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( perm(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( sumic(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( sumir(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Initialize Permutation with identity:
FORALL (i = 1:A%n)  newper(i) = i

!     Determine the permutation for the initial block partition
!     ('InUppr' and 'Perm' are used as work arrays!):
CALL redblack (nupp, A%N, blksiz, ar, newper)

!     Deallocate the storage for the CSR representation of the matrix:
CALL csrfree(ar)

!     Compute the sums in 'sumir' and 'sumic':
sumic = 0.0D0
sumir = 0.0D0

!     Store inverse of the permutation 'NewPer' in 'Perm':
FORALL (i = 1:A%n) perm(newper(i)) = i

DO newrow = 1, nupp
  oldrow = newper(newrow)
  DO nz = A%beg(oldrow), A%beg(oldrow+1)-1
    oldcol = A%jco(nz)
    IF (perm(oldcol) <= nupp) THEN
      sumic(oldcol) = sumic(oldcol) + ABS(A%co(nz))
      sumir(oldrow) = sumir(oldrow) + ABS(A%co(nz))
    END IF
  END DO
END DO


!     For each diagonal block in the preliminary upper part, determine
!     whether its addition will maintain diagonal dominance:
nlow   = nupp
DO fircol = 1, nupp, blksiz
  lascol = fircol + blksiz - 1
  
!        Check if current block is a candidate for addition:
  DO ncr = fircol, lascol
    ocr = newper(ncr)
    canadd = (sumir(ocr) < redfctr*rspace(ocr)) .ANd. (sumic(ocr) < redfctr*cspace(ocr))
    IF (.NOT. canadd) EXIT
  END DO

  IF (canadd) THEN
    InUppr(fircol:lascol) = .true.
  ELSE
    nlow = nlow - blksiz
    InUppr(fircol:lascol) = .false.
  END IF
END DO


!     Update the permutation 'NewPer' and the sums in 'sumic' and
!     'sumir':
IF (nlow < nupp) THEN
!        Save 'NewPer(1:Nupp)' into 'Perm':
   perm = newper
  
  newupp = 0
  DO ncr = 1, nupp
    IF ( InUppr(ncr) ) THEN
      newupp         = newupp + 1
      newper(newupp) = perm(ncr)
    ELSE
!              InUppr(ncr) .EQ. .FALSE.
      nlow         = nlow + 1
      newper(nlow) = perm(ncr)
    END IF
  END DO
  
!        Adjust value 'Nupp':
  nupp = newupp
END IF

DEALLOCATE( InUppr, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( perm, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( sumic, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( sumir, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

!     End of  order
END SUBROUTINE order

END MODULE