!#begindoc
 
#ifndef WITH_UNION

#define partmatrix  anymatrix

#endif

MODULE m_wapldu

CONTAINS

SUBROUTINE wapldu (n, blksiz, nnzu, nnzl, offset, x)

USE m_dump
USE m_build
USE m_wacsc
USE m_wacsr
USE m_wadia

INTEGER, INTENT(IN)                  :: n
INTEGER, INTENT(IN)                  :: blksiz
INTEGER, INTENT(IN)                  :: nnzu
INTEGER, INTENT(IN)                  :: nnzl
INTEGER, INTENT(IN)                  :: offset
TYPE (partmatrix), POINTER            :: x

!     Allocate a partition
!     descriptor and the segments for a PLDU type partition in a MLP
!     type matrix.
!     The descriptor is filled so that the requested segments can be
!     referenced through the descriptor at 'x'.

!     Arguments:
!     ==========
!     N     i   Number of unknowns and equations represented by the
!                  partition.
!     BlkSiz   i   Number of rows/columns in a diagonal block.
!                  'N' should be an integer multiple of 'BlkSiz'!
!     NnzU     i   Number of non-zeros in matrix  U
!     NnzL     i   Number of non-zeros in matrix  L
!     Offset   i   Offset of partition in MLP matrix  P.
!     x   o   Location of the Partition descriptor.

!#enddoc

CHARACTER (LEN=*), PARAMETER 	:: rounam = 'wapldu'
INTEGER				:: ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam

PRINT '(A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)', 'Allocating ', n, &
 'x', n, ' PLDU-matrix with ', nnzu, ' non-zeros in the U-matrix, '&
  , nnzl, ' non-zeros in the L-matrix, offset ', offset, &
   ' and blocksize', blksiz
#endif

!     Request segment for partition descriptor.
ALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Allocate storage for a sparse block diagonal matrix (type DIAtp)
CALL wadia (n, blksiz, x%dia)

!     Allocate storage for a CSR matrix (type: CSR)
CALL wacsr (n, nnzu, x%utr)

!     Allocate storage for a sparse matrix in CSC format
CALL wacsc (n, nnzl, x%ltr)

!     Save information in partition descriptor:
!     Partition Type, Offset for block partition in root matrix and the
!     locations of the matrix block partitions:

x%typ   = pldutp
x%off  = offset
NULLIFY(x%prev)
NULLIFY(x%next)

END SUBROUTINE wapldu

END MODULE
