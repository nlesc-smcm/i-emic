!#begindoc
 
#ifndef WITH_UNION

#define diamatrix  anymatrix

#endif

MODULE m_permbdi

CONTAINS
 
SUBROUTINE permbdi (n, numrb, Ad, Bd )

USE m_build

INTEGER						, INTENT(IN)            :: n
INTEGER, DIMENSION(1:n)				, INTENT(IN)            :: numrb
TYPE (diamatrix)				, POINTER         	:: Ad
TYPE (diamatrix)				, POINTER         	:: Bd

!     The permuted Block Diagonal matrix  Ad  is copied into the
!     Block Diagonal matrix  Bd, with

!        Bd  =  Prb Ad inv(Prb)

!     Arguments:
!     ==========
!     N			i   Number of rows/columns in DIAtp matrix Ad.
!     Ad%blksiz   	i   Number of rows/columns in the diagonal blocks in Ad  and  Bd.
!			   'N' should be an integer multiple of 'Ad%blksiz'!
!     numRB    		i   Permutationvector of the Red-Black permutation, Prb.
!                  	    numRB(newrow), 1<=newrow<=N, is the old row/column
!                  	    number corresponding with row/column  newrow  in the
!                  	    transformed matrix   Prb Ad inv(Prb).
!     Ad%com     	i   The values of the elements of the block diagonal of Ad, stored in column major order.
!     Bd%com     	o   The values of the elements of the block diagonal of Bd, stored in column major order.

!#enddoc

!     Local Variables:
!     ================

INTEGER :: i

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'permbdi'
#endif

#ifdef DEBUG

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif 

!     Store the diagonal blocks in B in the order prescribed in numRB.
!     The position of the blocks in B may be different from the position
!     in A, but the relative positions of the rows/columns in a block
!     remain unchanged.

DO i = 1, n, Ad%blksiz
  Bd%com(1:Ad%blksiz,i:i+Ad%blksiz-1) = Ad%com(1:Ad%blksiz,numrb(i):numrb(i)+Ad%blksiz-1)
END DO

END SUBROUTINE permbdi

END MODULE