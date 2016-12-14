!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix
#define csrdmatrix  anymatrix

#endif

MODULE m_reordrb


CONTAINS

SUBROUTINE reordrb (g, numrb, a, Ad, b)

USE m_build
USE m_glbpars
USE m_vispars
USE m_xfnminr
USE m_fstrlen
USE m_wrtmtd
USE m_ioerrmsg
USE m_hernumsch
USE m_permbdi

INTEGER						, INTENT(IN)            :: g
TYPE (csrmatrix) 				, POINTER		:: a
INTEGER, DIMENSION(1:A%n)			, INTENT(IN)            :: numrb
TYPE (diamatrix) 				, POINTER		:: Ad
TYPE (csrdmatrix)				, POINTER		:: b

!     Reorder the matrix  A  according to the Red-Black permutation
!     Prb  and store the transformed matrix into  B:
!        B := Prb A inv(Prb)

!     The matrix  B  can be partitioned as
!            ( B_11 | B_12 )
!        B = (------+------)
!            ( B_21 | B_22 )
!     with B_11 = B(1:G,1:G) is a block-diagonal submatrix of B.

!     Arguments:
!     ==========
!     A%N        i   Number of rows/columns in the matrix  A.
!     G        i   Number of rows/columns in the sub_matrix  B_11.
!     Ad%blksiz   i   Number of rows/columns in a diagonal block.
!                  Both 'A%N' and 'G' should be an integral multiple of
!                  'Ad%blksiz'.
!     numRB    i   Permutationvector of the Red-Black permutation, Prb.
!                  numRB(newrow), 1<=newrow<=A%N, is the old row/column
!                  number corresponding with row/column  newrow  in the
!                  transformed matrix   Prb A inv(Prb).
!     Ad%beg     i   Ad%beg(r), 1<=r<=A%N: Index in 'A%jco' and 'A%co' of the
!                  first non-zero off-diagonal element in row  r  of the
!                  matrix  A.
!                  Ad%beg(A%N+1): Index in 'A%jco' and 'A%co' of last nonzero
!                  off-diagonal element + 1.
!     A%jco     i   A%jco(nz), Ad%beg(r)<=nz<Ad%beg(r+1), 1<=r<=A%N: Column
!                  number of the nonzero off-diagonal element A%co(nz) in
!                  row  r  of the matrix  A.
!     A%co      i   A%co(Ad%beg(r):Ad%beg(r+1)-1), 1<=r<=A%N: Values of the
!                  nonzero off-diagonal elements in row  r  of the
!                  matrix  A.
!     Ad%co     i   The main (block-)diagonal of the matrix  A.
!                  Each block is stored in column major order.
!     B%offd%beg     o   B%offd%beg(r), 1<=r<=A%N: Index in 'B%offd%jco' and 'B%offd%co' of the
!                  first non-zero off-diagonal element in row  r  of the
!                  matrix  B.
!                  B%offd%beg(A%N+1): Index in 'B%offd%jco' and 'B%offd%co' of last nonzero
!                  off-diagonal element + 1.
!     B%offd%jco     o   B%offd%jco(nz), B%offd%beg(r)<=nz<B%offd%beg(r+1), 1<=r<=A%N: Column
!                  number of the nonzero off-diagonal element B%offd%co(nz) in
!                  row  r  of the matrix  B.
!     B%offd%co      o   B%offd%co(B%offd%beg(r):B%offd%beg(r+1)-1), 1<=r<=A%N: Values of the
!                  nonzero off-diagonal elements in row  r  of the
!                  matrix  B.
!     B%diA%co     o   The main (block-)diagonal of the matrix  B.
!                  Each block is stored in column major order.
!     B%lotr   o   B%OFFD%lotr(r), 1<=r<=G, index in 'B%offd%jco' and 'B%offd%co' of
!                  last nonzero in 'r'-th row of B12.
!                  B%OFFD%lotr(r), G+1<=r<=A%N, index in 'B%offd%jco' and 'B%offd%co' of
!                  last nonzero in 'r'-th row of B21.

!#enddoc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'reordrb'

!     Local Variables:
!     ================

CHARACTER (LEN=20) 	:: rofnm
INTEGER 		:: fnmlen, ier

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


!     Store reordered off-diagonal part of A, into the off-diagonal
!     part of B:
!        B := Prb A inv(Prb)
!     and return location of boundary between lower partition blocks.

!     Request workspace for INT Temporary array:

CALL hernumsch (g, numrb, a, b )

!     Store reordered Block Diagonal part of A, Ad, into the (block-)
!     Diagonal part of B, B%dia:
!        B%dia := P Ad inv(P)
CALL permbdi (A%n, numrb, Ad, B%dia)


IF (visnro > nrowrtn) THEN
!        Visualisation of reordered matrix  B  required.
  
  nrowrtn = nrowrtn + 1
  
  CALL xfnminr ('Breord', 2, nrowrtn, rofnm)
  fnmlen = fstrlen(rofnm)
  IF (outlev >= 1) THEN
    PRINT '(A, X, A)' , 'Visualize the reordered matrix with:  vsm', rofnm(1:fnmlen)
  END IF
  CALL wrtmtd (rofnm(1:fnmlen), .true., B%offd, B%dia, 2)
END IF


END SUBROUTINE reordrb

END MODULE