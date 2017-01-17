!#begindoc
 
#ifndef WITH_UNION

#define scdematrix  anymatrix
#define fmmatrix  anymatrix

#endif

MODULE m_scd2fm

CONTAINS
 
SUBROUTINE scd2fm (n, s, b)

USE m_build

INTEGER					, INTENT(IN)	:: n
TYPE (scdematrix)			, POINTER	:: s
TYPE (fmmatrix)				, POINTER	:: b

!     Copy  SCDE matrix to FMtp matrix.

!     Copy the  N x N  (sub-)matrix  S, stored in SCDE format
!        [s%offd%beg,s%offd%jco,s%offd%co,s%dia%com]
!     into the  N x N  matrix  B  in FMtp format.

!     Arguments:
!     ==========
!     N        i   Number of rows/columns in the matrices S and B.
!     S%dia%blksiz   i   Number of rows/columns in each diagonal block in
!                  block diagonal matrix S%dia.
!     s%offd%beg     i   s%offd%beg(r): Location in 's%offd%jco' and 's%offd%co' of first
!                  non-zero element in each row r of the off-diagonal
!                  part of matrix  S.
!     s%offd%jco     i   s%offd%jco(s%offd%beg(r):s%offd%beg(r+1)-1): Column numbers of the
!                  non-zero off-diagonal elements in row r of matrix S%OFFD.
!     s%offd%co      i   s%offd%co(s%offd%beg(r):s%offd%beg(r+1)-1): Values of the non-zero
!                  off-diagonal elements in row r of matrix  S%OFFD.
!     s%dia%com     i   The main (block-)diagonal entries of matrix S.
!                  Each block is stored in column major order.
!     B%COM        o   The matrix S in Full Matrix format.

!#enddoc

!     Local Variables:
!     ================
!     cnr          Actual column number in  S%OFFD  and  B%COM.
!     rnr          Actual row number in  S%OFFD  and  B%COM.
!     bnr          Actual block number in (block-)diagonal
!                  sub-matrix  S%dia.  Element  S(rnr,rnr)  is also an
!                  element of block  bnr  in S%dia.
!     firrnr       Row number in  S%OFFD  of first row in block  bnr.
!     lasrnr       Row number in  S%OFFD  of last row in block  bnr.

INTEGER :: bnr, cnr, rnr, firrnr, lasrnr, i

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'scd2fm'
#endif

!     Statement function:
!     ===================

INTEGER :: corblk, rownr

!     CorBlk(rownr) returns the block number corresponding with row
!     'rownr'.

corblk(rownr) = ((rownr - 1) / s%dia%blksiz + 1)

#ifdef DEBUG

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Initalise  B  with zeros:

!CALL dsprav (n**2, 0.0D0, b%com, 1)
b%com = 0.0D0

!     Copy  S  into  B:

DO rnr = 1, n
  
!        Copy the elements from the Block Diagonal in row 'rnr' first:
  
!        Number of diagonal block corresponding to row 'rnr' in S%dia:
  bnr  = corblk(rnr)
  
  firrnr = (bnr - 1) * s%dia%blksiz + 1
  lasrnr = (bnr - 1) * s%dia%blksiz + s%dia%blksiz
  
!        Copy row 'rnr' of block 'bnr' in  S%dia  into  B%COM:
  b%com(rnr,firrnr:lasrnr) = s%dia%com((rnr-firrnr+1),firrnr:lasrnr)
  
!        Copy the other non-zero off-diagonal elements of row 'rnr'
!        into  B%COM:

  FORALL (i = s%offd%beg(rnr):s%offd%beg(rnr+1)-1)  b%com(rnr,s%offd%jco(i))  = S%OFFD%CO(i)

END DO

END SUBROUTINE scd2fm

END MODULE