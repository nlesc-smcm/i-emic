!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define scdematrix  anymatrix

#endif

MODULE m_scd2csr

CONTAINS
 
SUBROUTINE scd2csr (s, b)

USE m_build

TYPE (scdematrix)				, POINTER	:: s
TYPE (csrmatrix)				, POINTER	:: b

!     Copy SCDE matrix to CSR matrix.

!     Copy the  N x N  (sub-)matrix  S, stored in SCDE format
!        [s%offd%beg,s%offd%jco,s%offd%co,s%dia%co]
!     into the  N x N  matrix  B  in CSR format
!        [b%beg,b%jco,b%co]

!     Arguments:
!     ==========
!     s%n        i   Number of rows/columns in the matrices S and B.
!     S%dia%blksiz   i   Number of rows/columns in each diagonal block in
!                  block diagonal matrix S%dia%.
!     s%offd%beg     i   s%offd%beg(r): Location in 's%offd%jco' and 's%offd%co' of first
!                  non-zero element in each row r of the off-diagonal
!                  part of matrix  S.
!     s%offd%jco     i   s%offd%jco(s%offd%beg(r):s%offd%beg(r+1)-1): Column numbers of the
!                  non-zero off-diagonal elements in row r of matrix S.
!     s%offd%co      i   s%offd%co(s%offd%beg(r):s%offd%beg(r+1)-1): Values of the non-zero
!                  off-diagonal elements in row r of matrix  S.
!     s%dia%co     i   The main (block-)diagonal entries of matrix S.
!                  Each block is stored in column major order.
!     b%beg     o   b%beg(r): Location in 'b%jco' and 'b%co' of first
!                  non-zero element in each row r of matrix  B.
!     b%jco     o   b%jco(b%beg(r):b%beg(r+1)-1): Column numbers of the
!                  non-zero elements in row r of matrix B.
!     b%co      o   b%co(b%beg(r):b%beg(r+1)-1): Values of the non-zero
!                  elements in row r of matrix  B.

!#enddoc

!     Local Variables:
!     ================
!     cnr          Actual column number in  S.
!     rnr          Actual row number in  S.
!     bnr          Actual block number in (block-)diagonal
!                  sub-matrix  S%dia%.  Element  S(rnr,rnr)  is also an
!                  element of block  bnr  in S%dia%.
!     firrnr       Row number in  S  of first row in block  bnr.

!     Local Variables:
!     ================

INTEGER :: bnr, cnr, rnr, firrnr, lasrnr, i
INTEGER :: nz, nnzb

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'scd2csr'
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

!     Copy  S  into  B:

nnzb   = 0
b%beg(1) = nnzb + 1

!     Copy each row r from  S  and  S%dia%  into  B:

DO rnr = 1, s%n
  
!        Copy the elements from the Block Diagonal first.
  
!        Number of diagonal block corresponding to row 'rnr' in  S%dia%:
  bnr    = corblk(rnr)
  
  firrnr = (bnr - 1) * s%dia%blksiz
  
! Copy row 'rnr' of block 'bnr' in  S%dia% into  B:

  FORALL (i = 1:s%dia%blksiz)
    b%jco(nnzb+i) = i+firrnr
    b%co(nnzb+i)  = s%dia%com(rnr-firrnr,i+firrnr)
  END FORALL
  nnzb       = nnzb +  s%dia%blksiz

! Copy the non-zero off-diagonal elements of row 'rnr' into
! row 'rnr' of  B:

  b%jco(nnzb+1:nnzb+s%offd%beg(rnr+1)-s%offd%beg(rnr)) = s%offd%jco(s%offd%beg(rnr):s%offd%beg(rnr+1)-1)
  b%co (nnzb+1:nnzb+s%offd%beg(rnr+1)-s%offd%beg(rnr))  = S%offd%co (s%offd%beg(rnr):s%offd%beg(rnr+1)-1)
  nnzb       = nnzb + s%offd%beg(rnr+1)-s%offd%beg(rnr)
  b%beg(rnr+1) = nnzb + 1
END DO

END SUBROUTINE scd2csr

END MODULE