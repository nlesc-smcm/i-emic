!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix

#endif

MODULE m_csr2csc

CONTAINS

!     Originally subroutine  csrcsc2  from  sparskit:
!     ===============================================

SUBROUTINE csr2csc (ma, na, job, a, b)

USE m_build

INTEGER				, INTENT(IN)		:: ma
INTEGER				, INTENT(IN)            :: na
INTEGER				, INTENT(IN)            :: job
TYPE (cscmatrix)		, POINTER		:: a
TYPE (csrmatrix)		, POINTER		:: b

!     Compressed Sparse Row     to      Compressed Sparse Column

!     Convert rectangular  mA x nA  matrix  A  from CSR format to CSC
!     format in the  nA x mA  matrix  B:
!        call csr2csc (mA, nA, 1, A, B)

!     This routine can also be used to convert the  m x n  matrix  B
!     in CSC format to the  n x m  matrix  A  in CSR format:
!        call csr2csc (n, m, 1, B, A)

!     Transposition of the  m x n  matrix  A  in CSR format with the
!     transposed matrix stored in the  n x m  matrix  B  in CSR format
!     is realized by:
!        call csr2csc (m, n, 1, A, B)
!     Transposition operation:   Not in place!

!     Arguments:
!     ==========
!     mA       i   Number of rows in matrix  A.
!     nA       i   Number of columns in matrix  A.
!     job      i   Indicate whether to fill the values or not.
!                  = 1  fill pattern in  b%beg, b%jco  and values in  B,
!                  = 2  fill pattern in  b%beg, b%jco  only.
!     a%beg    i   Array of length  mA+1.  a%beg(i) contains the position
!                  of the Beginning of Row segment i in  A  and  a%jco.
!     a%jco    i   Array of length nnz containing the column numbers
!                  of the corresponding elements in  A.
!     a%co     i   Array of length nnz (nnz=number of nonzero elements
!                  in input matrix) containing the nonzero elements.

!     b%beg    o   Array of length nA+1.  b%beg(j) contains the position
!                  of the Beginning of Column segment j in  B  and  b%jco.
!     b%jco    o   Array of size nnz containing the row numbers of the
!                  corresponding array elements in  B.
!                  Note: the row numbers in one column segment are in
!                  increasing order.
!     b%co     o   If job = 1:
!                  Array of size nnz containing the nonzero elements of
!                  the matrix  A, stored columnwise.
!                  If job = 2:
!                  Value is not changed!

!#enddoc

!     Local variables:
!     ================
INTEGER :: i, j, ij, next

!     Compute the lengths of the columns of A in b%beg(2:nA+1):

b%beg = 0

DO  i = 1, ma
  FORALL (ij = a%beg(i):a%beg(i+1) - 1) b%beg(a%jco(ij) + 1) = b%beg(a%jco(ij) + 1) + 1
END DO

!     {  A(j: 1 <= j <= nA:  b%beg(j+1) = "length j-th column in A")  }

!     Compute pointers to beginning of column segments from lengths:

b%beg(1) = 1
DO i = 1, na
  b%beg(i+1) = b%beg(i+1) + b%beg(i)
END DO
b%nnz = b%beg(na+1)-b%beg(1);

!     Now do the actual copying:

!     {  A(j: 1 <= j <= nA: b%beg(j) = "position next elmt in column j") }
IF ( job == 1 ) THEN
  DO  i = 1, ma
    DO  ij = a%beg(i), a%beg(i+1) - 1
      j         = a%jco(ij)
      next      = b%beg(j)
      b%co(next)   = a%co(ij)
      b%jco(next) = i
      b%beg(j)    = next + 1
    END DO
  END DO
ELSE
  DO  i = 1, ma
    DO  ij = a%beg(i), a%beg(i+1) - 1
      j         = a%jco(ij)
      next      = b%beg(j)
      b%jco(next) = i
      b%beg(j)    = next + 1
    END DO
  END DO
END IF

!     Shift down b%beg, b%beg(2:nA+1) = b%beg(1:nA), and restore b%beg(1):

DO  j = na, 1, -1
  b%beg(j+1) = b%beg(j)
END DO
b%beg(1) = 1
b%nnz = b%beg(na+1) - b%beg(1) 

!     End of  csr2csc
END SUBROUTINE csr2csc

END MODULE

