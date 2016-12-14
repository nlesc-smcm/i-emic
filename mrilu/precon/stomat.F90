!#begindoc
 
#ifndef WITH_UNION

#define csrdmatrix  anymatrix
#define csrmatrix  anymatrix
#define scbmmatrix  anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_stomat

CONTAINS
 
SUBROUTINE stomat (n, g, a, b)

USE m_build
USE m_wacsr
USE m_wascbm

INTEGER						, INTENT(IN)		:: n
INTEGER						, INTENT(IN)            :: g
TYPE (csrdmatrix)				, POINTER	        :: a
TYPE (anymatrix)				, POINTER		:: b

!     Store matrix 'A', [A, A%dia], into matrix 'B', [b].
!     The matrix 'B' will be stored either in CSR format (G = 0) or in SCBM format (G > 0).

!     If G = 0:
!     Stores the matrix 'A' into matrix 'B' in CSR format:
!        B = A + A%dia

!     If G > 0:
!     Stores the matrix 'A' into matrix 'B' in SCBM format:
!        B_11(1:G,1:G)     = B(1:G,1:G)     = A%dia(1:G,1:G)
!        B_12(1:G,1:N-G)   = B(1:G,G+1:N)   = A(1:G,G+1:N),
!        B_21(1:N-G,1:G)   = B(G+1:N,1:G)   = A(G+1:N,1:G),
!        B_22(1:N-G,1:N-G) = B(G+1:N,G+1:N)
!                          = A(G+1:N,G+1:N) + A%dia(G+1:N,G+1:N)

!     Arguments:
!     ==========
!     N        		i   Number of rows/columns in the matrix  A.
!     G        		i   Number of rows/columns in the left upper part of matrix A.
!     A%dia%blksiz   	i   Number of rows/columns in a diagonal block.
!                           Both 'N' and 'G' should be an integer multiple of 'A%dia%blksiz'.
!     a%offd%beg     	i   a%offd%beg(i): index in 'a%offd%jco' and 'a%offd%co' of the first
!                           off-diagonal element in row i of matrix A.
!     a%offd%jco     	i   a%offd%jco(nz): column number of non-zero off-diagonal
!                           element  a%offd%co(nz)  of matrix  A.
!     a%offd%co      	i   a%offd%co(nz): value of non-zero off-diagonal element of matrix  A.
!     a%dia%com     	i   a%dia%com(be): value of block-diagonal element of matrix  A.
!     A%lotr    	i   A%lotr(i), G+1 <= i <= N, index in 'a%offd%jco' and
!                           'a%offd%co' of last nonzero in 'i'-th row of [A21 A22].
!     b      		o   Location of the descriptor of a copy of the matrix A.
!                           The copy B of A is stored in CSR format if G = 0 and
!                           the copy B is stored in SCBM format if G > 0.

!#enddoc

!     Local Variables:
!     ================
!     Nnz        Number of non-zeros stored in  B, [begB,jcoB,coB]
!     BDval      Value of actual entry in block-diagonal matrix A%dia

TYPE (scbmmatrix), POINTER			:: xb
INTEGER 					:: nz12, nz21, nz22, x12nnz, x21nnz, x22nnz, BlkSiz
INTEGER 					:: col, fircol, nz, row
DOUBLE PRECISION 				:: bdval
TYPE (csrmatrix), POINTER			:: x22, x12, x21
TYPE (diamatrix), POINTER  			:: xb11d


#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'stomat'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


BlkSiz = a%dia%blksiz

IF (g == 0) THEN
  
!        Request storage for the CSR matrix B:
  nz22 = a%offd%beg(n+1) - a%offd%beg(1) + n * BlkSiz
  
  CALL wacsr (n, nz22, x22)
  b => csrtoany(x22)
  
!        Matrix B is treated as sub-matrix B22
  x22nnz = 0
  x22%beg(1) = x22nnz + 1
ELSE
!        0 < G < N
!        Request storage for the SCBM matrix B:
  nz12 = a%offd%beg(g+1) - a%offd%beg(1)
  nz21 = SUM(a%lotr(g+1:n)) - SUM(a%offd%beg(g+1:n)) + n - g
  nz22 = a%offd%beg(n+1) - a%offd%beg(g+1) - nz21 + (n - g) * BlkSiz
 
  CALL wascbm (n, g, BlkSiz, nz12, nz21, nz22, xb)

  b => scbmtoany(xb)
  
!        Store the, left upper, block diagonal, partition:
!           B_11 = B(1:G,1:G) = A%dia(1:G,1:G)
  xb11d => xb%a11d
  xb11d%com = a%dia%com(1:BlkSiz,1:g)
  
  
!        Store the remainder of matrix 'A' into matrix 'B':
  
!        Right upper partition:
!           B_12 = B(1:G,G+1:N) = A(1:G,G+1:N)
  
  x12 => xb%a12

  
  x12%beg(1:g+1) = a%offd%beg(1:g+1)
!  x12nnz         = x12%beg(g+1)-x12%beg(1)
  x12nnz         = SUM(x12%beg(2:g+1)-x12%beg(1:g))
  
  x12%jco(1:x12nnz) = a%offd%jco(1:x12nnz) - g
  x12%co(1:x12nnz)  = a%offd%co(1:x12nnz)
  
  
!        Lower partitions:
!           B_21 = B(G+1:N,1:G)   = A(G+1:N,1:G)
!           B_22 = B(G+1:N,G+1:N) = A(G+1:N,G+1:N) + A%dia(G+1:N,G+1:N)
  
  x21 => xb%a21
  x22 => xb%a22

  x21nnz = 0
  x22nnz = 0
  x21%beg(1) = x21nnz + 1
  x22%beg(1) = x22nnz + 1
END IF

DO row = g+1, n
  
!        Store off-diagonal column pointer and value
  DO nz = a%offd%beg(row), a%offd%beg(row+1)-1
    IF (a%offd%jco(nz) <= g) THEN
      x21nnz = x21nnz + 1
      x21%jco(x21nnz) = a%offd%jco(nz)
      x21%co(x21nnz)  = a%offd%co(nz)
    ELSE
      x22nnz = x22nnz + 1
      x22%jco(x22nnz) = a%offd%jco(nz) - g
      x22%co(x22nnz)  = a%offd%co(nz)
    END IF
  END DO
  
!        A%diad the non-zero elements in the block diagonal 'a%dia%com' into
!        CSR (sub-)matrix B22:
  
  fircol = (row-1)/BlkSiz*BlkSiz+1
  
  DO col = 1, BlkSiz
    bdval = a%dia%com(row-fircol+1,col+fircol-1)
    IF (bdval /= 0.0D0) THEN
      x22nnz = x22nnz + 1
      x22%jco(x22nnz) = col+fircol-1 - g
      x22%co(x22nnz)  = bdval
    END IF
  END DO
  
!        Store pointer(s) to first non-zero in the next row
  IF (g > 0)  x21%beg(row-g+1) = x21nnz + 1
  x22%beg(row-g+1) = x22nnz + 1
END DO


END SUBROUTINE stomat

END MODULE
