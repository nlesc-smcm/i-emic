!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_eblkdia

CONTAINS

SUBROUTINE eblkdia (a, Ad)

USE m_glbpars
USE m_build

TYPE (csrmatrix)				, POINTER		:: a
TYPE (diamatrix)				, POINTER		:: Ad

!     Extracts the a block diagonal from the sparse matrix A stored in CSR format.
!     The extracted block diagonal matrix is stored in  [Ad%com]
!     and the off-diagonal non-zero elements of matrix A are stored
!     back in [A%beg,A%jco,A%co].


!     Arguments:
!     ==========
!     A%N	i   Number of rows/columns in the matrix A.
!     Ad%blksiz	i   Number of rows/columns in a diagonal block.
!                   'A%N' should be an integer multiple of 'Ad%blksiz'!
!     A%beg     io  In:  A%beg(i): index in 'A%jco' and 'A%co' of the
!                        first non-zero element in row i of matrix A.
!                   Out: A%beg(i): index in 'A%jco' and 'A%co' of the
!                        first off-block-diagonal element in row i of matrix A.
!     A%jco     io  In:  A%jco(nz): column number of non-zero element A%co(nz).
!                   Out: A%jco(nz): column number of off-block-diagonal
!                        non-zero element A%co(nz).
!     A%co      io  In:  A%co(nz): value of non-zero element of  A.
!                   Out: A%co(nz): value of off-block-diagonal non-zero element.
!     Ad%com     o   The values of the elements in the diagonal blocks of
!                   matrix A, stored in column major order.

!#enddoc

!     Local Variables:
!     ================

INTEGER :: row, start, nza, nnza
INTEGER :: i, firrow

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'eblkdia'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Initialise all diagonal blocks with 0:
Ad%com = 0.0D0

nnza  = 0

DO firrow = 0, A%n-1, Ad%blksiz
  
  DO i = 1, Ad%blksiz
    row        = firrow + i
    start      = A%beg(row)
    A%beg(row) = nnza + 1
    
!           Put default diagonal element in diagonal block:
    Ad%com(i,row) = neglgbl
    
!           Extract diagonal block from A and shift the off-diagonal
!           elements in A:
    DO nza = start, A%beg(row+1)-1
      IF( (A%jco(nza) > firrow) .AND. (A%jco(nza) <= firrow + Ad%blksiz) ) THEN
!       Part of the block diagonal, store in block diagonal:
        Ad%com(i,A%jco(nza)) = A%co(nza)
      ELSE
!       Shift off-block-diagonal value:
        nnza       = nnza + 1
        A%co (nnza) = A%co(nza)
        A%jco(nnza) = A%jco(nza)
      END IF
  END DO
END DO
END DO
A%beg(A%n+1) = nnza + 1

END SUBROUTINE eblkdia

END MODULE