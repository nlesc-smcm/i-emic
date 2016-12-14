!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix
#define diamatrix  anymatrix
#define partmatrix anymatrix
#define mlpmatrix  anymatrix

#endif

MODULE m_chkcnt

CONTAINS

SUBROUTINE chkcnt (xprc, nnzd, nnzoff, nnzlas)

USE m_dump
USE m_build
USE m_mterrmsg

TYPE (mlpmatrix), POINTER                :: xprc
INTEGER, INTENT(OUT)                     :: nnzd
INTEGER, INTENT(OUT)                     :: nnzoff
INTEGER, INTENT(OUT)                     :: nnzlas

!     Check the structure of the MLP preconditioner matrix  P,
!     descriptor xprc, and count the number of nonzeros in  P.

!     Arguments:
!     ==========
!     xprc    i   Location of Matrix descriptor.
!     NnzD     o   Number of non-zeros in the diagonal blocks of all
!                  partions except the last.
!     NnzOff   o   Number of non-zeros in the off-diagonal parts of all
!                  partions except the last.
!     NnzLas   o   Number of non-zeros in the last partition.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'chkcnt'

!     Local Variables:
!     ================

INTEGER 					:: npar
TYPE (csrmatrix), POINTER			:: ixupp
TYPE (cscmatrix), POINTER			:: ixlow
TYPE (partmatrix), POINTER			:: xpar
TYPE (diamatrix), POINTER  			:: xd

!     Initialise counts of non-zero elements:
nnzd   = 0
nnzoff = 0
nnzlas = 0

!     Actual partition is first one:
xpar => xprc%first

!     Loop over linked-list to consider each partition, but the last:
DO WHILE ( ASSOCIATED( xpar) )
!DO WHILE ( .NOT. ASSOCIATED( xpar, xprc%last) )

  SELECT CASE (xpar%typ)
    CASE (pldutp)

!     Partition is of type PLDU
  
!     Extract the matrix pointer and check the storage type
!     (diagonal matrix expected):

      xd => xpar%dia
  
!     Increment number of non-zeros:
  
      nnzd = nnzd + xd%n * xd%blksiz
   
! ###    Lower triangular part of partition:
!        Lower triangular part stored in CSC format
 
      ixlow => xpar%ltr
  
!     Check compatibility of diagonal and off-diagonal blocks

      IF (xd%n /= ixlow%n) CALL dump(__FILE__,__LINE__, 'Incompatible dimensions in partition!')
    
!     Increment number of non-zeros:
 
      nnzoff = nnzoff + ixlow%nnz
  
! ###    Upper triangular part of partition:
!        Upper triangular part stored in CSR format
 
      ixupp => xpar%utr
  
!     Extract and check the matrix pointers

!     Check compatibility of diagonal and off-diagonal blocks:

      IF (xd%n /= ixupp%n) CALL dump(__FILE__,__LINE__, 'Incompatible dimensions in partition!')
    
!           Increment number of non-zeros:

      nnzoff = nnzoff + ixupp%nnz

    CASE (pffptp)

!     Last partition is a Full Factored partition with Pivots:
  
!     Descriptor Full Matrix:
!     Matrix Size:
!     Number of elements in last partition:

       nnzlas = nnzlas + xprc%last%fm%n**2

    CASE (psfptp) 

!     Last partition is a Sparse Factored partition with Pivots:
  
!     Diagonal part of partition:
!     Extract the matrix pointer and check the storage type
!     (diagonal matrix expected):

      xd => xprc%last%dia
  
!     Off-diagonal matrix:
!     Check compatibility of diagonal and off-diagonal dimensions

      IF (xd%n /= xprc%last%offd%n) CALL dump(__FILE__,__LINE__, 'Incompatible dimensions in partition!')
  
!     Number of non-zeros in last partition:

      nnzlas = nnzlas + xd%n*xd%blksiz + xprc%last%offd%nnz

    CASE DEFAULT

!     Issue error message and return (ier < 0):
      CALL mterrmsg (rounam, xprc%typ, pldutp)
      CALL dump(__FILE__,__LINE__,'Illegal matrix type')

    END SELECT
  
!        Next partition
  xpar => xpar%next

END DO

END SUBROUTINE chkcnt

END MODULE
