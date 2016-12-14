!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_permrc

CONTAINS
 
SUBROUTINE permrc (a, perm)

USE m_build

TYPE (csrmatrix)		, POINTER		:: a
INTEGER, DIMENSION(1:)		, INTENT(IN)            :: perm

!     Permute Row or Column numbers.
!     Reorder the column(/row) numbers of a CSR(/CSC) stored matrix
!     in a%jco according the permutation given in Perm.

!     Arguments:
!     ==========
!     A%N       i   Number of rows(/columns) in CSR(/CSC) matrix A.
!     a%beg	i   a%beg(r): Location in 'a%jco' and 'coA' of first
!                   non-zero off-diagonal element in each row/column
!                   r of matrix  A.
!     a%jco     i   a%jco(a%beg(r):a%beg(r+1)-1): Column numbers of the
!                   non-zero off-diagonal elements in row r of matrix A.
!     a%jco     io  Input:  a%jco(a%beg(r):a%beg(r+1)-1): Column/Row numbers
!                           of the non-zero off-diagonal elements in
!                           row/column r of matrix A.
!                   Output: a%jco(a%beg(r):a%beg(r+1)-1): Permuted
!                           Column/Row numbers of the non-zero
!                           off-diagonal elements in row/column r of
!                           matrix A.
!     Perm     	i   Permutation vector.

!#enddoc

INTEGER :: nz

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'permrc'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

FORALL (nz  = a%beg(1):a%beg(a%n+1) - 1)  a%jco(nz) = perm(a%jco(nz))

!     End of  permrc
END SUBROUTINE permrc

END MODULE