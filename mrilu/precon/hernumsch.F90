!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define csrdmatrix  anymatrix

#endif

MODULE m_hernumsch

CONTAINS

SUBROUTINE hernumsch (g, numrb, a, b )

USE m_build
USE m_dump
 
INTEGER					, INTENT(IN)            :: g
TYPE (csrmatrix)			, POINTER		:: a
INTEGER, DIMENSION(1:a%n)		, INTENT(IN)            :: numrb
TYPE (csrdmatrix)			, POINTER		:: b

!     Computes the renumbered matrix  B  from the off-diagonal
!     matrix  A, with:

!                                 (  O  | B12 )
!        B  =  Prb A inv(Prb)  =  (-----+-----)
!                                 ( B21 | B22 )

!     Prb  is a permutation matrix characterized by the permutation
!     vector 'numRB'


!     Arguments:
!     ==========
!     a%N        i   Number of rows/columns in the matrix A.
!     G        i   Number of rows in upper partition, B11, of matrix B.
!     numRB    i   Permutationvector of the Red-Black permutation, Prb.
!                  numRB(newrow), 1<=newrow<=N, is the old row/column
!                  number corresponding with row/column  newrow  in the
!                  transformed matrix   Prb A inv(Prb).
!     a%beg     i   a%beg(r), 1<=r<=a%N: Index in 'a%jco' and 'a%co' of the
!                  first non-zero off-diagonal element in row  r  of
!                  matrix  A.
!                  a%beg(a%N+1): Index in 'a%jco' and 'a%co' of last nonzero
!                  element + 1.
!     a%jco     i   a%jco(nz), a%beg(r)<=nz<a%beg(r+1), 1<=r<=a%N: Column
!                  number of the non-zero off-diagonal element  a%co(nz)
!                  in row  r  of matrix  A.
!     a%co      i   a%co(a%beg(r):a%beg(r+1)-1), 1<=r<=a%N: Values of the non-
!                  zero off-diagonal elements in row  r  of matrix  A.
!     b%offd%beg     o   b%offd%beg(r), 1<=r<=a%N: Index in 'b%offd%jco' and 'b%offd%co' of the
!                  first non-zero off-diagonal element in row  r  of
!                  matrix  B.
!                  b%offd%beg(a%N+1): Index in 'b%offd%jco' and 'b%offd%co' of last nonzero
!                  element + 1.
!     b%offd%jco     o   b%offd%jco(nz), b%offd%beg(r)<=nz<b%offd%beg(r+1), 1<=r<=N: Column
!                  number of the non-zero off-diagonal element  b%offd%co(nz)
!                  in row  r  of matrix  B.
!     b%offd%co      o   b%offd%co(b%offd%beg(r):b%offd%beg(r+1)-1), 1<=r<=a%N: Values of the non-
!                  zero off-diagonal elements in row  r  of matrix  B.
!     B%lotr   o   B%lotr(r), 1<=r<=G, index in 'b%offd%jco' and 'b%offd%co' of
!                  last nonzero in 'r'-th row of B12.
!                  B%lotr(r), G+1<=r<=a%N, index in 'b%offd%jco' and 'b%offd%co' of
!                  last nonzero in 'r'-th row of B21.

!#enddoc

!     Local Variables:
!     ================

INTEGER					:: ier
INTEGER 				:: newrow, newcol, oldrow, lasnzb, nzb, nz, bnnz
INTEGER, ALLOCATABLE, DIMENSION(:)	:: invper

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'hernumsch'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( invper(1:a%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Construct permutationvector for permutation inv(Prb) in 'InvPer':
FORALL (newrow = 1:a%n) invper(numrb(newrow)) = newrow

bnnz  = 0
b%offd%beg(1) = 1

!     The upper partition:
DO newrow = 1, g
  oldrow = numrb(newrow)
  
!        Store non-zero off-diagonal elements in 'OldRow' of  A into
!        'NewRow' of  B:
  DO nz = a%beg(oldrow), a%beg(oldrow+1)-1
    bnnz       = bnnz + 1
    b%offd%jco(bnnz) = invper(a%jco(nz))
    b%offd%co(bnnz)  = a%co(nz)
  END DO
  b%lotr(newrow) = bnnz
  b%offd%beg(newrow+1) = bnnz + 1
END DO

!     The lower partition:
DO newrow = g+1, a%n
  oldrow = numrb(newrow)
  
  lasnzb = bnnz + a%beg(oldrow+1) - a%beg(oldrow)
  nzb    = lasnzb
  
!        Store non-zero off-diagonal elements in 'OldRow' of  A into
!        'NewRow' of  B:
  DO nz = a%beg(oldrow), a%beg(oldrow+1)-1
    newcol = invper(a%jco(nz))
    
    IF (newcol <= g) THEN
!              Element in block  B21, store in first part of 'NewRow':
      bnnz       = bnnz + 1
      b%offd%jco(bnnz) = newcol
      b%offd%co(bnnz)  = a%co(nz)
    ELSE
!              Element in block  B22, store in last part of 'NewRow':
      b%offd%jco(nzb) = newcol
      b%offd%co(nzb)  = a%co(nz)
      nzb       = nzb - 1
    END IF
  END DO
  
  b%lotr(newrow) = bnnz
  bnnz         = lasnzb
  b%offd%beg(newrow+1) = bnnz + 1
END DO

DEALLOCATE( invper, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')


END SUBROUTINE hernumsch

END MODULE