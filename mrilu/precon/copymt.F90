!#begindoc
 
#ifndef WITH_UNION

#define scdematrix  anymatrix

#endif

MODULE m_copymt

CONTAINS

SUBROUTINE copymt (A, DropTol, B)

USE m_build
USE m_wascde
USE m_wcompr
USE m_glbpars

TYPE (scdematrix)	, POINTER		:: A
DOUBLE PRECISION	, INTENT(IN)         	:: DropTol
TYPE (scdematrix)	, POINTER		:: B

!     Copy the matrix [A] into a new matrix [B].  Only the
!     (off-diagonal) elements of [A] which are greater than 'DropTol' are copied into matrix [B].

!     Arguments:
!     ==========
!     A      	i   Location of descriptor of the matrix A.
!     DropTol  	i   The (off-diagonal) elements of matrix [A], of which
!                   the absolute value is greater than 'DropTol, are copied into [B].
!     B      	o   Location of descriptor of the matrix B.

!#enddoc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'copymt'

!     Local variables:
!     ================

INTEGER 					:: i, j, nnza, nnzb
DOUBLE PRECISION				:: val

  nnza = csrnnz(A%offd)
  CALL wascde (A%n, A%dia%blksiz, nnza, B)
  
  IF (DropTol == 0.0D0) THEN
    B%offd%beg = A%offd%beg
    B%offd%jco = A%offd%jco
    B%offd%co  = A%offd%co
  ELSE
    nnzb          = 0
    B%offd%beg(1) = 1
    DO i = 1, A%offd%n
      DO j = A%offd%beg(i), A%offd%beg(i+1)-1
        val = A%offd%co(j)
        IF ( ABS(val) > DropTol) THEN
          nnzb             = nnzb + 1
          B%offd%co(nnzb)  = val 
          B%offd%jco(nnzb) = A%offd% jco(j)
        END IF
      END DO
      B%offd%beg(i+1) = nnzb + 1
    END DO
    IF (nnzb < nnza) THEN
      IF (outlev >= 5) THEN
        PRINT '(/, A, X, A, /, 2(3X, A, I11, /))' , 'In', rounam,  &
              'Nr nonzeros in A (org):', nnza, 'Nr nonzeros in B (new):', nnzb
      END IF
    END IF
  END IF
    
  B%dia%com = A%dia%com

!   Compress the storage by shifting the representation of 'B' and adjust the required size:

  CALL wcompr (scdetoany(B))
  
END SUBROUTINE copymt

END MODULE