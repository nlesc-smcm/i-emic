!#begindoc
 
#ifndef WITH_UNION

#define prcmatrix  anymatrix
#define scbmmatrix  anymatrix

#endif

MODULE m_solprc

CONTAINS

SUBROUTINE solprc (a, x, b)

USE m_dump
USE m_build
USE m_possred
USE m_chkcnt
USE m_cmpsol
USE m_glbpars
USE m_presred
USE m_dperv

TYPE (prcmatrix)			, POINTER		:: a
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN OUT)     	:: x
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(IN)     	:: b

!     Solves  the linear system
!        Aorig x = b
!     using a CG-like iterative method with a preconditioner.

!     The scaled and Red-Black reordered version,
!        A (= Prb Scale Aorig inv(Prb))
!     of the original matrix  Aorig, is partitioned as:
!            ( A_11 | A_12 )
!        A = (------+------)
!            ( A_21 | A_22 )
!     where  A_11  is a non-singular block-diagonal sub-matrix of  A.

!     The preconditioner, Prc, is stored and can be reached through the descriptor 
!     at 'a'

!     The inverse of the (block-)diagonal submatrix, A_11, is stored
!     and can be reached through the descriptor at 'a'.
!     The rest of the matrix  A  is stored in [begA,jcoA,coA].


!     Arguments:
!     ==========
!     A%N        i   Number of rows/columns in the matrix  A  or  Aorig.
!     A%BlkSiz   i   Number of rows/columns in a diagonal block.
!                  The value of 'A%N' should be an integer multiple of
!                  'A%BlkSiz'.
!     a    i   Location of the
!                  complete preconditioner of matrix A, computed in a
!                  previous call of subroutine 'cmpprc'.
!     x        io  In:  Initial guess for the solution of  Aorig x = b
!                  Out: The solution of the linear system  Aorig x = b
!     b        i   In:  Right-hand side of the matrix equation
!                       Aorig x = b

!     It is assumed that the global variables in the common blocks have
!     been initialised by calling the subroutine 'inisol'.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'solprc'

!     Local Variables:
!     ================

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: rhs, rhs2, x2
TYPE (scbmmatrix), POINTER			:: Aaro
INTEGER 					:: nnzd, nnzoff, nnzlas, ier
DOUBLE PRECISION 				:: begtim, endtim


#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Compute the solution in 'x':
!     ----------------------------

CALL CPU_TIME(begtim)

ALLOCATE( rhs(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Scale the right-hand side vector, b = Scale b:

rhs = A%scale * b

!     Permute the right-hand side vector, b = Prb Scale b:
CALL dperv (.false., A%perrb, rhs, rhs)

!     Permute the initial starting vector, x = Prb x:
CALL dperv (.false., A%perrb, x, x)

IF (outlev >= 4) THEN
!        Check the preconditioner partition storage and count number
!        of non-zeros:
  CALL chkcnt (A%mlp, nnzd, nnzoff, nnzlas)
  
  PRINT '(2(/, A, 1X, I9), /, A, 4X, F9.2 )' ,  &
      'Number of non-zeros in the LDU-decomposition:', nnzd + nnzoff + nnzlas,  &
      'Number of non-zeros in last block:           ', nnzlas,  &
      'Average number of non-zeros per row in L+D+U:',  &
      DBLE(nnzd+nnzoff+nnzlas)/DBLE(A%nschur)
END IF

IF (A%g > 0) THEN
!        Solve reduced problem first.
  
!        Request workspace for right hand side, rhs2, and solution
!        vector, x2, of reduced system:

  ALLOCATE (rhs2(1:A%nschur), STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
  
  ALLOCATE (x2(1:A%nschur), STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!        Extract location of the submatrices of Aro:
  Aaro => anytoscbm(A%aro)
  
!        Construct the right-hand side 'rhs2' and starting vector 'x2'
!        of the reduced problem
!           S x2 = rhs2
!        with
!           rhs2 = b_2 - A_21 (1/A_11) b_1
!           x2   = x_2
  CALL presred (Aaro, rhs, x, rhs2, x2)

!        Solve the reduced problem.
  
!        Compute the solution of  S x2 = rhs  and store the result
!        into 'x2':
  CALL cmpsol (a, rhs2, x2)
  
!        Compute the solution of the complete problem into
!        x = (x_1 | x_2)', with:
!        x_1 := (1/A_11) ( rhs1 - A_12 x_2 )  and
!        x_2 := x2
  CALL possred (Aaro, x2, rhs, x)
!        Free workspace of DP Temporary arrays:

  DEALLOCATE (rhs2, STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

  DEALLOCATE (x2, STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

ELSE
!        No reduced problem, Solve complete problem.
  
!        Compute the solution of  A x = b  and store the result in 'x':
  CALL cmpsol (a, rhs, x)
  
END IF

DEALLOCATE( rhs, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

!     Permute back the solution vector, x = inv(Prb) x:
CALL dperv (.true., A%perrb, x, x)

CALL CPU_TIME(endtim)

IF (outlev >= 3) PRINT '(A, X, F7.2, /)' , 'Time Solution:', endtim - begtim

RETURN

END SUBROUTINE solprc

END MODULE
