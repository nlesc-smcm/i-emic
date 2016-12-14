!#begindoc
 
#ifndef WITH_UNION

#define prcmatrix  anymatrix
#define scbmmatrix  anymatrix

#endif

MODULE m_applprc

CONTAINS

SUBROUTINE applprc (a, x, b)

USE m_dump
USE m_build
USE m_possred
USE m_chkcnt
USE m_applmlp
USE m_glbpars
USE m_presred
USE m_dperv

TYPE (prcmatrix)			, POINTER		:: a
DOUBLE PRECISION, DIMENSION(1:A%n)	, INTENT(OUT)     	:: x
DOUBLE PRECISION, DIMENSION(1:A%n)    	, INTENT(IN)        	:: b

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


!     Arguments:
!     ==========
!     A%N	i   Number of rows/columns in the matrix  A  or  Aorig.
!     A%blksiz	i   Number of rows/columns in a diagonal block.
!                   The value of 'A%N' should be an integer multiple of
!                   'A%blksiz'.
!     a        	i   Location of the
!                   complete preconditioner of matrix A, computed in a
!                   previous call of subroutine 'cmpprc'.
!     x        	o   Out: The solution of the linear system  Aorig x = b
!     b        	i   In:  Right-hand side of the matrix equation
!                       Aorig x = b

!     It is assumed that the global variables in the common blocks have
!     been initialised by calling the subroutine 'inisol'.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'applprc'

!     Local Variables:
!     ================

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: rhs, x2
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

!     Scale the input vector, x = Scale b:

x = A%scale * b


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
  
!        Request workspace for x2
  
  ALLOCATE( rhs(1:A%n), STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

  ALLOCATE (x2(1:A%nschur), STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

  Aaro => anytoscbm(A%aro)

!        Construct the right-hand side 'rhs2' and starting vector 'x2'
!        of the reduced problem
!           S x2 = rhs2
!        with
!           rhs2 = b_2 - A_21 (1/A_11) b_1
!           x2   = x_2
  CALL presred (Aaro, x, x, x2, x2)

!        Solve the reduced problem.
  
!        Compute the solution of  S x2 = x2  and store the result
!        into 'x2':
  CALL applmlp (a, x2)
  
!        Compute the solution of the complete problem into
!        x = (x_1 | x_2)', with:
!        x_1 := (1/A_11) ( rhs1 - A_12 x_2 )  and
!        x_2 := x2
  CALL possred (Aaro, x2, x, x)
!        Free workspace of DP Temporary arrays:

  DEALLOCATE (x2, STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

  DEALLOCATE (rhs, STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

ELSE
!        No reduced problem, Solve complete problem.
  
!        Compute the solution of  A x = b  and store the result in 'x':
!  x=b
  CALL applmlp (a, x)
  
END IF

!     Permute back the solution vector, x = inv(Prb) x:
CALL dperv (.true., A%perrb, x, x)

CALL CPU_TIME(endtim)

IF (outlev >= 3) THEN
  PRINT '(A, X, F7.2, /)' , 'Time Solution:', endtim - begtim
END IF


END SUBROUTINE applprc

END MODULE
