!#begindoc

#ifndef WITH_UNION

#define prcmatrix  anymatrix

#endif

MODULE m_applmlp
 
CONTAINS

SUBROUTINE applmlp (a, x)

USE m_build
USE m_solve
USE m_dperv

TYPE (prcmatrix)			, POINTER 		:: a
DOUBLE PRECISION, DIMENSION(1:a%nschur)	, INTENT(IN OUT)        :: x

!     CoMPute the SOLution of the linear equation
!        S x = b
!     using the Multi Level Preconditioner [a%mlp].

!     The matrix 'S' is the Schur-complement
!        S = A_22 - A_21 (1/A_11) A_1
!     of the sub-matrix 'A_11' in the matrix 'A' of the original
!     linear system
!        A x = b

!     The original equation can be partitioned as follows:
!        ( A_11 | A_12 ) ( x_1 )     ( b_1 )
!        (------+------) (-----)  =  (-----)
!        ( A_21 | A_22 ) ( x_2 )     ( b_2 )
!     ==
!        x_1                              =  (1/A_11) (b_1 - A_12 x_2)   ,
!        (A_22 - A_21 (1/A_11) A_12) x_2  =  (b_2 - A_21 (1/A_11) b_1)
!     ==
!        x_1    =  (1/A_11) (b_1 - A_12 x_2)   ,
!        S x  =  b
!     with
!        A_11   is a A%G x A%G block-diagonal sub-matrix of matrix A.
!        b  = P (b_2 - A_21 (1/A_11) b_1)
!        S    = A_22 - A_21 (1/A_11) A_12
!               is the Schur-complement of  A_11  in  A.

!     Arguments:
!     ==========
!     A%N       i   Number of rows/columns in matrix  A.
!     A%G       i   Number of rows/columns in left upper partition A_11.
!     a%aro 	i   Location of the descriptor of the matrix of the
!                   linear system to be solved.
!     a%mlp 	i   Location of the descriptor of the preconditioner Prc.
!     b      	io  In:  Right Hand Side of the equation.
!                   Out: --- undefined.
!     x        	io  In:  Starting value for solution of the equation.
!                   Out: Solution of the equation  A x = b.

!#enddoc

!     Local Variables:
!     ================

CHARACTER (LEN=*), PARAMETER :: rounam = 'applmlp'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Because
!        S x = b
!     ==
!        P S (1/P) P x = P b
!     and the preconditioner for  P S (1/P)  is stored in [a%mlp], we
!     will solve (P x) from
!        (P S (1/P)) (P x) = (P b)


! Permute the starting value for the solution  x := P x:

  CALL dperv (.false.,  a%mlp%perm, x, x)

! Apply preconditioner to x, i.e. solve a system using mlp

  CALL solve (a, x)
      
! Unpermute the solution  x := inv(P) x:
  
  CALL dperv (.true., a%mlp%perm, x, x)
      
  END SUBROUTINE applmlp

END MODULE
