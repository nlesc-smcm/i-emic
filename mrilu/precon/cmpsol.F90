!#begindoc

#ifndef WITH_UNION

#define prcmatrix  anymatrix

#endif

MODULE m_cmpsol
 
CONTAINS

SUBROUTINE cmpsol (a, b, x)

USE m_build
USE m_dump
USE m_bicgstab
USE m_bicgstabr
USE m_cg
USE m_solpars
USE m_gmres
USE m_dperv

TYPE (prcmatrix)			, POINTER 		:: a
DOUBLE PRECISION, DIMENSION(1:a%nschur)	, INTENT(IN OUT)        :: b
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

CHARACTER (LEN=*), PARAMETER :: rounam = 'cmpsol'

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

      IF (a%nschur == 0) RETURN

!     Permute the right hand side  b := P b:
      CALL dperv (.false., a%mlp%perm, b, b)

!     Permute the starting value for the solution  x := P x:
      CALL dperv (.false., a%mlp%perm, x, x)

      IF (cgtype == 0) THEN
!        Reduced system to be solved by Conjugate Gradient method:
  
!        Solve the reduced system: 
         CALL cg (a, x, b)
  
      ELSE IF (cgtype == 1) THEN
!        Reduced system to be solved by Bi-CGSTAB method:
    
!        Solve the reduced system:
         CALL bicgstab (a, x, b)
    
      ELSE IF (cgtype == 2) THEN
!        Reduced system to be solved by Bi-CGSTAB method:
      
!        Solve the reduced system:
         CALL bicgstabr (a, x, b)
      
      ELSE IF (cgtype == 3) THEN
!        Reduced system to be solve by GMRES(M) method:
        
!        Solve the reduced system:
         CALL gmres (a, x, b)
        
      ELSE
         PRINT '(A, 2X, A, /, 3X, A, I11)' , 'Internal error in', rounam, 'Unknown solution method, CGType =', cgtype
         CALL dump(__FILE__,__LINE__,'Unknown solution method')
      END IF
      
!     Unpermute the solution  x := inv(P) x:
      CALL dperv (.true., a%mlp%perm, x, x)

  END SUBROUTINE cmpsol

END MODULE