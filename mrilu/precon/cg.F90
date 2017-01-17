!#begindoc

#ifndef WITH_UNION

#define prcmatrix  anymatrix

#endif

#ifdef MKL

#else

#define dot(x,y) DOT_PRODUCT(x,y)

#endif

MODULE m_cg

CONTAINS
 
SUBROUTINE cg (a, x, b)

use m_dump
USE m_build
USE m_chkcnt
USE m_wennz
USE m_glbpars
USE m_solpars
USE m_matvecp
USE m_solve

#ifdef MKL
USE mkl95_precision, ONLY: wp => dp
USE mkl95_blas, ONLY: dot, iamax
#endif

TYPE (prcmatrix)			, POINTER		:: a
DOUBLE PRECISION, DIMENSION(1:a%nschur)	, INTENT(IN OUT)        :: x
DOUBLE PRECISION, DIMENSION(1:a%nschur)	, INTENT(IN OUT)        :: b

!     Implementation of the CG method, as described in
!        "Templates for the Solution of Linear Systems:
!         Building Blocks for Iterative Methods" by Richard Barrett,...,
!     to solve the linear system
!         A x = b

!     The (N-G)x(N-G)-matrix  A  is the Schur-complement of the
!     submatrix  A_11  of the NxN-matrix  A;  A_11 = A(1:G,1:G).
!     The preconditioner, Prc, is stored and can be reached through the descriptor
!     at 'a%mlp'.

!     The maximum number of iterations, MaxNits,
!     the absolute tolerance, AbsTol,
!     and the reduction tolerance, RedTol, are stored in
!     common block /solpars/.

!     The iteration process stops as soon as one of the following
!     conditions holds:
!     . 2-norm residual preconditioned system <= AbsTol for the
!       solution in  x,
!     . 2-norm residual preconditioned system, for the solution in x,
!       has been reduced by a factor <= RedTol, or
!     . the number of iterations >= MaxNits; in this case a warning
!       message is written to standard output.

!     Arguments:
!     ==========
!     a%N   	i   Number of rows/columns in original matrix  A.
!     a%G   	i   Number of rows/columns in left upper part A_11.
!     a%aro 	i   Location of descriptor of the
!                   matrix A, of the linear system.
!     a%mlp 	i   Location of descriptor of the
!                   Preconditioner matrix Prc, stored in MLP format.
!     a%mlp%mlp%perm   i   Preconditioner permutation vector.
!     x        	io  In:  Initial guess for the solution of  A x = b
!                   Out: The solution of the linear system  A x = b_in
!     b        	io  In:  Right hand side of linear system
!                   Out: Residual  b_in - A x_out

!#enddoc

!     Local Variables:
!     ================

!     nflops           Number of floating point operations in
!                      double precision to prevent integer overflow.

INTEGER 					:: i, nza, nzp, ier
INTEGER 					:: nnzd, nnzoff, nnzlas
DOUBLE PRECISION 				:: nflops
DOUBLE PRECISION 				:: nrmz0, nrmz, tol
DOUBLE PRECISION 				:: alpha, rho, rhoold, pAp
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)  	:: p, Ap, z

CHARACTER (LEN=*), PARAMETER :: rounam = 'cg'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


!     Check the preconditioner partition storage and count number of
!     non-zeros:
CALL chkcnt (a%mlp, nnzd, nnzoff, nnzlas)

nflops = 0
nzp    = nnzd + nnzoff + nnzlas

!     Calculate the number of nonzeros in [a%aro]:
nza = wennz (a%aro)
  
ALLOCATE( p(1:a%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( Ap(1:a%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( z(1:a%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

IF (outlev >= 4) PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual', 'Flops/unknown'
  
!     Compute  r^(0) = b - A x^(0)  for the initial guess x^(0) in x:
  
  CALL matvecp (a, x, Ap)

!     {  Ap = q^(0) = A x^(0)  }
  
  b = b - Ap

!     {  b = r^(0) = b - A x^(0)  }
  
  nflops = nflops + 2 * nza + a%nschur
  
!     Solve preconditioned system  Prc z^(0) = r^(0), returning the solution  z^(0)  in  z:

z = b

!     {  z = r^(0)  }

CALL solve (a, z)

!     {  z = z^(0) = inv(Prc) r^(0)  ,  Ap = ?  }

!     Save  ||z^(0)||  in  nrmz0:

nrmz0 = MAXVAL(ABS( z ))
nrmz  = nrmz0

!     {  nrmz = nrmz0 = ||z^(0)||  }

nflops = nflops + 2 * nzp + a%nschur

IF (outlev >= 4) PRINT 9002 , 0, nrmz, nrmz/nrmz0, nflops/DBLE(a%n)

tol = MAX(abstol, nrmz0*redtol)

IF (nrmz0 > tol) THEN

! iteration should start

rho    = dot( b, z )

!     {  rho = rho^(0) = <r^(0) , z^(0)>  }

nflops = nflops + 2 * a%nschur

DO   i = 1, maxnits

! {  i >= 2  ==>  rhoold = rho^(i-2)        ,
! rho = rho^(i-1) = <r^(i-1) , z^(i-1)>  }
  
  IF ( i == 1 ) THEN
!   Initialise search direction:

    p = z

!   {  p^(1) = z^(0)  }
  ELSE
!   Update the search direction:
    
!   {  beta = beta^(i-1) = rho^(i-1) / rho^(i-2)  }
!   {  p^(i) = z^(i-1) + beta^(i-1) * p^(i-1)  }
    
    p = z + (rho / rhoold)*p
    
    nflops = nflops + 2 * a%nschur
  END IF
    
! Compute  Ap := A p:

  CALL matvecp (a, p, Ap)

! {  Ap = q^(i) = A p^(i)  }
  
  nflops = nflops + 2 * nza
  
! Calculate alpha^(i):
  
  pAp = dot( p, Ap )

! {  pAp = <p^(i),q^(i)> = <p^(i), A p^(i)>  }
  
  alpha = rho / pAp

! {  alpha = alpha^(i) = rho^(i-1) / <p^(i), A p^(i)>  }
  
  nflops = nflops + 2 * a%nschur
  
! Update the solution vector  x^(i)  in  x:
  
  x = x + alpha * p

! {  x = x^(i) = x^(i-1) + alpha^(i) * p^(i)  }
  
  nflops = nflops + 2 * a%nschur
  
! Update the residual vector  r^(i)  in  b:
  
  b = b - alpha * Ap

! {  b = r^(i) = r^(i-1) - alpha^(i) * q^(i)  }
  
  nflops = nflops + 2 * a%nschur
  
! Solve preconditioned system  Prc z^(i) = r^(i), returning the solution  z^(i)  in  z:
  
  z = b

! {  z = r^(i)  }
  
  CALL solve (a, z)

! {  z = z^(i) = inv(Prc) r^(i)  ,  Ap = ?  }
  
  nflops = nflops + 2 * nzp
  
  rhoold = rho
  rho    = dot( b, z )

! {  rhoold = rho^(i-1)  ,
! {  rho = rho^(i) = <r^(i) , z^(i)>  }
  
  nflops = nflops + 2 * a%nschur
  
! {  rho = <r^(i),inv(Prc) r^(i)> <= tol**2  }

! Calculate  ||z^(i)||  in  nrmz:

   nrmz = MAXVAL(ABS( z ))

! {  nrmz = ||inv(Prc) r^(i)|| = ||z^(i)||  }
  
  nflops = nflops + a%nschur
  
  IF (outlev >= 4) PRINT 9002 , i, nrmz,  nrmz/nrmz0, nflops/DBLE(a%n)
  
! EXIT LOOP if following condition holds:

  IF ( nrmz <= tol ) EXIT

! Maximum number of iterations reached:

  IF ((i == maxnits) .AND. (outlev >= 1)) THEN
    PRINT '(/, A, 2X, A, A, /, 3X, A, 1P, 2(/, 3X, A, E12.5), /)',  &
      'Warning from subroutine', rounam, '!',  &
      'Residual too large after maximum number of iterations.',  &
      'Inf-norm residual original system:       ',  &
      MAXVAL(ABS( b )),  &
      'Inf-norm residual preconditioned system: ', nrmz
  END IF

END DO

END IF

DEALLOCATE( p, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( Ap, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( z, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')


RETURN

!     Fatal Error, issue error message and return:

!     Format statements:
9000 FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002 FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)

END SUBROUTINE cg

END MODULE