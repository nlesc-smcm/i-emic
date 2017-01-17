!#begindoc

#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define prcmatrix  anymatrix
#define scdematrix anymatrix

#endif

#ifdef MKL

#define norm(x) nrm2(x)

#else

#define dot(x,y) DOT_PRODUCT(x,y)
#define norm(x) SQRT(SUM(x**2))

#endif

MODULE m_gmres

CONTAINS
 
SUBROUTINE gmres (A, x, b)

USE m_dump
USE m_build
USE m_chkcnt
USE m_wennz
USE m_glbpars
USE m_solpars
USE m_matvecp
USE m_solve

#ifdef MKL

USE mkl95_precision, ONLY: wp => dp
USE mkl95_blas, ONLY: dot, nrm2, gemv

#endif

TYPE (prcmatrix)			, POINTER		:: a
DOUBLE PRECISION, DIMENSION(1:A%nschur)	, INTENT(IN OUT)        :: x
DOUBLE PRECISION, DIMENSION(1:A%nschur)	, INTENT(IN)            :: b

!     Implementation of the GMRES(M) method, as described in
!        "Templates for the Solution of Linear Systems:
!         Building Blocks for Iterative Methods" by Richard Barrett,...,
!     to solve the linear system
!         A x = b

!     The (N-A%G)x(N-A%G)-matrix  A  is the Schur-complement of the
!     submatrix  A_11  of the A%NxA%N-matrix  A;  A_11 = A(1:A%G,1:A%G).
!     The preconditioner, Prc, is stored and can be reached through the descriptor
!     at 'A%mlp'.

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

!     M        	i   Dimension of the subspace of R^(A%NSCHUR),  M <= A%NSCHUR.
!     A%N       i   Number of rows/columns in original matrix  A.
!     A%G       i   Number of rows/columns in left upper part A_11.
!     A%nschur  i   Size of 1st Schur-complement, S, of A_11 in A.
!     A%aro     i   Location of descriptor of the
!                   matrix A, of the linear system.
!     A%mlp    	i   Location of descriptor of the
!                   Preconditioner matrix Prc, stored in MLP format.
!     x        	io  In:  Initial guess for the solution of  A x = b
!                   Out: The solution of the linear system  A x = b
!     b        	i   Right-hand side of the matrix equation  A x = b
!                         the solution as indicated above.

!#enddoc

!     Global Parameters:
!     ==================

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'gmres'

!     Local Variables:
!     ================

!     nflops           Number of floating point operations in
!                      double precision to prevent integer overflow.
!     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
!     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
!     tol              = MAX (AbsTol, ||r^(0)||*RedTol)

INTEGER 					:: nschur, i, j, k, m, iter, ier
INTEGER 					:: nnzd, nnzoff, nnzlas, nza, nzp
DOUBLE PRECISION 				:: nflops
DOUBLE PRECISION 				:: r, factor, hki, nrmr, nrmr0, tol
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: w, jc, js, s, y
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)	:: v, h

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

nschur = A%nschur
m      = MIN(mgmres, nschur)

ALLOCATE( w(1:nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( v(1:nschur,1:m+1), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( h(1:m+1,1:m+1), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( jc(1:m), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( js(1:m), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( s(1:m+1), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( y(1:m+1), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Check the preconditioner partition storage and count number of
!     non-zeros:
CALL chkcnt (A%mlp, nnzd, nnzoff, nnzlas)

nflops = 0
nzp    = nnzd + nnzoff + nnzlas

!     Calculate the number of nonzeros in [A%aro]:
nza = wennz (A%aro)

IF (outlev >= 4) THEN
  PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual', 'Flops/unknown'
END IF

factor = 1.0D0
nrmr=tol+1

!     Starting point new iteration of GMRES(M)
iter = 0
DO WHILE ((iter < maxnits) .and. (nrmr > tol))
! iter=1,maxnits
  
!        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
!        for the initial guess   x^(iter)  in  x:
  
  CALL matvecp (A, x, w)

!        {  w = A x^(iter)  }
  
!        Calculate residual
  w = b - w
!        {  w = b - A x^(iter)  }
  
  CALL solve (A, w)

!        {  w = r^(iter) = inv(Prc) (b - A x^(iter)) }
  
  nflops = nflops + 2 * nza + nschur + 2 * nzp
  
  
!        Calculate 2-norm of preconditioned residual:

  nrmr  = norm( w )

!        {  nrmr = ||w|| = ||r^(iter)||  }
  
  nflops = nflops + 2 * nschur
  
  
  IF (iter == 0) THEN

    nrmr0 = nrmr

!   {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
    
    IF (outlev >= 4) THEN
      PRINT 9002 , iter, nrmr, nrmr/nrmr0, nflops/DBLE(A%n)
    END IF

    tol = MAX(abstol, nrmr0*redtol)

  END IF
  
  
  IF (nrmr > tol)  THEN
  
  
!        Initialise the Basis for the Generalised Minimal RESidual
!        Method:
  
  s(1) = nrmr
  s(2:m+1) = 0.0D0

! {  s = ||r^(0)|| e_1  }
  
  
  v(:,1) = w *  1.0D0 / nrmr

!        {  V(:,1) = w / ||w||  }
  
  nflops = nflops + nschur
  
  DO  i = 1, MIN(m,maxnits-iter+1)
    
!   Solve  w  from  Prc w = A V(:,i):
    
    CALL matvecp (A, v(:,i), w)

!   {  w = A V(:,i)  }
    
    CALL solve (A, w)

!   {  w = r^(i) = inv(Prc) A v^(i))  }
    
    nflops = nflops + 2 * nza + 2 * nzp

!   Start the modified Gram-Schmidt process:

#   ifdef MKL
      CALL gemv( v(:,1:i), w, h(1:i,i), TRANS='T'  )
#   else
      h(1:i,i) = MATMUL( w, v(:,1:i) )
#   endif

!   {  H(k,i) = <w,V(:,k)>  }

#   ifdef MKL
      CALL gemv( v(:, 1:i), h(1:i,i), w, ALPHA = -1.0D0, BETA = 1.0D0 )
#   else
      w = w - MATMUL( v(:, 1:i), h(1:i,i) )
#   endif

!   {  w := w - H(k,i) V(:,k)  }


    nflops = nflops + i * 4 * nschur
    
    
    h(i+1,i) = norm( w )
    
    v(:,i+1) = w * 1.0D0 / h(i+1,i)

!   {  H(i+1,i) = ||w||  ,  V(:,i+1) = w / H(i+1,i)  }
    
    nrmr = h(i+1,i)
    
    nflops = nflops + 3 * nschur
    
!   Apply the Givens rotations J_1, J_2, ..., J_(i-1) on  H(1:i,i):

    DO k = 1, i-1
      hki      = h(k,i)
      h(k  ,i) = jc(k)*hki - js(k)*h(k+1,i)
      h(k+1,i) = js(k)*hki + jc(k)*h(k+1,i)
    END DO
    
!   Construct Givens rotation  J_i  acting on H(i:i+1,i), such that (i+1)st component of J_i H(:,i) is 0:

    r = SQRT (h(i,i)**2 + h(i+1,i)**2)
    jc(i) =  h(i  ,i) / r
    js(i) = -h(i+1,i) / r
    
!           Apply  J_i  on  H(i:i+1,i):
    h(i  ,i) = jc(i) * h(i,i) - js(i) * h(i+1,i)
    h(i+1,i) = 0
    
!           Apply  J_i  on  s(i:i+1):
    r     = s(i)
    s(i  ) = jc(i) * r - js(i) * s(i+1)
    s(i+1) = js(i) * r + jc(i) * s(i+1)
    
    nflops = nflops + (i-1) * 6 + 15
    
    factor = factor * ABS(js(i))
    nrmr   = factor * nrmr0
    iter = iter + 1
    IF (nrmr <= tol  .Or. i == m  .Or. iter == maxnits) THEN

!     Compute the solution in  x:
      
!     Compute  y  as the solution of  H y = s(1:i), in which
!     the upper  i x i  triangular part of  H  has  H(j,k) as its elements.

      DO j = i, 1, -1

        y(j) = ( s(j) - dot( h(j,j+1:i), y(j+1:i) ) ) / h(j,j)
        
!       Update the solution in  x:

        x = x + v(:,j)* y(j)

      END DO

!     {  x = x^(0) + V(:,j:i) y(j:i)  }
      
      nflops = nflops + i*(i-1) + i * (1 + 2 * nschur)
    END IF
    
    IF (outlev >= 4) THEN
      PRINT 9002 , iter, nrmr, nrmr/nrmr0, nflops/DBLE(A%n)
    END IF
    
    IF (nrmr <= tol) EXIT

  END DO
  END IF
!   Maximum number of iterations reached?:

  IF ( (iter == maxnits) .AND. (outlev >= 1) ) THEN
    PRINT '(/, A, 2X, A, A, /, 3X, A, 1P, (/, 3X, A, E12.5), /)',  &
          'Warning from subroutine', rounam, '!',  &
          'Residual too large after maximum number of iterations.',  &
          '2-norm residual preconditioned system: ', nrmr
  END IF  

END DO    

  DEALLOCATE( w, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  DEALLOCATE( v, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  DEALLOCATE( h, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  DEALLOCATE( jc, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  DEALLOCATE( js, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  DEALLOCATE( s, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  DEALLOCATE( y, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')


    RETURN
    
    
!     Format statements:
    9000 FORMAT (/, a6, 2(2X, a18), 2X, a14)
    9002 FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
    
  END SUBROUTINE gmres

END MODULE
