!#begindoc

#ifndef WITH_UNION

#define prcmatrix  anymatrix

#endif

#ifdef MKL

#define norm(x) nrm2(x)

#else

#define dot(x,y) DOT_PRODUCT(x,y)
#define norm(x) SQRT(SUM(x**2))

#endif

MODULE m_bicgstabr

CONTAINS
 
SUBROUTINE bicgstabr (a, x, r)

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
USE mkl95_blas, ONLY: dot, nrm2
#endif

TYPE (prcmatrix)			, POINTER		:: a
DOUBLE PRECISION, DIMENSION(1:a%nschur)	, INTENT(IN OUT)	:: x
DOUBLE PRECISION, DIMENSION(1:a%nschur)	, INTENT(IN OUT)        :: r

!     Implementation of the Bi-CGSTAB method, using a Right
!     preconditioner, as described in
!        "Templates for the Solution of Linear Systems:
!         Building Blocks for Iterative Methods" by Richard Barrett,...,
!     to solve the linear system
!         A x = b

!     The (N-A%G)x(N-A%G)-matrix  A  is the Schur-complement of the
!     submatrix  A_11  of the NxN-matrix  A;  A_11 = A(1:A%G,1:A%G).
!     The preconditioner, Prc, is stored and can be reached through the descriptor
!     at 'a%mlp'.

!     The maximum number of iterations, MaxNits,
!     the absolute tolerance, AbsTol,
!     and the reduction tolerance, RedTol, are stored in
!     common block /solpars/.

!     The iteration process stops as soon as one of the following
!     conditions holds:
!     . relative 2-norm residual original system <= AbsTol for the
!       solution in  x,
!     . 2-norm residual original system, for the solution in x,
!       has been reduced by a factor <= RedTol, or
!     . the number of iterations >= MaxNits; in this case a warning
!       message is written to standard output.

!     Arguments:
!     ==========
!     N        	i   Number of rows/columns in original matrix  A.
!     A%G     	i   Number of rows/columns in left upper part A_11.
!     A%NSchur  i   Size of 1st Schur-complement, S, of A_11 in A.
!     a%aro     i   Location of descriptor of the
!                   matrix A, of the linear system.
!     a%mlp    	i   Location of descriptor of the
!                   Preconditioner matrix Prc, stored in MLP format.
!     a%mlp%mlp%perm   i   Preconditioner permutation vector.
!     x        	io  In:  Initial guess for the solution of  A x = b
!                   Out: The solution of the linear system  A x = b
!     r        	io  In:  Right hand side of linear system, b.
!                   Out: Residual  b - A x_out

!#enddoc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'bicgstabr'

!     Local Variables:
!     ================
!     nflops           Number of floating point operations in
!                      double precision to prevent integer overflow.
!     nrmb             = ||b||
!     nrmr             = ||b - A x^(i)|| = ||r^(i)|| in i-th iteration
!     nrmr0            = ||b - A x^(0)||
!     tol              = MAX (||b||*AbsTol, ||r^(0)||*RedTol)

INTEGER 					:: i, ier
INTEGER 					:: nnzd, nnzoff, nnzlas, nza, nzp
DOUBLE PRECISION 				:: nflops
DOUBLE PRECISION 				:: alpha, beta, rho1, rho, omega, tt, rtv
DOUBLE PRECISION 				:: nrmb, nrmr, nrmr0, tol
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: p, rt, s, psh, v

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
  
  IF (outlev >= 4) THEN
    PRINT 9000 , 'Iter.', 'Rel. residual', 'Reduction residual', 'Flops/unknown'
  END IF
  
ALLOCATE( p(1:a%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( rt(1:a%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( s(1:a%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( psh(1:a%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( v(1:a%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
  
  nrmb = norm( r )
!     {  nrmb = ||b||  }
  
  nflops = nflops + 2 * a%nschur
  
  
!     Compute  r^(0) = b - A x^(0)  for the initial guess
!     x^(0)  in  x:
  
  CALL matvecp (a, x, p)
!     {  p = A x^(0)  }
  
  r = r - p
!     {  r = r^(0) = b - A x^(0)  }
  
  nrmr0 = norm( r )
!     {  nrmr0 = ||b - A x^(0)||  }
  
  nflops = nflops + 2 * nza + 3 * a%nschur
  
  
  tol = MAX(nrmb*abstol, nrmr0*redtol)
  
  
  IF (nrmb <= neglgbl) THEN
    x = 0.0D0
!        {  x = 0  ,  r = b = b - A x  }
    nrmr = nrmb
    
    IF (outlev >= 4) THEN
      PRINT 9001 , 0, nrmr, '---', nflops/DBLE(a%n)
    END IF
    
  ELSE
    nrmr  = nrmr0
    
    IF (outlev >= 4) THEN
      PRINT 9002 , 0, nrmr/nrmb, nrmr/nrmr0, nflops/DBLE(a%n)
    END IF
  
  IF ( nrmr > tol)  THEN
  
! iteration should start
  
!     Choose  rt  such that  <rt,r^(0)> != 0:
  rt = r
!     {  rt = r^(0)  }
  
  
  alpha  = 0.0D0
  omega  = 1.0D0
  rho    = 1.0D0
  
  DO   i = 1, maxnits
    rho1 = rho
    rho  = dot ( rt, r )

!   {  rho = rho^(i) = <rt,r^(i-1)>  ,  rho1 = rho^(i-1)  }
    
    nflops = nflops + 2 * a%nschur
    
    IF (ABS(rho) <= (machprc*machprc)) THEN
!           ZERO (machine precision reached), ALGORITHM FAILS
      IF (outlev >= 1) THEN
        PRINT '(/, A, 2X, A, A, /, 3X, A, X, I3, A, 1P, 2(/, 3X, A, E12.5), /)', &
          'Warning from subroutine', rounam, '!',  &
          'Stopped at iteration', i, ', because a too small scalar value.',  &
          '2-norm residual original system:          ', nrmr,  &
          'relative 2-norm residual original system: ', nrmr/nrmb
      END IF
      EXIT
    END IF
    
    
    IF (i == 1) THEN
      p = r
!           {  p = p^(i) = r^(i-1)  }
    ELSE
!           {  i >= 2  }
      beta = (rho/rho1)*(alpha/omega)
!           {  beta = (rho^(i)/rho^(i-1)) * (alpha^(i-1)/omega^(i-1))
      
!           p := r + beta * (p - omega * v)
      p = r + beta * (p - omega * v)
!           {  p = p^(i) = r^(i-1) +
!                          beta * (p^(i-1) - omega^(i-1) * v^(i-1))  }
      
      nflops = nflops + 4 * a%nschur
    END IF
    
    
    psh = p
    
    CALL solve (a, psh)

!        {  psh = phat = inv(Prc) p^(i)  ,  y = ?  }
    
    CALL matvecp (a, psh, v)
!        {  v^(i) = A inv(Prc) p^(i)  ,  psh = phat  }
    
    nflops = nflops + 2 * nza + 2 * nzp
    
    
    rtv = dot ( rt, v )

!   {  rtv = <rt,v^(i)>  }
    
    IF (ABS(rtv) <= (machprc*machprc)) THEN
!           ZERO (machine precision reached), ALGORITHM FAILS
      IF (outlev >= 1) THEN
        PRINT '(/, A, 2X, A, A, /, 3X, A, X, I3, A, 1P, 2(/, 3X, A, E12.5), /)', &
          'Warning from subroutine', rounam, '!',  &
          'Stopped at iteration', i, ', because a too small scalar value.',  &
          '2-norm residual original system:          ', nrmr,  &
          'relative 2-norm residual original system: ', nrmr/nrmb
      END IF
      EXIT
    END IF
    
    alpha = rho / rtv
!        {  alpha^(i) = rho^(i) / <rt,v^(i)>  }
    
    s = r - alpha * v
!        {  s = s^(i) = r^(i-1) - alpha^(i) * v^(i)  }
    
    nflops = nflops + 4 * a%nschur
    
    !f
    x = x - alpha * psh
!        {  x = xhalf = x^(i-1) - alpha^(i) * phat  }
    
    nflops = nflops + 2 * a%nschur
    
    
    psh = s
!        {  psh = s^(i)  }
    
    CALL solve (a, psh)

!        {  psh = shat = inv(Prc) s^(i)  ,  y = ?  }
    
    CALL matvecp (a, psh, r)
!        {  r =  t = A inv(Prc) s^(i)  ,  psh = shat  }
    
    nflops = nflops + 2 * nza + 2 * nzp
    
    
    tt = SUM( r**2 )
!        {  tt = <t,t>  }
    
    IF (ABS(tt) <= (machprc*machprc)) THEN
!           ZERO (machine precision reached), ALGORITHM FAILS
      IF (outlev >= 1) THEN
        PRINT '(/, A, 2X, A, A, /, 3X, A, X, I3, A, 1P, (/, 3X, A, E12.5), /)', &
          'Warning from subroutine', rounam, '!',  &
          'Stopped at iteration', i, ', because a too small scalar value.',  &
          '2-norm residual original system:          ', nrmr,  &
          'relative 2-norm residual original system: ', nrmr/nrmb
      END IF
      EXIT
    END IF
    
    omega = dot ( r, s ) / tt

!   {  omega = omega^(i) = <t,s^(i)> / <t,t>  }
    
    nflops = nflops + 4 * a%nschur
    
    
    x = x + omega * psh
    r = s - omega * r
!        {  x^(i) = x^(i-1) + alpha^(i) * phat + omega^(i) * shat  ,
!           r^(i) = s^(i) - omega^(i) * t                          }
    
    nflops = nflops + 4 * a%nschur
    
    
!        Calculate 2-norm residual in nrmr:
    
    nrmr = norm( r )
    
    nflops = nflops + 2 * a%nschur
    
    
    IF (outlev >= 4) THEN
      PRINT 9002 , i, nrmr/nrmb, nrmr/nrmr0, nflops/DBLE(a%n)
    END IF
    
    IF ( nrmr <= tol ) EXIT

!   If maximum number of iterations reached:

    IF (( i == maxnits ) .AND. (outlev >= 1)) THEN
      PRINT '(/, A, 2X, A, A, /, 3X, A, 1P, 2(/, 3X, A, E12.5), /)',  &
          'Warning from subroutine', rounam, '!',  &
          'Residual too large after maximum number of iterations.',  &
          '2-norm residual original system:          ', nrmr,  &
          'relative 2-norm residual original system: ', nrmr/nrmb
    END IF

  END DO

  END IF
  END IF
  
DEALLOCATE( p, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( rt, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( s, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( psh, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( v, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

  
  RETURN
  
!     Format statements:
  9000 FORMAT (/, a6, 2(2X, a18), 2X, a14)
  9001 FORMAT (i6, 1P, (2X, 6X, e12.6), (2X, a18), 2X, 6X, 0PF8.2)
  9002 FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
  
END SUBROUTINE bicgstabr

END MODULE