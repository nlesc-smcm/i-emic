!#Begindoc
 
      MODULE m_bgskrylov
        ! contains all variants of krylov methods used by bgsprec

      CONTAINS
        
!*********************************************************************************
       SUBROUTINE bgsbicgstab (n, ixa, ixbgs, numbgs, x, r,variant,type)
          
          USE m_build
          USE m_chkcnt
          USE m_wennz
          USE m_glbpars
          USE m_solpars
          USE m_matvec

          USE m_bgsprec
          USE m_bgskit
          
          INTEGER				, INTENT(IN)            :: n
          TYPE (anymatrix)			, POINTER		:: ixa
          TYPE (bgsprec)                 	, POINTER		:: ixbgs
          INTEGER, DIMENSION(1:n) 		, INTENT(IN)            :: numbgs
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: r
          INTEGER :: variant ! 1-3
          INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
         
          !     Implementation of the Bi-CGSTAB method, using a Left
          !     preconditioner.  A variation of the algorithm described in
          !        "Templates for the Solution of Linear Systems:
          !         Building Blocks for Iterative Methods" by Richard Barrett,...,
          !     to solve the linear system
          !         A x = b
          
          !     Arguments:
          !     ==========
          !     N        i   Number of rows/columns in original matrix  A.
          !     ixA      i   Location of descriptor of the
          !                  matrix A, of the linear system.
          !     ixBGS    i   Location of descriptor of the
          !                  Preconditioner matrix Prc, stored in MLP format.
          !     numBGS   i   Preconditioner permutation vector.
          !     x        io  In:  Initial guess for the solution of  A x = b
          !                  Out: The solution of the linear system  A x = b
          !     r        io  In:  Right hand side of linear system, b.
          !                  Out: Residual  b - A x_out
          
          !     The array arguments
          !        p, rt, s, v, y
          !     are work arrays, and are used as local variables inside the
          !     subroutine.
          
          !#enddoc
          
          !     Local Parameters:
          !     =================
          
          CHARACTER (LEN=8), PARAMETER :: rounam = 'bgsbicgstab'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 				:: i, nschur
          DOUBLE PRECISION 			:: alpha, beta, rho1, rho, omega, tt, rtv
          DOUBLE PRECISION 			:: nrmr0, nrmr, tol
          DOUBLE PRECISION, DIMENSION(1:n)	:: p, rt, s, v
          
          IF (outlev >= 2) THEN
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual'
          END IF
          !     Compute residual preconditioned system,
          !     r^(0) = inv(Prc) (b - A x^(0)),  for the initial guess
          !     x^(0)  in  x:
          p = 0
          CALL matvec (n, 1D0, ixa, x, p)
          !     {  p = A x^(0)  }
          r = r - p
          !     {  r = b - A x^(0)  }
          CALL apply_bilu (ixbgs, r, numbgs,variant, type)
          !     {  r = r^(0) = inv(Prc) (b - A x^(0))  ,  y = ?  }
          !     Calculate 2-norm of residual preconditioned system:
          nrmr0 = DSQRT ( DOT_PRODUCT( r, r ))
          nrmr = nrmr0
          !     {  nrmr = nrmr0 = ||r^(0)||  }
          IF (outlev >= 4) THEN
             PRINT 9002 , 0, nrmr, nrmr/nrmr0
          END IF
          tol = MAX(abstol, nrmr0*redtol)
          IF (nrmr0 <= tol)  GO TO 200
          !     Choose  rt  such that  <rt,r^(0)> != 0:
          rt = r
          !     {  rt = r^(0)  }
          alpha  = 0.0D0
          omega  = 1.0D0
          rho    = 1.0D0
          DO   i = 1, maxnits
             rho1 = rho
             rho  = DOT_PRODUCT( rt, r )
             !        {  rho = rho^(i) = <rt,r^(i-1)>  ,  rho1 = rho^(i-1)  }
             IF (DABS(rho) <= (machprc*machprc)) THEN
                !           ZERO (machine precision reached), ALGORITHM FAILS
                GO TO 110
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
             END IF
             v = 0 
             CALL matvec (n, 1D0, ixa, p, v)
             !        {  v^(i) = A p^(i)  ,  p = phat^(i)  }
             !        Solve  Prc z = v ,  the solution is stored in  v
             CALL apply_bilu (ixbgs, v, numbgs,variant,type)
             !        {  v^(i) = inv(Prc) A p^(i)  ,  p = phat^(i)  ,  y = ?  }
             rtv = DOT_PRODUCT( rt, v )
             !        {  rtv = <rt,v^(i)>  }
             IF (DABS(rtv) <= (machprc*machprc)) THEN
                !           ZERO (machine precision reached), ALGORITHM FAILS
                GO TO 110
             END IF
             alpha = rho / rtv
             !        {  alpha^(i) = rho^(i) / <rt,v^(i)>  }
             s = r - alpha * v
             !        {  s = s^(i) = r^(i-1) - alpha^(i) * v^(i)  }
             x = x + alpha * p
             !        {  x = xhalf = x^(i-1) + alpha^(i) * phat^(i)  }
             r = 0 
             CALL matvec (n, 1D0, ixa, s, r)
             !        {  r = A s^(i)  ,  s = shat^(i)  }
             CALL apply_bilu (ixbgs, r, numbgs,variant,type)
             !        {  r = t = inv(Prc) A s^(i)  ,  s = shat^(i)  ,  y = ?  }
             tt = DOT_PRODUCT( r, r )
             !        {  tt = <t,t>  }
             IF (DABS(tt) <= (machprc*machprc)) THEN
                !           ZERO (machine precision reached), ALGORITHM FAILS
                GO TO 110
             END IF
             omega = DOT_PRODUCT( r, s ) / tt
             !        {  omega = omega^(i) = <t,s^(i)> / <t,t>  }
             x = x + omega * s
             r = s - omega * r
             !        {  x^(i) = x^(i-1) +
             !                   alpha^(i) * phat^(i) + omega^(i) * shat^(i)  ,
             !           r^(i) = s^(i) - omega^(i) * t                        }
             !        Calculate 2-norm preconditioned residual in nrmr:
             nrmr = DSQRT ( DOT_PRODUCT( r, r ))
             IF (outlev >= 2) THEN
                PRINT 9002 , i, nrmr, nrmr/nrmr0
             END IF
             IF ( nrmr <= tol )  GO TO 200
          END DO
          !     Maximum number of iterations reached:
          IF (outlev >= 1) THEN
             WRITE(6,9004)  'Warning from subroutine', rounam, '!'
             WRITE(6,*)  'Residual too large after maximum number of iterations.'
             WRITE(6,9005)  '2-norm residual preconditioned system: ', nrmr
          END IF
          GO TO 200
          !     Method fails:
110       CONTINUE
          IF (outlev >= 1) THEN
             WRITE(6,9004) 'Warning from subroutine', rounam, '!'
             WRITE(6,*) 'Stopped at iteration', i, ', because a too small scalar value.'
             WRITE(6,9005)  '2-norm residual preconditioned system: ', nrmr
          END IF
200       CONTINUE
          WRITE(99,9003) i,rounam,mgmres,nrmr
          WRITE(6,9003) i,rounam,mgmres,nrmr
          IF (outlev >= 1) THEN
             WRITE(6,9006) INT(REAL(sppiter)/REAL(i+1))
          END IF
          !     Normal Return:
          RETURN
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i3,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
9006      FORMAT (' Average number of iterations on depth-averaged spp:',i4)

1000      CONTINUE
          
          !     End of  bgsbicgstab
        END SUBROUTINE bgsbicgstab
        
!*********************************************************************************

        SUBROUTINE bgsgmres (m, n, ixa, ixbgs, numbgs, x, b,variant, type)
          
          USE m_build
          USE m_chkcnt
          USE m_wennz
          USE m_glbpars
          USE m_solpars
          USE m_matvecp
          USE m_matvec
          USE m_solve

          USE m_bgsprec
          USE m_bgskit
          
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (anymatrix)       		, POINTER		:: ixa
          TYPE (bgsprec)			, POINTER		:: ixbgs
          INTEGER, DIMENSION(1:n) 		, INTENT(IN)            :: numbgs
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: b
          INTEGER :: variant ! 1-3
          INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
          
          !     Implementation of the GMRES(M) method, as described in
          !        "Templates for the Solution of Linear Systems:
          !         Building Blocks for Iterative Methods" by Richard Barrett,...,
          !     to solve the linear system
          !         A x = b
          
          !     The NxN-matrix  A  is the Schur-complement of the
          !     submatrix  A_11  of the NxN-matrix  A;  A_11 = A(1:G,1:G).
          !     The Block-Gauss-Seidel preconditioner, bgsprc, is stored and can be reached
          !     through the descriptor at 'ixbgs'.
          
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
          !     M        i   Dimension of the subspace of R^(N),  M <= N.
          !     N        i   Number of rows/columns in original matrix  A.
          !     ixA      i   Location of descriptor of the
          !                  matrix A, of the linear system.
          !     ixBGS    i   Location of descriptor of the
          !                  Preconditioner matrix BGS, stored in bgsprec format.
          !     numBGS   i   Preconditioner permutation vector.
          !     x        io  In:  Initial guess for the solution of  A x = b
          !                  Out: The solution of the linear system  A x = b
          !     b        i   Right-hand side of the matrix equation  A x = b
          
          !#enddoc
          
          !     Global Parameters:
          !     ==================
          
          !     Local Parameters:
          !     =================
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'bgsgmres'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 					:: i, j, k, iter
          INTEGER 					:: nnzd, nnzoff, nnzlas, na, nza, nzp
          INTEGER, DIMENSION(1:n) 	                :: numprc2

          DOUBLE PRECISION 				:: nflops
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w
          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v
          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y
          
#ifdef debug
          !     TRACE INFORMATION
          PRINT '(A, X, A)' , 'Entry:', rounam
#endif
          
          
          !     Check the preconditioner partition storage and count number of
          !     non-zeros:
          nnzd = 1
          nnzoff = 1
          nnzlas = 1
          nflops = 0
          nzp    = nnzd + nnzoff + nnzlas
          
          !     Calculate the number of nonzeros in [ixA]:
          na = 1; nza = 1;
          
          IF (outlev >= 2) THEN
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual', 'Flops/unknown'
          END IF
          
          
          sppiter = 0 ! set counter for subiterations on saddle point problem to zero
          iter   = 0
          factor = 1.0D0
          DO i = 1,n
             numprc2(i) = i
          END DO
          
          
          !     Starting point new iteration of GMRES(M)
100       CONTINUE
          
          !        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
          !        for the initial guess   x^(iter)  in  x:
          
          w = 0
          CALL matvec (n, 1D0,ixa, x, w)
          !        {  w = A x^(iter)  }
          
          !        Calculate residual
          w = b - w
          !        {  w = b - A x^(iter)  }
          
          CALL apply_bilu (ixbgs, w, numbgs,variant,type)

          !        {  w = r^(iter) = inv(Prc) (b - A x^(iter))  }
          
          nflops = nflops + 2 * nza + n + 2 * nzp
          
          
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          
          nflops = nflops + 2 * n
          

          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             
             IF (outlev >= 4) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0, nflops/DBLE(n)
             END IF
             
             tol = MAX(abstol, nrmr0*redtol)
          END IF
          
          IF (nrmr <= tol)  GO TO 200
          
          
          !        Initialise the Basis for the Generalised Minimal RESidual
          !        Method:

          s(1) = nrmr
          s(2:m+1) = 0.0D0
          !        {  s = ||r^(0)|| e_1  }

          ff = 1.0D0 / nrmr
          v(:,1) = w * ff
          
          !        {  V(:,1) = w / ||w||  }
          
          nflops = nflops + n

          DO  i = 1, MIN(m,maxnits-iter)
             
             !           Solve  w  from  Prc w = A V(:,i):
             
             w=0
             CALL matvec (n,1D0, ixa, v(:,i), w)
             !           {  w = A V(:,i)  }
             
             CALL apply_bilu (ixbgs, w, numbgs,variant,type)
             !           {  w = r^(i) = inv(Prc) A v^(i))  }
             
             nflops = nflops + 2 * nza + 2 * nzp
             
             
             !           Start the modified Gram-Schmidt process:
             DO k = 1, i
                h(k,i) = DOT_PRODUCT( w, v(:,k) )
                !              {  H(k,i) = <w,V(:,k)>  }
                w = w - h(k,i) * v(:,k)
                !              {  w := w - H(k,i) V(:,k)  }
             END DO
             
             nflops = nflops + i * 4 * n
             
             
             h(i+1,i) = DSQRT ( DOT_PRODUCT( w, w ))
             
             ff = 1.0D0 / h(i+1,i)
             v(:,i+1) = w * ff
             !           {  H(i+1,i) = ||w||  ,  V(:,i+1) = w / H(i+1,i)  }
             
             nrmr = h(i+1,i)
             
             nflops = nflops + 3 * n
             
             
             !           Apply the Givens rotations J_1, J_2, ..., J_(i-1)
             !           on  H(1:i,i):
             DO k = 1, i-1
                hki      = h(k,i)
                h(k  ,i) = jc(k)*hki - js(k)*h(k+1,i)
                h(k+1,i) = js(k)*hki + jc(k)*h(k+1,i)
             END DO
             
             !           Construct Givens rotation  J_i  acting on H(i:i+1,i), such
             !           that (i+1)st component of J_i H(:,i) is 0:
             ff = DSQRT (h(i,i)**2 + h(i+1,i)**2)
             jc(i) =  h(i  ,i) / ff
             js(i) = -h(i+1,i) / ff
             
             !           Apply  J_i  on  H(i:i+1,i):
             h(i  ,i) = jc(i) * h(i,i) - js(i) * h(i+1,i)
             h(i+1,i) = 0
             
             !           Apply  J_i  on  s(i:i+1):
             ff     = s(i)
             s(i  ) = jc(i) * ff - js(i) * s(i+1)
             s(i+1) = js(i) * ff + jc(i) * s(i+1)
             
             nflops = nflops + (i-1) * 6 + 15
             
             iter   = iter + 1
             factor = factor * DABS(js(i))
             nrmr   = factor * nrmr0
             
             IF (nrmr <= tol  .Or. i == m  .Or. iter == maxnits) THEN
                !              Compute the solution in  x:
                
                !              Compute  y  as the solution of  H y = s(1:i), in which
                !              the upper  i x i  triangular part of  H  has  H(j,k)
                !              as its elements.
                DO j = i, 1, -1
                   ff = s(j) - DOT_PRODUCT( h(j,j+1:i), y(j+1:i) )
                   y(j) = ff / h(j,j)
                   
                   !                 Update the solution in  x:
                   x = x + v(:,j) *  y(j)
                   !                 {  x = x^(0) + V(:,j:i) y(j:i)  }
                END DO
                
                nflops = nflops + i*(i-1) + i * (1 + 2 * n)
             END IF
             
             IF (outlev >= 2) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0, nflops/DBLE(n)
             END IF
             
             IF (nrmr <= tol)  GO TO 200
          END DO
          
          IF (iter < maxnits)  GO TO 100
          
          
          !     Maximum number of iterations reached:
          IF (outlev >= 1) THEN
             WRITE(6,9004)  'Warning from subroutine', rounam, '!'
             WRITE(6,*)  'Residual too large after maximum number of iterations.'
             WRITE(6,9005)  '2-norm residual preconditioned system: ', nrmr
          END IF
          
200       CONTINUE
         
          WRITE(99,9003) iter,rounam,mgmres,nrmr
          WRITE(6,9003) iter,rounam,mgmres,nrmr

          IF (outlev >= 1) THEN
             WRITE(6,9006) INT(REAL(sppiter)/REAL(iter+1))
          END IF
          !     Normal Return:
          RETURN
          
          
          !     Format statements:
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i2,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
9006      FORMAT (' Average number of iterations on depth-averaged spp:',i4)

1000      CONTINUE
          
          !     End of  gmres
        END SUBROUTINE bgsgmres

!*********************************************************************************
        SUBROUTINE bgsfgmres (m, n, ixa, ixbgs, numprc, x, b,variant, type)
          
          USE m_build
          USE m_glbpars
          USE m_solpars
          USE m_matvec
          USE m_bgsprec
          USE m_bgskit
          
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (anymatrix)			, POINTER		:: ixa
          TYPE (bgsprec)			, POINTER		:: ixbgs
          INTEGER, DIMENSION(1:n) 		, INTENT(IN)            :: numprc
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: b
          INTEGER :: variant ! 1-3
          INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
          
          !     Implementation of the GMRES(M) method, as described in
          !        "Templates for the Solution of Linear Systems:
          !         Building Blocks for Iterative Methods" by Richard Barrett,...,
          !     to solve the linear system
          !         A x = b
          
          !     Arguments:
          !     ==========
          !     M        i   Dimension of the subspace of R^(N),  M <= N.
          !     N        i   Number of rows/columns in original matrix  A.
          !     ixA      i   Location of descriptor of the
          !                  matrix A, of the linear system.
          !     ixBGS    i   Location of descriptor of the
          !                  Preconditioner matrix BGS, stored in bgsprec format.
          !     numPrc   i   Preconditioner permutation vector.
          !     x        io  In:  Initial guess for the solution of  A x = b
          !                  Out: The solution of the linear system  A x = b
          !     b        i   Right-hand side of the matrix equation  A x = b
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'bgsfgmres'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 					:: i, j, k, iter
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w, xsum
!          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v, z
!          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: v, z
          DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y
          
          ALLOCATE(v(n,m+1), z(n,m+1), h(m+1,m+1))
          IF (outlev >= 1) THEN
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual'
          END IF
          sppiter = 0 ! set counter for subiterations on saddle point problem to zero
          iter   = 0
          factor = 1.0D0
          !     Starting point new iteration of GMRES(M)
100       CONTINUE
          !        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
          !        for the initial guess   x^(iter)  in  x:
          w = 0
          CALL matvec (n, 1D0, ixa, x, w)
          !        {  w = A x^(iter)  }
          !        Calculate residual
          w = b - w
          !        {  w = b - A x^(iter)  }
          ! Right-preconditioning so no prec application
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||(b - A x^(0))||  }
             IF (outlev >= 1) THEN
                IF (nrmr>0) THEN
                   PRINT 9002 , iter, nrmr, nrmr/nrmr0
                ELSE
                   PRINT 9002 , iter, nrmr, nrmr
                END IF
             END IF
             tol = MAX(abstol, nrmr0*redtol)
          END IF
          ! do at least one iteration
          ! IF (nrmr <= tol)  GO TO 200
          !        Initialise the Basis for the Generalised Minimal RESidual
          !        Method:
          s(1) = nrmr
          s(2:m+1) = 0.0D0
          !        {  s = ||r^(0)|| e_1  }
          ff = 1.0D0 / nrmr
          v(:,1) = w * ff
          !        {  V(:,1) = w / ||w||  }
          DO  i = 1, MIN(m,maxnits-iter)
             z(:,i) = v(:,i)
             CALL apply_bilu (ixbgs, z(:,i), numprc, variant,type)
             !           {  z^(i) = inv(Prc) v^(i)  }
             w = 0
             CALL matvec (n, 1D0, ixa, z(:,i), w)
             !           {  w = A z^(i)  }
             !           Start the modified Gram-Schmidt process:
             DO k = 1, i
                h(k,i) = DOT_PRODUCT( w, v(:,k) )
                !              {  H(k,i) = <w,V(:,k)>  }
                w = w - h(k,i) * v(:,k)
                !              {  w := w - H(k,i) V(:,k)  }
             END DO
             h(i+1,i) = DSQRT ( DOT_PRODUCT( w, w ))
             ff = 1.0D0 / h(i+1,i)
             v(:,i+1) = w * ff
             !           {  H(i+1,i) = ||w||  ,  V(:,i+1) = w / H(i+1,i)  }
             nrmr = h(i+1,i)
             !           Apply the Givens rotations J_1, J_2, ..., J_(i-1)
             !           on  H(1:i,i):
             DO k = 1, i-1
                hki      = h(k,i)
                h(k  ,i) = jc(k)*hki - js(k)*h(k+1,i)
                h(k+1,i) = js(k)*hki + jc(k)*h(k+1,i)
             END DO
             !           Construct Givens rotation  J_i  acting on H(i:i+1,i), such
             !           that (i+1)st component of J_i H(:,i) is 0:
             ff = DSQRT (h(i,i)**2 + h(i+1,i)**2)
             jc(i) =  h(i  ,i) / ff
             js(i) = -h(i+1,i) / ff
             !           Apply  J_i  on  H(i:i+1,i):
             h(i  ,i) = jc(i) * h(i,i) - js(i) * h(i+1,i)
             h(i+1,i) = 0
             !           Apply  J_i  on  s(i:i+1):
             ff     = s(i)
             s(i  ) = jc(i) * ff - js(i) * s(i+1)
             s(i+1) = js(i) * ff + jc(i) * s(i+1)
             iter   = iter + 1
             factor = factor * DABS(js(i))
             nrmr   = factor * nrmr0
             IF (nrmr <= tol  .Or. i == m  .Or. iter == maxnits) THEN
                !              Compute the solution in  x:
                !              Compute  y  as the solution of  H y = s(1:i), in which
                !              the upper  i x i  triangular part of  H  has  H(j,k)
                !              as its elements.
                DO j = i, 1, -1
                   ff = s(j) - DOT_PRODUCT( h(j,j+1:i), y(j+1:i) )
                   y(j) = ff / h(j,j)
                   !                 Update the solution in  x:
                   x = x + z(:,j) *  y(j)
                   !                 {  x = x^(0) + Z(:,j:i) y(j:i)  }
                END DO
                w = 0
                CALL matvec (n, 1D0, ixa, x, w)
                !        {  w = A x^(iter)  }
                !        Calculate residual
                w = b - w
             END IF
             IF (outlev >= 1) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
             END IF
             IF (nrmr <= tol)  GO TO 200
          END DO
          IF (iter < maxnits)  GO TO 100
          !     Maximum number of iterations reached:
          IF (outlev >= 1) THEN
             WRITE(6,9004)  'Warning from subroutine', rounam, '!'
             WRITE(6,*)  'Residual too large after maximum number of iterations.'
             WRITE(6,9005)  '2-norm residual preconditioned system: ', nrmr
          END IF
          DEALLOCATE(v, z, h)
          
200       CONTINUE
          WRITE(99,9003) iter,rounam,mgmres,nrmr
          WRITE(6,9003) iter,rounam,mgmres,nrmr
          IF (outlev >= 1) THEN
             WRITE(6,9006) INT(REAL(sppiter)/REAL(iter+1))
          END IF
          !     Normal Return:
          GOTO 1000
!          RETURN
          !     Format statements:
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i3,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
9006      FORMAT (' Average number of iterations on depth-averaged spp:',i4)

1000      CONTINUE
          !     End of  bgsfgmres
        END SUBROUTINE bgsfgmres

!****************************************************************************

      END MODULE m_bgskrylov
      
