!#Begindoc
 
      MODULE m_bgsgmres
        ! contains all variants of krylov methods used by bgsprec

      CONTAINS
        
      SUBROUTINE sppgmres (m, n, ixbgs, x, b, prectype)
          ! modules from mrilu libraries
          USE m_build
          USE m_chkcnt
          USE m_wennz
          USE m_glbpars
          USE m_solpars
          ! modules of bgs library
          USE m_bgsprec
          USE m_sparsekit
 
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (bgsprec)			, POINTER		:: ixbgs
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: b
          CHARACTER (LEN = 2)                   , INTENT(IN)            :: prectype          
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
          !     x        io  In:  Initial guess for the solution of  A x = b
          !                  Out: The solution of the linear system  A x = b
          !     b        i   Right-hand side of the matrix equation  A x = b
          !     prectype i   type of saddle point preconditioner used
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'sppgmres'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 					:: i, j, k, iter
          INTEGER, DIMENSION(1:n) 	                :: numprc2
          DOUBLE PRECISION 				:: nflops
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w
          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v
          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y
          
          IF (outlev >= 3) THEN
             WRITE(6,*) 'Solving SPP iteratively'
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual'
          END IF
          
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
          CALL matvec_spp (ixbgs, x, w)
          !        {  w = A x^(iter)  }
          !        Calculate residual
          w = b - w
          !        {  w = b - A x^(iter)  }
          CALL apply_sppprec (ixbgs, w, prectype)
          !        {  w = r^(iter) = inv(Prc) (b - A x^(iter))  }
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             IF (outlev >= 3) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
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
          DO  i = 1, MIN(m,maxnits-iter)
             !           Solve  w  from  Prc w = A V(:,i):
             w = 0
             CALL matvec_spp (ixbgs, v(:,i), w)
             CALL apply_sppprec (ixbgs, w, prectype)
             !           {  w = r^(i) = inv(Prc) A v^(i))  }
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
                   x = x + v(:,j) *  y(j)
                   !                 {  x = x^(0) + V(:,j:i) y(j:i)  }
                END DO
             END IF
             IF (outlev >= 3) THEN
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
200       CONTINUE
          !     Normal Return:
          sppiter = sppiter+iter
          IF (outlev>=2) THEN
             WRITE (6,9003) iter, rounam, mgmres, nrmr
          END IF
          RETURN
          
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i2,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
       
1000      CONTINUE
          !     End of  sppgmres
        END SUBROUTINE sppgmres

!*********************************************************************************
        
      SUBROUTINE sppfgmres (m, n, ixbgs, x, b, prectype)
          ! modules from mrilu libraries
          USE m_build
          USE m_chkcnt
          USE m_wennz
          USE m_glbpars
          USE m_solpars
          ! modules of bgs library
          USE m_bgsprec
          USE m_sparsekit
 
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (bgsprec)			, POINTER		:: ixbgs
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: b
          CHARACTER (LEN = 2)                   , INTENT(IN)            :: prectype          

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
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'sppfgmres'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 					:: i, j, k, iter
          DOUBLE PRECISION 				:: nflops
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w, xsum
          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v, z
          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y
          
          IF (outlev >= 3) THEN
             WRITE(6,*) ' Solving system spp iteratively' 
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual'
          END IF
          iter   = 0
          factor = 1.0D0
          !     Starting point new iteration of GMRES(M)
100       CONTINUE
          !        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
          !        for the initial guess   x^(iter)  in  x:
          w = 0
          CALL matvec_spp (ixbgs, x, w)
          !        {  w = A x^(iter)  }
          !        Calculate residual
          w = b - w
          !        {  w = b - A x^(iter)  }
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             IF (outlev >= 3) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
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
          DO  i = 1, MIN(m,maxnits-iter)
             z(:,i) = v(:,i)
             CALL apply_sppprec (ixbgs, z(:,i), prectype)
             !           {  z^(i) = inv(Prc) v^(i)  }
             w = 0 
             CALL matvec_spp (ixbgs, z(:,i), w)
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
                CALL matvec_spp (ixbgs, z(:,i), w)
                !        {  w = A x^(iter)  }
                !        Calculate residual
                w = b-w
             END IF
             IF (outlev >= 3) THEN
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
200       CONTINUE
          !     Normal Return:
          IF (outlev>=2) THEN
             WRITE(6,9003) iter, rounam, mgmres, nrmr
          END IF
          RETURN
          RETURN
          
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i2,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
         
1000      CONTINUE
          !     End of  sppfgmres
        END SUBROUTINE sppfgmres

!*********************************************************************************
            
        SUBROUTINE matvec_spp(ixBGS,x,b)
          ! Computes b=Ax, where A is the depth-averaged saddle point problem
          USE m_bgsprec
          USE m_csrvec
          IMPLICIT none
          ! IMPORT/EXPORT
          TYPE (bgsprec), POINTER         :: ixBGS
          REAL, DIMENSION(:), INTENT(IN)  :: x
          REAL, DIMENSION(:), INTENT(OUT) :: b
          CHARACTER (LEN = 2)             :: prectype
          ! LOCAL
          INTEGER :: n, m
          !
          b = 0
          n = ixBGS%MAuv%n
          m = ixBGS%MDuv%n
          CALL csrvec(1D0, ixBGS%MAuv, x(1:n),     b(1:n))
          CALL csrvec(1D0, ixBGS%MGuv, x(n+1:n+m), b(1:n))
          CALL csrvec(1D0, ixBGS%MDuv, x(1:n),     b(n+1:n+m))
        END SUBROUTINE matvec_spp

!*********************************************************************************
        
        SUBROUTINE apply_sppprec(ixBGS,x,prectype)
          ! computes y = inv(P)x, where P is the prectype preconditioner
          ! for the depth averaged saddle point problem
          USE m_bgsprec
          USE m_applprc
          USE m_csrvec
          IMPLICIT none
          ! IMPORT/EXPORT
          TYPE (bgsprec), POINTER           :: ixBGS
          REAL, DIMENSION(:), INTENT(INOUT) :: x
          CHARACTER (LEN = 2)               :: prectype
          ! LOCAL
          INTEGER :: n, m,ier
          REAL, DIMENSION(:), ALLOCATABLE   :: y,z,w
          !
          n = ixBGS%MAuv%n
          m = ixBGS%MDuv%n
          SELECT CASE (prectype)
          CASE ('SI')
             ALLOCATE(y(n),z(m))
             y = x(1:n); x(1:n) = 0
             CALL applprc(ixBGS%PMAuv,x(1:n),y)
             CALL csrvec(-1D0,ixBGS%MDuv,x(1:n),x(n+1:n+m))
	     z = - x(n+1:n+m); x(n+1:n+m) = 0
             CALL applprc(ixBGS%PMCp,x(n+1:n+m),z)
             y = 0
             CALL csrvec(1D0,ixBGS%MGuv,x(n+1:n+m),y)
             CALL csrvec(-1D0,ixBGS%PDMAuv,y,x(1:n))
             DEALLOCATE(y,z)
          CASE ('SR')
             ALLOCATE(y(n),z(n+m),w(m))
             z = x
             y = 0
             CALL csrvec(1D0,ixBGS%PDMAuv,z(1:n),y)
             CALL csrvec(-1D0,ixBGS%MDuv,y,z(n+1:n+m))
             w = - z(n+1:n+m); z(n+1:n+m) = 0
             CALL applprc(ixBGS%PMCp,z(n+1:n+m),w)
             CALL csrvec(-1D0,ixBGS%MGuv,z(n+1:n+m),z(1:n))
             y = z(1:n); z(1:n) = 0
             CALL applprc(ixBGS%PMAuv,z(1:n),y)
             ! compute correction on x
             CALL csrvec(-1D0,ixBGS%MAuv,z(1:n),x(1:n))
             CALL csrvec(-1D0,ixBGS%MGuv,z(n+1:n+m),x(1:n))
             CALL csrvec(-1D0,ixBGS%MDuv,z(1:n),x(n+1:n+m))
             ! apply SIMPLE-iteration to corrected x
             y = x(1:n); x(1:n) = 0
             CALL applprc(ixBGS%PMAuv,x(1:n),y)
             CALL csrvec(-1D0,ixBGS%MDuv,x(1:n),x(n+1:n+m))
	     w = - x(n+1:n+m); x(n+1:n+m) = 0
             CALL applprc(ixBGS%PMCp,x(n+1:n+m),w)
             y = 0
             CALL csrvec(1D0,ixBGS%MGuv,x(n+1:n+m),y)
             CALL csrvec(-1D0,ixBGS%PDMAuv,y,x(1:n))
             ! end of SIMPLE-iteration
             x = x + z
             DEALLOCATE(y,z,w)
          CASE ('WS')
             x(n+1:n+m) = x(n+1:n+m)*ixBGS%omega
             ALLOCATE(y(n))
             y = x(1:n); x(1:n) = 0
             CALL applprc(ixBGS%PMAuv,x(1:n),y)
             DEALLOCATE(y)
          CASE ('ES')
             x(n+1:n+m) = x(n+1:n+m)*ixBGS%omega
             CALL csrvec(-1D0,ixBGS%MGuv,x(n+1:n+m),x(1:n))
             ALLOCATE(y(n))
             y = x(1:n); x(1:n) = 0
             CALL applprc(ixBGS%PMAuv,x(1:n),y)
             DEALLOCATE(y)
          CASE ('GD')
             x(n+1:n+m) = x(n+1:n+m)*ixBGS%omega
             ALLOCATE(y(n))
             y = x(1:n); x(1:n) = 0
             CALL applprc(ixBGS%PMAuv,x(1:n),y)
             DEALLOCATE(y)
          CASE ('AC')
             CALL csrvec(ixBGS%omega, ixBGS%MGuv ,x(n+1:n+m),x(1:n))
             x(n+1:n+m) = -ixBGS%omega*x(n+1:n+m)
             ALLOCATE(y(n))
             y = x(1:n); x(1:n) = 0
             CALL applprc(ixBGS%PMAuv,x(1:n),y)
             DEALLOCATE(y)
             CALL csrvec(ixBGS%omega,ixBGS%MDuv,x(1:n),x(n+1:n+m))
          CASE DEFAULT
             STOP 'No spp preconditioner selected'
          END SELECT
        END SUBROUTINE apply_sppprec
        
!*********************************************************************************
        
        SUBROUTINE prcgmres (m, n, ixa, ixprc, numprc, x, b)
          
          USE m_build
          USE m_chkcnt
          USE m_wennz
          USE m_glbpars
          USE m_solpars
          USE m_matvecp
          USE m_matvec
          USE m_csrvec
          USE m_solve
          USE m_applprc
          
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (csrmatrix)			, POINTER		:: ixa
          TYPE (prcmatrix)		        , POINTER		:: ixprc
          INTEGER, DIMENSION(1:n) 		, INTENT(IN)            :: numprc
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: b
          
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
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'prcgmres'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 					:: i, j, k, iter
          DOUBLE PRECISION 				:: nflops
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w,z
          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v
          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y
          
          IF (outlev >= 3) THEN
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual'
          END IF
          iter   = 0
          factor = 1.0D0
          !     Starting point new iteration of GMRES(M)
100       CONTINUE
          !        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
          !        for the initial guess   x^(iter)  in  x:
          w = 0
          CALL csrvec (1D0, ixa, x, w)
          !        {  w = A x^(iter)  }
          !        Calculate residual
          z = b - w
          !        {  w = b - A x^(iter)  }
          CALL applprc (ixprc, w, z)
          !        {  w = r^(iter) = inv(Prc) (b - A x^(iter))  }
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             IF (outlev >= 3) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
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
          DO  i = 1, MIN(m,maxnits-iter)
             !           Solve  w  from  Prc w = A V(:,i):
             z = 0 
             CALL csrvec ( 1D0, ixa, v(:,i), z)
             !           {  w = A V(:,i)  }
             CALL applprc (ixprc, w, z)
             !           {  w = r^(i) = inv(Prc) A v^(i))  }
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
             END IF
             IF (outlev >= 3) THEN
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
200       CONTINUE
          !     Normal Return:
          IF (outlev>=2) THEN
             WRITE(6,9003) iter, rounam, mgmres, nrmr
          END IF
          RETURN
          RETURN
          
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i2,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)

1000      CONTINUE
          !     End of  prgmres
        END SUBROUTINE prcgmres

!*********************************************************************************
        
        SUBROUTINE prcfgmres (m, n, ixa, ixprc, numprc, x, b)
          
          USE m_build
          USE m_chkcnt
          USE m_wennz
          USE m_glbpars
          USE m_solpars
          USE m_matvecp
          USE m_matvec
          USE m_csrvec
          USE m_solve
          USE m_applprc
          
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (csrmatrix)			, POINTER		:: ixa
          TYPE (prcmatrix)		        , POINTER		:: ixprc
          INTEGER, DIMENSION(1:n) 		, INTENT(IN)            :: numprc
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: b
          
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
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'prcfgmres'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 					:: i, j, k, iter
          DOUBLE PRECISION 				:: nflops
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w, xsum
          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v, z
          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y

          IF (outlev >= 3) THEN
             WRITE(6,*) ' Solving system ATS iteratively' 
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual', 'Flops/unknown'
          END IF
          iter   = 0
          factor = 1.0D0
          !     Starting point new iteration of GMRES(M)
100       CONTINUE
          !        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
          !        for the initial guess   x^(iter)  in  x:
          w = 0
          CALL csrvec (1D0, ixa, x, w)
          !        {  w = A x^(iter)  }
          !        Calculate residual
          w = b - w
          !        {  w = b - A x^(iter)  }
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             IF (outlev >= 3) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
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
          DO  i = 1, MIN(m,maxnits-iter)
             CALL applprc (ixprc, z(:,i), v(:,i))
             !           {  z^(i) = inv(Prc) v^(i)  }
             w = 0 
             CALL csrvec(1D0, ixA, z(:,i),w)
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
                CALL csrvec(1D0, ixA, z(:,i), w)
                !        {  w = A x^(iter)  }
                !        Calculate residual
                w = b-w
             END IF
             IF (outlev >= 3) THEN
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
200       CONTINUE
          !     Normal Return:
          IF (outlev>=2) THEN
             WRITE(6,9003) iter, rounam, mgmres, nrmr
          END IF
          RETURN
          RETURN
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i2,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
         
1000      CONTINUE
          !     End of  prcfgmres
        END SUBROUTINE prcfgmres

!*********************************************************************************
        
        SUBROUTINE SCgmres (m, n, ixbgs, x, b, var, type)
          
          USE m_build
          USE m_glbpars
          USE m_solpars
          USE m_bgsprec
          
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (bgsprec)		        , POINTER		:: ixbgs
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: b
          INTEGER :: var  ! 1 = TS, 2 = uv, 3 = w/p - schurcomplement
          INTEGER :: type ! 2 = depth averaged spp, 3 = modified simpler

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
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'SCgmres'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 					:: i, j, k, iter
          DOUBLE PRECISION 				:: nflops
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w,z
          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v
          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y
          
          IF (outlev >= 3) THEN
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual'
          END IF
          iter   = 0
          factor = 1.0D0
          !     Starting point new iteration of GMRES(M)
100       CONTINUE
          !        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
          !        for the initial guess   x^(iter)  in  x:
          w = 0
          CALL matvecSC(ixbgs, x, w, var, type)
          !        {  w = A x^(iter)  }
          !        Calculate residual
          z = b - w
          !        {  w = b - A x^(iter)  }
          CALL precSC(ixbgs, z, w, var, type)
          !        {  w = r^(iter) = inv(Prc) (b - A x^(iter))  }
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             IF (outlev >= 3) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
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
          DO  i = 1, MIN(m,maxnits-iter)
             !           Solve  w  from  Prc w = A V(:,i):
             z = 0 
             CALL matvecSC(ixbgs, v(:,i), z, var, type)
             !           {  w = A V(:,i)  }
             CALL precSC(ixbgs, z, w, var, type)
             !           {  w = r^(i) = inv(Prc) A v^(i))  }
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
             END IF
             IF (outlev >= 3) THEN
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
200       CONTINUE
          !     Normal Return:
          IF (outlev>=2) THEN
             WRITE(6,9003) iter, rounam, mgmres, nrmr
          END IF
          RETURN
          RETURN
          
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i2,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)

1000      CONTINUE
          !     End of  SCgmres
        END SUBROUTINE SCgmres

!*********************************************************************************
        SUBROUTINE SCfgmres (m, n, ixbgs, x, b, var, type)
          
          USE m_build
          USE m_chkcnt
          USE m_wennz
          USE m_glbpars
          USE m_solpars
          USE m_matvecp
          USE m_matvec
          USE m_csrvec
          USE m_solve
          USE m_applprc
          USE m_bgsprec
          USE m_sparsekit
          
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (bgsprec)		        , POINTER		:: ixbgs
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: b
          INTEGER :: var  ! 1 = TS, 2 = uv, 3 = w/p - schurcomplement
          INTEGER :: type ! 2 = depth averaged spp, 3 = modified simpler
          
          !     Implementation of the GMRES(M) method, as described in
          !        "Templates for the Solution of Linear Systems:
          !         Building Blocks for Iterative Methods" by Richard Barrett,...,
          !     to solve the linear system
          !         A x = b
          !     used for the Schur complement

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
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'SCfgmres'
          
          !     Local Variables:
          !     ================
          !     NSchur           Size of 1st Schur-complement, S, of A_11 in A.
          !     nflops           Number of floating point operations in
          !                      double precision to prevent integer overflow.
          !     nrmr             = ||r^(i)|| = ||inv(Prc) (b - A x^(i))||
          !     nrmr0            = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||
          !     tol              = MAX (AbsTol, ||r^(0)||*RedTol)
          
          INTEGER 					:: i, j, k, iter
          DOUBLE PRECISION 				:: nflops
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w, xsum
          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v, z
          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y
          
          IF (outlev >= 3) THEN
             WRITE(6,*) ' Solving system SC iteratively, variant', var
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual'
          END IF
          iter   = 0
          factor = 1.0D0
          !     Starting point new iteration of GMRES(M)
100       CONTINUE
          !        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
          !        for the initial guess   x^(iter)  in  x:
          w = 0
          CALL matvecSC(ixbgs, x, w, var, type)
          !        {  w = A x^(iter)  }
          !        Calculate residual
          w = b - w
          !        {  w = b - A x^(iter)  }
          !        {  w = r^(iter) = inv(Prc) (b - A x^(iter))  }
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             IF (outlev >= 3) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
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
          DO  i = 1, MIN(m,maxnits-iter)
             CALL precSC(ixbgs, v(:,i), z(:,i),var,type)
             !           {  z^(i) = inv(Prc) v^(i)  }
             w = 0 
             CALL matvecSC(ixbgs, z(:,i), w,var,type)
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
                CALL matvecSC(ixbgs, z(:,i), w,var,type)
                !        {  w = A x^(iter)  }
                !        Calculate residual
                w = b-w
             END IF
             IF (outlev >= 3) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
             END IF
             IF (nrmr <= tol)  GO TO 200
          END DO
          IF (iter < maxnits)  GO TO 100
          !     Maximum number of iterations reached:
          IF (outlev >= 1) THEN
!             WRITE(6,9004)  'Warning from subroutine', rounam, '!'
!             WRITE(6,*)  'Residual too large after maximum number of iterations.'
!             WRITE(6,9005)  '2-norm residual preconditioned system: ', nrmr
          END IF
200       CONTINUE
          !     Normal Return:
          IF (outlev>=2) THEN
             WRITE(6,9003) iter, rounam, mgmres, nrmr
          END IF
          RETURN
          RETURN
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i2,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
         
1000      CONTINUE
          !     End of  SCfgmres
        END SUBROUTINE SCfgmres

!*********************************************************************************
        SUBROUTINE matvecSC(ixBGS,x,w,var,type)
          USE m_bgsprec
          ! computes w = (ATS+...)*x
          ! the action of the ATS schur complement of the bgs matrix 
          TYPE (bgsprec), POINTER :: ixBGS
          DOUBLE PRECISION, DIMENSION(:), INTENT(IN OUT) :: x,w
          INTEGER :: var  ! 1 = TS, 2 = uv, 3 = w/p - schurcomplement
          INTEGER :: type ! 2 = depth averaged spp, 3 = modified simpler
          ! 
          IF (var.eq.1) THEN
             CALL matvecSCTS(ixBGS,x,w,type)
          ELSEIF (var.eq.2) THEN
             CALL matvecSCspp(ixBGS,x,w,type)
          ELSEIF ((var.eq.3).and.(.not.ixBGS%transposed)) THEN
             CALL matvecSCp(ixBGS,x,w,type)
          ELSEIF ((var.eq.3).and.(ixBGS%transposed)) THEN
             CALL matvecSCw(ixBGS,x,w,type)
          ELSE
             WRITE(99,*) " error: var bigger than 3 in matvecSC " 
             STOP
          END IF
        END SUBROUTINE matvecSC

!*********************************************************************************
        SUBROUTINE matvecSCTS(ixBGS,x,w,type)
          USE m_bgsprec
          USE m_csrvec
          USE m_applprc
          USE m_sparsekit
          USE m_bgspars
          USE m_inisol
          ! computes w = (ATS+...)*x
          ! the action of the ATS-schur complement of the bgs matrix 
          TYPE (bgsprec), POINTER :: ixBGS
          DOUBLE PRECISION, DIMENSION(:), INTENT(IN OUT) :: x,w
          INTEGER :: type ! 2 = depth averaged spp, 3 = modified simpler
          ! 
          INTEGER :: m,n
          DOUBLE PRECISION, DIMENSION(:), POINTER :: y,z,v
          !
          IF (((ixBGS%BwTS%nnz==0).AND.(ixBGS%BuvTS%nnz==0)).OR.((ixBGS%BTSw%nnz==0).AND.(ixBGS%BTSuv%nnz==0))) THEN
             w = 0
             CALL csrvec(1D0,ixBGS%ATS,x,w)
          ELSE
             nTS = ixBGS%ATS%n; nw = ixBGS%BwTS%n;
             nuv = ixBGS%Auv%n; nzp = ixBGS%M2Duv%n;
             ALLOCATE(y(nw),z(nw)); y=0; z=0; w=0;
             CALL csrvec(1D0,ixBGS%BwTS,x,y)
             CALL usolve(ixBGS%Ap,y,z)        ! z = Ap\BwTS*xTS
             DEALLOCATE(y)
             ALLOCATE(y(nuv)); y=0;
             CALL csrvec(1D0,ixBGS%BuvTS,x,y) ! y = [BuvTS*xTS]
             CALL csrvec(-1D0,ixBGS%Guv1,z,y) ! y = [BuvTS*xTS-Guv1*z]
             DEALLOCATE(z)
             IF (type.eq.1) THEN
                ALLOCATE(z(nuv)); z=0
                CALL applprc(ixBGS%PAuv,z,y)  ! z = PAuv\y
                DEALLOCATE(y)
             ELSEIF (type.eq.2) THEN
                ALLOCATE(z(nuv+nzp)); z=0;
                ! CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
                ! CALL bigsppfgmres(mgmres,nuv+nzp,ixBGS,y,z) ! solve bigSPP: y = SPP\z
                z(1:nuv) = y; z(nuv+1:nuv+nzp) = 0;
                CALL apply_bigsppprec (ixbgs, z, 'SR') ! z = K\[y;0]
                ! CALL inisol(cgtype,mgmres,maxnits,loctolabs,loctolred)
                DEALLOCATE(y)
             END IF

             ALLOCATE(y(nw),v(nw)); y = 0; v = 0
             CALL csrvec(1D0,ixBGS%Duv1,z(1:nuv),v)
             CALL lsolve(ixBGS%Aw,v,y)
             w = 0
             CALL csrvec(1D0,ixBGS%ATS,x,w)
             CALL csrvec(-1D0,ixBGS%BTSuv,z(1:nuv),w)
             CALL csrvec(+1D0,ixBGS%BTSw,y,w)
             ! y = ATS*x-BTSuv*z+BTSw*(Aw\(Duv1*z));
             DEALLOCATE(y,z,v)
          END IF

        END SUBROUTINE matvecSCTS

!*********************************************************************************
        SUBROUTINE matvecSCspp(ixBGS,x,w,type)
          USE m_bgsprec
          USE m_csrvec
          USE m_applprc
          USE m_sparsekit
          USE m_bgspars
          USE m_inisol
          IMPLICIT NONE
          ! computes w = (SPP+...)*x
          ! the action of the SPP-schur complement of the bgs matrix 
          TYPE (bgsprec), POINTER :: ixBGS
          DOUBLE PRECISION, DIMENSION(:), INTENT(IN OUT), TARGET :: x,w
          INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
          ! 
          INTEGER :: m,n,nTS,nuv,nw,nzp,nz
          DOUBLE PRECISION, DIMENSION(:), POINTER :: y,z,v, xuv,xzp
          !
          w = 0
          CALL matvec_bigspp(ixBGS,x,w)

          IF (((ixBGS%BwTS%nnz==0).AND.(ixBGS%BuvTS%nnz==0)).OR.((ixBGS%BTSw%nnz==0).AND.(ixBGS%BTSuv%nnz==0))) THEN
             CONTINUE
          ELSE
             ! NB need to be implemented!
             nTS = ixBGS%ATS%n; nw = ixBGS%BwTS%n;
             nuv = ixBGS%Auv%n; nzp = ixBGS%M2Duv%n; nz = nuv+nzp;
             xuv => x(1:nuv); xzp => x(nuv+1:nzp);
             ALLOCATE(y(nw),z(nw)); y=0; z=0;
             CALL csrvec(1D0,ixBGS%Duv1,xuv,z)
             CALL lsolve(ixBGS%Aw,y,z);            ! y = Aw\(Duv1*xuv)
             DEALLOCATE(z); ALLOCATE(z(nTS)); z=0;
             CALL csrvec(1D0,ixBGS%BTSuv,xuv,z) 
             CALL csrvec(-1D0,ixBGS%BTSw,y,z)      ! z = BTSuv*xuv - BTSw*(Aw\Duv1*xuv)
             DEALLOCATE(y); ALLOCATE(y(nTS)); y=0;
             CALL applprc(ixBGS%PATS,y,z)          ! y = PATS\z
             DEALLOCATE(z); ALLOCATE(z(nw));  z=0;
             CALL csrvec(-1D0,ixBGS%BuvTS,y,w(1:nuv)) ! add first part to w
             !
             CALL csrvec(1D0,ixBGS%BwTS,y,z);      ! z = BwTS*y
             DEALLOCATE(y); ALLOCATE(y(nw));  y=0;
             CALL usolve(ixBGS%Ap,z,y)             ! y = Ap\z
             CALL csrvec(1D0,ixBGS%Guv1,y,w(1:nuv)) ! add second part to w
             DEALLOCATE(y,z)
          END IF
        END SUBROUTINE matvecSCspp

!*********************************************************************************
        SUBROUTINE matvecSCp(ixBGS,x,w,type)
          USE m_bgsprec
          USE m_csrvec
          USE m_applprc
          USE m_sparsekit
          USE m_bgspars
          USE m_inisol
          IMPLICIT NONE
          ! computes w = (Ap+...)*x
          ! the action of the Ap schur complement of the bgs matrix 
          TYPE (bgsprec), POINTER :: ixBGS
          DOUBLE PRECISION, DIMENSION(:), INTENT(IN OUT) :: x,w
          INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
          ! 
          INTEGER :: m,n, nTS, nw,nuv,nzp, np
          DOUBLE PRECISION, DIMENSION(:), POINTER :: y,z,v,u
          !
          IF (((ixBGS%BwTS%nnz==0).AND.(ixBGS%BuvTS%nnz==0)).OR.((ixBGS%BTSw%nnz==0).AND.(ixBGS%BTSuv%nnz==0))) THEN
             w = 0
             CALL csrvec(1D0,ixBGS%Ap,x,w)
          ELSE
             nTS = ixBGS%ATS%n; nw = ixBGS%BwTS%n;
             nuv = ixBGS%Auv%n; nzp = ixBGS%M2Duv%n; np = ixBGS%Ap%n;
             ALLOCATE(y(nuv),z(nuv),v(nuv+nzp)); y=0; z=0; v=0; w=0;
             CALL csrvec(1D0,ixBGS%Guv1,x,y)  ! y = Guv*x
             IF (type.eq.1) THEN
                CALL applprc(ixBGS%PAuv,z,y)  ! z = PAuv\y
             ELSEIF (type.eq.2) THEN
                v = 0; v(1:nuv) = y;
                !         CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
                !         CALL inisol(cgtype,mgmres,maxnits,1D-6,1D-6)
                !         CALL bigsppfgmres(mgmres,nuv+nzp,ixBGS,yuvzp,buvzp) ! solve bigSPP
                CALL apply_bigsppprec (ixbgs, v, 'SR')
                z = v(1:nuv)                  ! v = SPP\[y;0]
             END IF
             DEALLOCATE(y,v)
             ALLOCATE(y(nTS),v(nw),u(nw)); y=0; v=0; u=0;
             CALL csrvec(1D0,ixBGS%Duv1,z,v)  ! v = Duv1*z 
             CALL lsolve(ixBGS%Aw,v,u)        ! u = Aw\v = Aw\Duv1*z 
             CALL csrvec(-1D0,ixBGS%BTSw,u,y) ! y = - BTSw*(Aw\Duv1*z)
             CALL csrvec(1D0,ixBGS%BTSuv,z,y) ! y = BTSuv*z - BTSw*(Aw\Duv1*z)
             DEALLOCATE(z,v,u)
             AllOCATE(z(nTS),v(nw)); z=0; v=0;
             CALL applprc(ixBGS%PATS,z,y)     ! z = ATS\y
             CALL csrvec(1D0,ixBGS%BwTS,z,w)  ! w = BwTS*z
             CALL csrvec(1D0,ixBGS%Ap,x,w)    ! w = Ap*x + BwTS*(ATS\y)
             ! w = Ap*x + BwTS*(ATS\ (BTSuv-BTSw*(Aw\Duv1)))*(SPP\Guv1*x))
             DEALLOCATE(z,v)
          END IF
        END SUBROUTINE matvecSCp
        
!*********************************************************************************
        SUBROUTINE matvecSCw(ixBGS,x,w,type)
          USE m_bgsprec
          USE m_csrvec
          USE m_applprc
          USE m_sparsekit
          USE m_bgspars
          USE m_inisol
          IMPLICIT NONE
          ! computes w = (Ap+...)*x
          ! the action of the Ap schur complement of the bgs matrix 
          TYPE (bgsprec), POINTER :: ixBGS
          DOUBLE PRECISION, DIMENSION(:), INTENT(IN OUT) :: x,w
          INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
          ! 
          INTEGER :: m,n, nTS, nw,nuv,nzp, np
          DOUBLE PRECISION, DIMENSION(:), POINTER :: y,z,v,u
          !
          IF (((ixBGS%BwTS%nnz==0).AND.(ixBGS%BuvTS%nnz==0)).OR.((ixBGS%BTSw%nnz==0).AND.(ixBGS%BTSuv%nnz==0))) THEN
             w = 0
             CALL csrvec(1D0,ixBGS%Ap,x,w)
          ELSE
             ! NB need to be implemented!
             nTS = ixBGS%ATS%n; nw = ixBGS%BwTS%n;
             nuv = ixBGS%Auv%n; nzp = ixBGS%M2Duv%n; np = ixBGS%Ap%n;
             ALLOCATE(y(nuv),z(nuv),v(nuv+nzp)); y=0; z=0; v=0; w=0;
             CALL csrvec(1D0,ixBGS%Guv1,x,y)  ! y = Guv*x
             IF (type.eq.1) THEN
                CALL applprc(ixBGS%PAuv,z,y)  ! z = PAuv\y
             ELSEIF (type.eq.2) THEN
                v = 0; v(1:nuv) = y;
                !         CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
                !         CALL inisol(cgtype,mgmres,maxnits,1D-6,1D-6)
                !         CALL bigsppfgmres(mgmres,nuv+nzp,ixBGS,yuvzp,buvzp) ! solve bigSPP
                CALL apply_bigsppprec (ixbgs, v, 'SR')
                z = v(1:nuv)                  ! v = SPP\[y;0]
             END IF
             DEALLOCATE(y,v)
             ALLOCATE(y(nTS),v(nw),u(nw)); y=0; v=0; u=0;
             CALL csrvec(1D0,ixBGS%Duv1,z,v)  ! v = Duv1*z 
             CALL lsolve(ixBGS%Aw,v,u)        ! u = Aw\v = Aw\Duv1*z 
             CALL csrvec(-1D0,ixBGS%BTSw,u,y) ! y = - BTSw*(Aw\Duv1*z)
             CALL csrvec(1D0,ixBGS%BTSuv,z,y) ! y = BTSuv*z - BTSw*(Aw\Duv1*z)
             DEALLOCATE(z,v,u)
             AllOCATE(z(nTS),v(nw)); z=0; v=0;
             CALL applprc(ixBGS%PATS,z,y)     ! z = ATS\y
             CALL csrvec(1D0,ixBGS%BwTS,z,w)  ! w = BwTS*z
             CALL csrvec(1D0,ixBGS%Ap,x,w)    ! w = Ap*x + BwTS*(ATS\y)
             ! w = Ap*x + BwTS*(ATS\ (BTSuv-BTSw*(Aw\Duv1)))*(SPP\Guv1*x))
             DEALLOCATE(z,v)
          END IF
        END SUBROUTINE matvecSCw
        
!*********************************************************************************
        SUBROUTINE precSC(ixBGS, x, w, var,type)
   !          CALL precSC(ixbgs, v(:,i), z(:,i),var)
          USE m_bgsprec
          USE m_csrvec
          USE m_applprc
          USE m_sparsekit
          USE m_bgspars
          USE m_inisol
          IMPLICIT NONE
          ! computes w = Prec\x
          TYPE (bgsprec), POINTER :: ixBGS
          DOUBLE PRECISION, DIMENSION(:), INTENT(IN OUT) :: x,w
          INTEGER :: var  ! 1 = TS, 2 = uv, 3 = w/p - schurcomplement
          INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
          DOUBLE PRECISION, DIMENSION(:), POINTER :: z
          ! 
          IF (var.eq.1) THEN
             CALL applprc (ixbgs%PATS, w , x)
          ELSEIF (var.eq.2) THEN             
             w = x
             CALL apply_bigsppprec(ixBGS,w,'SR')
          ELSEIF ((var.eq.3).and.(.not.ixBGS%transposed)) THEN
             CALL usolve(ixBGS%Ap,x,w)
          ELSEIF ((var.eq.3).and.(ixBGS%transposed)) THEN
             CALL lsolve(ixBGS%Aw,x,w)
          ELSE
             WRITE(99,*) " error: var bigger than 3 in precSC " 
             STOP
          END IF
        END SUBROUTINE precSC

!*********************************************************************************

      SUBROUTINE bigsppfgmres (m, n, ixbgs, x, b)
          ! modules from mrilu libraries
          USE m_build
          USE m_chkcnt
          USE m_wennz
          USE m_glbpars
          USE m_solpars
          ! modules of bgs library
          USE m_bgsprec
          USE m_sparsekit
 
          INTEGER				, INTENT(IN)		:: m
          INTEGER				, INTENT(IN)            :: n
          TYPE (bgsprec)			, POINTER		:: ixbgs
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: x
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: b

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
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'bigsppfgmres'
          
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
          DOUBLE PRECISION, DIMENSION(1:n,1:m+1)	:: v, z
          DOUBLE PRECISION, DIMENSION(1:m+1,1:m+1)	:: h
          DOUBLE PRECISION, DIMENSION(1:m)		:: jc, js
          DOUBLE PRECISION, DIMENSION(1:m+1)		:: s,y
          
          IF (outlev >= 3) THEN
             WRITE(6,*) ' Solving system bigSPP iteratively' 
             PRINT 9000 , 'Iter.', 'Residual', 'Reduction residual', 'Flops/unknown'
          END IF
          iter   = 0
          factor = 1.0D0
          !     Starting point new iteration of GMRES(M)
100       CONTINUE
          !        Compute residual  w = r^(iter) = inv(Prc) (b - A x^(iter))
          !        for the initial guess   x^(iter)  in  x:
          w = 0
          CALL matvec_bigspp (ixbgs, x, w)
          !        {  w = A x^(iter)  }
          !        Calculate residual
          w = b - w
          !        {  w = b - A x^(iter)  }
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             IF (outlev >= 3) THEN
                PRINT 9002 , iter, nrmr, nrmr/nrmr0
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
          DO  i = 1, MIN(m,maxnits-iter)
             z(:,i) = v(:,i)
             CALL apply_bigsppprec (ixbgs, z(:,i), 'SR')
             !           {  z^(i) = inv(Prc) v^(i)  }
             w = 0 
             CALL matvec_bigspp (ixbgs, z(:,i), w)
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
                CALL matvec_bigspp (ixbgs, z(:,i), w)
                !        {  w = A x^(iter)  }
                !        Calculate residual
                w = b-w
             END IF
             IF (outlev >= 3) THEN
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
200       CONTINUE
          !     Normal Return:
          IF (outlev>=2) THEN
             WRITE(6,9003) iter, rounam, mgmres, nrmr
          END IF
          RETURN
          RETURN
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i2,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
         
1000      CONTINUE
          !     End of  bigsppfgmres
        END SUBROUTINE bigsppfgmres

!*********************************************************************************
            
        SUBROUTINE matvec_bigspp(ixBGS,x,b)
          ! Computes b=Ax, where A is the big saddle point problem
          ! with depth-averaged pressure
          USE m_bgsprec
          USE m_csrvec
          IMPLICIT none
          ! IMPORT/EXPORT
          TYPE (bgsprec), POINTER         :: ixBGS
          REAL, DIMENSION(:), INTENT(IN)  :: x
          REAL, DIMENSION(:), INTENT(OUT) :: b
          ! LOCAL
          INTEGER :: n, m
          !
          b = 0
          n = ixBGS%Auv%n
          m = ixBGS%M2Duv%n
          CALL csrvec(1D0, ixBGS%Auv,     x(1:n),     b(1:n))
          CALL csrvec(1D0, ixBGS%GuvM1T, x(n+1:n+m), b(1:n))
          CALL csrvec(1D0, ixBGS%M2Duv,  x(1:n),     b(n+1:n+m))
        END SUBROUTINE matvec_bigspp

!*********************************************************************************
        
        SUBROUTINE apply_bigsppprec(ixBGS,x,prectype)
          ! computes y = inv(P)x, where P is the prectype preconditioner
          ! for the depth averaged saddle point problem
          USE m_bgsprec
          USE m_applprc
          USE m_csrvec
          IMPLICIT none
          ! IMPORT/EXPORT
          TYPE (bgsprec), POINTER           :: ixBGS
          REAL, DIMENSION(:), INTENT(INOUT) :: x
          CHARACTER (LEN = 2)               :: prectype
          ! LOCAL
          INTEGER :: n, m,ier
          REAL, DIMENSION(:), ALLOCATABLE   :: y,z,w
          !
          n = ixBGS%Auv%n
          m = ixBGS%M2Duv%n
          SELECT CASE (prectype)
          CASE ('SI')
             ALLOCATE(y(n),z(m))
             y = x(1:n); x(1:n) = 0
             CALL applprc(ixBGS%PAuv,x(1:n),y)
             CALL csrvec(-1D0,ixBGS%M2Duv,x(1:n),x(n+1:n+m))
	     z = - x(n+1:n+m); x(n+1:n+m) = 0
             CALL applprc(ixBGS%PCp,x(n+1:n+m),z)
             y = 0
             CALL csrvec(1D0,ixBGS%GuvM1T,x(n+1:n+m),y)
             CALL csrvec(-1D0,ixBGS%PDAuv,y,x(1:n))
             DEALLOCATE(y,z)
          CASE ('SR')
             ALLOCATE(y(n),z(n+m),w(m))
             z = x
             y = 0
             CALL csrvec(1D0,ixBGS%PDAuv,z(1:n),y)
             CALL csrvec(-1D0,ixBGS%M2Duv,y,z(n+1:n+m))
             w = - z(n+1:n+m); z(n+1:n+m) = 0
             CALL applprc(ixBGS%PCp,z(n+1:n+m),w)
             CALL csrvec(-1D0,ixBGS%GuvM1T,z(n+1:n+m),z(1:n))
             y = z(1:n); z(1:n) = 0
             CALL applprc(ixBGS%PAuv,z(1:n),y)
             ! compute correction on x
             CALL csrvec(-1D0,ixBGS%Auv,z(1:n),x(1:n))
             CALL csrvec(-1D0,ixBGS%GuvM1T,z(n+1:n+m),x(1:n))
             CALL csrvec(-1D0,ixBGS%M2Duv,z(1:n),x(n+1:n+m))
             ! apply SIMPLE-iteration to corrected x
             y = x(1:n); x(1:n) = 0
             CALL applprc(ixBGS%PAuv,x(1:n),y)
             CALL csrvec(-1D0,ixBGS%M2Duv,x(1:n),x(n+1:n+m))
	     w = - x(n+1:n+m); x(n+1:n+m) = 0
             CALL applprc(ixBGS%PCp,x(n+1:n+m),w)
             y = 0
             CALL csrvec(1D0,ixBGS%GuvM1T,x(n+1:n+m),y)
             CALL csrvec(-1D0,ixBGS%PDAuv,y,x(1:n))
             ! end of SIMPLE-iteration
             x = x + z
             DEALLOCATE(y,z,w)
          CASE DEFAULT
             STOP 'No spp preconditioner selected'
          END SELECT
        END SUBROUTINE apply_bigsppprec
        
!*********************************************************************************

      END MODULE m_bgsgmres

