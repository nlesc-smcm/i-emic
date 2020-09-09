      MODULE m_mrilusvp
!
      CONTAINS

!**************************************************************
      SUBROUTINE compute_singular_vectors(ixA,svp1,svp2)
!     Computes the singular vectors in case of MRILU preconditioner
      USE m_build
      USE m_sparsekit
      USE m_wfree
      USE m_bgskit
      USE m_bgsprec
      IMPLICIT none
!      INCLUDE 'par.com'
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER          :: ixA
      REAL, DIMENSION(:)                 :: svp1,svp2
!     LOCAL
      TYPE (bgsprec), POINTER            :: ixBGS
      INTEGER, DIMENSION(:), ALLOCATABLE :: ord,iord
      INTEGER                            :: dumw, dump, ier
!     construct subblocks and depth-averaged system
      WRITE (6,*) "Computing block ordering for matrix"
      ! 1) construct block-ordering and reoderder matrix
      ! 1a) construct general block ordering
      WRITE (6,*) "Constructing block ordering"
      ALLOCATE(ord(ixA%n))
      CALL block_ordering2(ord)
      ! 1b) detect dummy w and p points
      WRITE (6,*) "Detecting dummy w and p points"
      CALL detect_dummies(ixA,ord,dumw,dump)
      ! 1c) construct inverse ordering
      WRITE (6,*) "Constructing inverse ordering"
      ALLOCATE(iord(ixA%n))
      CALL invert_ordering(ord,iord)
      ! 1d) reorder matrix
      WRITE (6,*) "Reordering matrix"
      CALL reorder_matrix(ixA,ord,iord)
      ! 2) extract gradient and divergence blocks
      WRITE (6,*) "Splitting matrix"
      CALL extract_grad_div(ixA,ixBGS,dumw,dump)
      ! 2c) compute depth-averaging matrices
      WRITE (6,*) "Computing depth-averaged equations"
      CALL depth_averaging(ixBGS,3)
      ! 2f) construct singular vectors of the pressure
      WRITE (6,*) "Constructing singular vectors for the pressure"
      IF (ixBGS%transposed) THEN
         CALL build_svp3(ixBGS)
      ELSE
         CALL build_svp(ixBGS)
      END IF
!      CALL frs_svps(ixBGS,1)
      ! 3) insert
      CALL compute_big_svps(ixBGS,svp1,svp2,dump)      
      ! 4) restore original ordering of the matrix
      WRITE (6,*) "Reorder matrix back"
      CALL reorder_matrix(ixA,iord,ord)
      svp1 = svp1(iord)
      svp2 = svp2(iord)
      CALL delete_BGSprec_svp(ixBGS)
      DEALLOCATE(ord,iord)

      END SUBROUTINE compute_singular_vectors

!**************************************************************
      SUBROUTINE block_ordering2(ord)
!     Builds the blockordering (uv,w,p,T,S) and its inverse
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      INTEGER, DIMENSION(:) :: ord
!     LOCAL
      integer NN,i
!
      NN = ndim/6
      DO i = 0, NN-1
         ord(2*i + UU)  = nun*i + UU
         ord(2*i + VV)  = nun*i + VV
         ord(i+1+(WW-1)*NN) = nun*i + WW
         ord(i+1+(PP-1)*NN) = nun*i + PP
         ord(2*i+(TT-1)*NN+1) = nun*i + TT
         ord(2*i+(TT-1)*NN+2) = nun*i + SS
      ENDDO
      END SUBROUTINE block_ordering2

!**************************************************************
      SUBROUTINE extract_grad_div(ixA,ixBGS,dumw,dump)
!     Stores the following blocks of the matrix in 
!     the block-Gauss-Seidel preconditioner
!     
!     [ Auv    0       0     0      Guv     BuvTS ]
!     [ 0      Adumw   0     0      Bdumwp  0     ]
!     [ 0      0       0     0      Gw      BwTS  ]
!     [ 0      0       0     Adump  0       0     ]
!     [ Duv    Bpdumw  Dw    0      0       0     ]
!     [ BTSuv  0       BTSw  0      0       ATS   ]
!     
      USE m_build
      USE m_bgsprec
      USE m_sparsekit
      USE m_wfree
      IMPLICIT none
      INCLUDE 'par.com'
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER    :: ixA,ixB
      TYPE (bgsprec), POINTER      :: ixBGS
      INTEGER                      :: dumw,dump
!     LOCAL
      INTEGER :: nr, nnzA,nA,ier

      ! allocate preconditioner
      ALLOCATE(ixBGS, STAT=ier)
      nr = n*m*(l+la)
      nA = ixA%n
      ! extract the convection-diffusion-Coriolis equation
      CALL extract_submatrix(ixA,ixBGS%Guv ,1          ,2*nr,3*nr+dump+1,4*nr)
      ! extract the hydrostatic pressure equation
      CALL extract_submatrix(ixA,ixBGS%Gw   ,2*nr+dumw+1,3*nr,3*nr+dump+1,4*nr)
      ! extract the continuity equation
      CALL extract_submatrix(ixA,ixBGS%Duv  ,3*nr+dump+1,4*nr,1     ,2*nr)
      CALL extract_submatrix(ixA,ixBGS%Dw   ,3*nr+dump+1,4*nr,2*nr+dumw+1,3*nr)
      ! multiply continuity equation with -1
      ixBGS%Duv%co = - ixBGS%Duv%co
      ixBGS%Dw %co = - ixBGS%Dw %co

      ! extract dummy blocks to check tranposed or not
      CALL extract_submatrix(ixA,ixBGS%Bpdumw,3*nr+dump+1,4*nr,2*nr+1,2*nr+dumw)
      CALL extract_submatrix(ixA,ixBGS%Bdumwp,2*nr+1,2*nr+dumw,3*nr+dump+1,4*nr)
      
      IF (ixBGS%Bpdumw%nnz.GT.0) THEN
         ixBGS%transposed = .false.
      ELSEIF (ixBGS%Bdumwp%nnz.GT.0) THEN
         ixBGS%transposed = .true.
         WRITE(*,*) "Solving transposed system"
      ELSE
         ixBGS%transposed = .false.
      END IF

      END SUBROUTINE extract_grad_div

!**************************************************************
      SUBROUTINE delete_BGSprec_svp(ixBGS)
!     frees the memory claimed to store the BGS preconditioner
      USE m_bgsprec
      USE m_wfree
      USE m_bgspars
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER :: ixBGS
!
      !  subblocks:
      CALL csrfree(ixBGS%Guv)
      CALL csrfree(ixBGS%Gw )
      CALL csrfree(ixBGS%Duv)
      CALL csrfree(ixBGS%Dw )
      !  dummy matrices
      CALL csrfree(ixBGS%Bpdumw)
      CALL csrfree(ixBGS%Bdumwp)
      !  depth-averaging
      CALL csrfree(ixBGS%M1)
      CALL csrfree(ixBGS%M2)
      CALL csrfree(ixBGS%M1T)
      CALL csrfree(ixBGS%GuvM1T)
      CALL csrfree(ixBGS%M2Duv)
!
      DEALLOCATE(ixBGS%svp1)
      DEALLOCATE(ixBGS%svp2)
      DEALLOCATE(ixBGS%svp3)
      DEALLOCATE(ixBGS%svp4)
!
      DEALLOCATE(ixBGS)
      NULLIFY(ixBGS)
!
      END SUBROUTINE delete_BGSprec_svp

!**************************************************************
        
        SUBROUTINE gmres_svp (m, n, ixa, ixprc, numprc, x, b,svp1,svp2)
          
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
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: svp1
          DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN)            :: svp2
          
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
          !     numPrc   i   Preconditioner permutation vector.
          !     x        io  In:  Initial guess for the solution of  A x = b
          !                  Out: The solution of the linear system  A x = b
          !     b        i   Right-hand side of the matrix equation  A x = b
          
          !#enddoc
          
          !     Global Parameters:
          !     ==================
          
          !     Local Parameters:
          !     =================
          
          CHARACTER (LEN=*), PARAMETER :: rounam = 'gmres_svp'
          
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

          DOUBLE PRECISION 				:: nflops
          DOUBLE PRECISION 				:: ff, factor, hki, nrmr, nrmr0, tol
          
          DOUBLE PRECISION, DIMENSION(1:n) 		:: w,z
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

          nnzd = 0; nnzoff = 0; nnzlas = 0;

          nflops = 0
          nzp    = nnzd + nnzoff + nnzlas
          
          !     Calculate the number of nonzeros in [ixA]:
          na = 1; nza = 1;

          
          
          IF (outlev >= 3) THEN
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
          z = b - w
          !        {  w = b - A x^(iter)  }
          
          CALL applprc (ixprc, w, z)
          !        {  w = r^(iter) = inv(Prc) (b - A x^(iter))  }
          w = w - dot_product(w,svp1)*svp1
          w = w - dot_product(w,svp2)*svp2

          nflops = nflops + 2 * nza + n + 2 * nzp
          
          
          !        Calculate 2-norm of preconditioned residual:
          nrmr  = DSQRT ( DOT_PRODUCT( w, w ))
          !        {  nrmr = ||w|| = ||r^(iter)||  }
          
          nflops = nflops + 2 * n
          
          
          IF (iter == 0) THEN
             nrmr0 = nrmr
             !           {  nrmr0 = ||r^(0)|| = ||inv(Prc) (b - A x^(0))||  }
             
             IF (outlev >= 3) THEN
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
             
             z = 0 
             CALL csrvec ( 1D0, ixa, v(:,i), z)
             !           {  w = A V(:,i)  }
             
             CALL applprc (ixprc, w, z)
             !           {  w = r^(i) = inv(Prc) A v^(i))  }
             w = w - dot_product(w,svp1)*svp1
             w = w - dot_product(w,svp2)*svp2
             
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
             
             IF (outlev >= 3) THEN
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
          
          
          !     Normal Return:
          IF (outlev>=2) THEN
             WRITE(6,9003) iter, rounam, mgmres, nrmr
             WRITE(99,9003) iter, rounam, mgmres, nrmr
          END IF
          RETURN
          RETURN
          
          
          !     Format statements:
9000      FORMAT (/, a6, 2(2X, a18), 2X, a14)
9002      FORMAT (i6, 1P, 2(2X, 6X, e12.6), 2X, 6X, 0PF8.2)
9003      FORMAT (' After ',i3,' iterations of ',a10,'(',i3,') 2-norm residual is',e12.5)
9004      FORMAT (a, a10,a)
9005      FORMAT (a,e12.5)
         

1000      CONTINUE
          
          !     End of  gmres

        END SUBROUTINE gmres_svp

!*********************************************************************************

      END MODULE m_mrilusvp
      
