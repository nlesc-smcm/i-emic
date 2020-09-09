      MODULE m_gsolve
!
      CONTAINS
!**************************************************************
      SUBROUTINE g_solve(rl,det,id)
!     This routines solves the equation Ax=rl for non-symmetric
!     matrix A using the block-Gauss-Seidel preconditioner in 
!     (f)gmres
      ! location of the matrix
      USE m_mat
      ! modules from F90mrilu libraries
      USE m_build
      USE m_wacsr
      USE m_wfree
      USE m_inisol
      USE m_iniglb
      USE m_iniprc
      USE m_cmpprc
      USE m_solprc
      USE m_matvec
      USE m_wrtmt
      USE m_wrtvec
      USE m_csrvec
      ! module that defines location of the matrix and its preconditioner
      USE m_matloc
!      USE m_bgsloc
      ! modules that contain block-Gauss-Seidel preconditioner routines
      USE m_bgsprec
      USE m_bgskit
      USE m_bgspars
      ! modules that contain (flexible) iterative Krylov subspace methods
      ! for the block-Gauss-Seidel preconditioner
      USE m_bgskrylov
      USE m_sparsekit
      ! module containing mrilu parameters
      USE m_mrilupars
      USE m_mrilusvp
      ! module for scaling
      USE m_scaling
      USE m_start
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      REAL, DIMENSION(1:ndim), INTENT(INOUT) :: rl
      REAL, INTENT(OUT)                      :: det
      INTEGER, INTENT(INOUT)                 :: id
!     LOCAL
      REAL                    :: timeg1,timeg2
      INTEGER                 :: ier, nnzA, errinf, NNN, restarts
      REAL XS(ndim), XEX(ndim), RL2(ndim)
      REAL     :: redtol, abstol
      TYPE (csrmatrix), POINTER          :: ixB
      
      CALL cpu_time(timeg1)
      ! set tolerances
      abstol = tolabs
      redtol = tolred
      restarts = 1       ! number of restarts in (f)gmres
      call set_prec_parameters(prec,spptype,epsw_mrilu,epsw_uv,epsw_muv,epsw_p,epsw_TS,maxiter,restarts)
      IF (id == 0) THEN
         IF (prec == 1) THEN      ! mrilu
            IF (associated(ixPrc)) THEN
               CALL csrfree(ixA)
               CALL prcfree(ixPrc)
               DEALLOCATE(svp1,svp2)
            END IF
         ELSEIF (prec == 2) THEN  ! bilu
            IF (associated(ixA)) THEN
               CALL csrfree(ixA)
               DEALLOCATE(ord,svp1,svp2)
               CALL delete_BGSprec(ixBGS,spptype)
            END IF
         END IF
         ! allocate space for ixA and put matrix A in this structure
         nnzA = begA(ndim+1)-1
         ! scale equations
         IF (scale) THEN
            IF (.not.allocated(cs)) THEN
               ALLOCATE(cs(nun),rs(nun),csa(nun),rsa(nun))
            END IF
            CALL scaling_land(begA,jcoA,coA,rs,cs)
            rsa = (/1.,1.,1.,1.,.3,1./)
            csa = (/1.,1.,1.,1.,.3,1./)
            CALL rcscale_a(rl,rs,cs,rsa,csa)
         ELSE
            ! compute and show averaged block only
            IF (.not.allocated(cs)) THEN
               ALLOCATE(cs(nun),rs(nun),csa(nun),rsa(nun))
            END IF
            CALL scaling_land(begA,jcoA,coA,rs,cs)
            WRITE(*,*) "This scaling is not used!"
         END IF
         CALL wacsr(ndim,nnzA,ixA)
         ! put matrix in csr-structure of m_build
         ixA%beg = begA
         ixA%jco = jcoA(1:nnzA)
         ixA%co  = coA(1:nnzA)
         ixA%nnz = ixA%beg(ndim+1)-1
         id = 0
         IF (prec == 1) THEN ! mrilu
            ALLOCATE(svp1(ndim),svp2(ndim))
            CALL compute_singular_vectors(ixA,svp1,svp2)
            CALL iniglb(4)
            CALL set_mrilu_pars('d')
            CALL copy_csr(ixA,ixB)
            CALL cmpprc(6,ixB,ixPrc)
            CALL cpu_time(timeg2)
            WRITE(99,*) 'TIME TO CONSTRUCT MRILU-PRECONDITIONER: ',timeg2 - timeg1
         ELSEIF (prec == 2) THEN ! bilu
            ALLOCATE(ord(ndim),svp1(ndim),svp2(ndim))
            CALL build_BGSprec(ixA,ixBGS,ord,svp1,svp2,spptype)
            CALL cpu_time(timeg2)
            WRITE(99,*) 'TIME TO CONSTRUCT BGS-PRECONDITIONER:   ',timeg2 - timeg1
         ELSE
            STOP 'no preconditioner defined'
         END IF
         CALL cpu_time(timeg1)
         XS = 0
      ELSEIF (id == 1) THEN
         IF (scale) CALL scalesol_a(rl,rs,rsa)
         ! second system solved
         XS = rb2old
         IF (scale) CALL scalesol_a(XS,1/cs,1/csa)
      ELSE
         IF (scale) CALL scalesol_a(rl,rs,rsa)
         ! first system solved with old matrix and preconditioner
         XS = 0
      END IF

      ! solve the system with (f)gmres
      IF (prec == 1) THEN ! mrilu
         CALL iniglb(outlev)
         CALL inisol(cgtype1,mgmres1,maxnits1,redtol1,abstol1)
!         CALL solprc (ixPrc, XS, RL)
         CALL gmres_svp(mgmres1,ndim,ixA,ixPrc,ixPrc%mlp%perm,XS,RL,svp1,svp2)
      ELSEIF (prec == 2) THEN ! bilu
         CALL iniglb(outlev)
         CALL inisol(cgtype,mgmres,maxnits,redtol,abstol)
         IF (cgtype==1) THEN
            CALL bgsbicgstab(ndim,csrtoany(ixA),ixBGS,ord,XS,RL,variant,spptype)
         ELSEIF (cgtype==2) THEN
            CALL bgsgmres(mgmres,ndim,csrtoany(ixA),ixBGS,ord,XS,RL,variant,spptype)
         ELSEIF (cgtype==3) THEN
            CALL bgsfgmres(mgmres,ndim,csrtoany(ixA),ixBGS,ord,XS,RL,variant,spptype)
         ELSE
            STOP 'error in gsolve: cgtype not defined'
         END IF
      END IF
      ! Set flag to use old preconditioning and old matrix
      id = id + 1
      det = 0.0D0

      IF (scale) CALL scalesol_a(XS,cs,csa)
      RL = XS
      call cpu_time(timeg2)
      write(99,*) 'TIME TO SOLVE SYSTEM WITH GMRES:        ',timeg2 - timeg1
!      write(99,*) 'after precon, id = ', id
      END SUBROUTINE g_solve

!**************************************************************

      END MODULE m_gsolve
