      SUBROUTINE eigen
      USE m_bgspars
      USE m_gsolve
      implicit none
      include 'bag.com'
      save eivec
*     CONSTANT
      integer  jmax,jmin,kmax,mm,ll
*     parameter(jmax=20,jmin=10,kmax=nf,mm=5,ll=2)
      parameter(jmax=30,jmin=2*nf,mm=10,ll=2)
      integer  lwork
      parameter(lwork= 10 + mm+2*ll+5*jmax+3*nf)
*      parameter(lwork= 1)
*     LOCAL
      integer   j,maxstep,method, maxnmv,order,testspace
      logical   wanted
      real      tau,tol,lock
      complex   alpha(jmax),beta(jmax), target
      complex   zwork(ndim,lwork)
      complex   tins(ndim)
      complex   res(ndim),tmp(ndim)
      complex   eivec(ndim,nf)
      logical   skip(nf)
      integer   id
      real      diag(ndim),eprec
      complex   rl(ndim)
      logical   basis

**<< HO
      real      err(nf)
**      include 'scale.com'
**      cs = 1
**      rs = 1
**>>
*
*...  Parameters for JDQZ
*
      tolabs = 1e-6
      tolred = 1e-6
      kmax = nf
      print*,"kmax = ",kmax," jmax = ",jmax," m = ",mm," lwork ",lwork 
      tau  = 1.0e+00
      target = cmplx(tau,0.0)
*      tau  = 0.05
*      target = cmplx(tau,0.45)
      tol    = 1.e-5 !1.e-6 !3.e-5
      maxstep=   50
      lock   = 1.e-7
      wanted = .true.
*     wanted=.true. : de eigenvectoren worden weggeschreven in eivec
*     wanted=.false.: de Schurvectoren worden weggeschreven in eivec
      order = 0
*     order =  0: nearest to target
*     order = -1: smallest real part
*     order =  1: largest real part
*     order = -2: smallest complex part
*     order =  2: largest complex part
      method = 1
*     method = 1: gmres(m)
*     method = 2: cgstab(l)
      maxnmv = 100    ! maximum number of matvecs in cgstab or gmres
      testspace = 3
*     Testspace 1: w = "Standard Petrov" 
*     Testspace 2: w = "Standard 'variable' Petrov" 
*     Testspace 3: w = "Harmonic Petrov"
*      
      tmp = 0
      id = 0
*LtR 22/11/01      call matrix(u,1.0)
      call matrix(u,tau)
      diag = 0.01
      call g_solve(diag,det,id)
*     eprec = 1.0e-03
*     call shuflmat
*     call Ascale(rl,diag)
*     call precon1(diag,eprec)
      call matrix(u,0.0)
*LtR 22/11/01      call matrix(u,tau)
**<< HO
**      call rcscaleAB(rs,cs)
**      skip = .false.
**      DO j=1,kmax
**       IF (skip(j)) cycle
**       IF ((abs(sig(j,2)).GT.1.0e-01).AND.(j.LT.nf)) THEN
**          eivec(:,j)=cmplx(w(:,j),w(:,j+1))
**          skip(j+1) = .true.
**       ELSE
**          eivec(:,j)=cmplx(w(:,j),0)
**       ENDIF
**      ENDDO
** >>     
      write(99,*) 'Stabjdqz begun'
      w   = 0.0
      sig = 0.0
      basis = .false.
!     call jdqz(alpha, beta, eivec, wanted, ndim, target, tol, kmax,
!    +     jmax, jmin, method, mm, ll, maxnmv, maxstep,
!    +     lock, order, testspace, basis)
!     call jdqz(alpha, beta, eivec, wanted, ndim, target, tol, kmax,
!    &     jmax, jmin, method, mm, ll, maxnmv, maxstep,
!    &     lock, order, testspace, basis,
!    &     zwork, lwork,tins )
!      call jdqz(alpha, beta, eivec, wanted, ndim, target, tol, kmax,
!    &     jmax, jmin, method, mm, ll, maxnmv, maxstep,
!    &     lock, order, testspace, basis,
!    &     zwork, lwork)
      call jdqz(alpha, beta, eivec, wanted, ndim, target, tol, kmax,
     &     jmax, jmin, method, mm, ll, maxnmv, maxstep,
     &     lock, order, testspace,
     &     zwork, lwork )

      write(99,*) 'jdqz done'
**<< HO
**      call scale_eivec(eivec,cs)
      if( kmax .eq. 0 ) write(99,*) 'jdqz did not converge'
**>>
      basis = .true.
      DO j = 1, kmax
        sig(j,1) = real(alpha(j)/beta(j))
        sig(j,2) = aimag(alpha(j)/beta(j))
        call amul(ndim,eivec(:,j),res) 
        call bmul(ndim,eivec(:,j),tmp)
        res = beta(j)*res - alpha(j)*tmp
        err(j) = real(sqrt(dot_product(res,res)) )
        write(99,999) sig(j,1),sig(j,2),err(j)
      END DO
      skip = .false.
      DO j=1,kmax
       IF (skip(j)) cycle
       IF ((abs(sig(j,2)).GT.1.0e-01).AND.(j.LT.nf)) THEN
         w(:,j)  = real(eivec(:,j))
         w(:,j+1)= aimag(eivec(:,j))
         skip(j+1) = .true.
       ELSE
         w(:,j)  = real(eivec(:,j))
       ENDIF
      ENDDO
      call sort(sig,w)
*
 999  format('sig= ',e12.4,'  +i ',e12.4,' res: ',e12.4)
      END
****************************************************************
      SUBROUTINE Amul(n,q,r)
*     r = A(q)
      implicit none
*     IMPORT/EXPORT
      integer,             intent(in) :: n
      complex,dimension(n),intent(in) :: q
      complex,dimension(n),intent(out):: r
*     LOCAL
      real,   dimension(n) :: rr,ri
*
      call matAvec(real(q) ,rr)
      call matAvec(aimag(q),ri)
      r = cmplx(rr,ri)
*
      END
****************************************************************
      SUBROUTINE Bmul(n,q,r)
*     r = B(q)
      implicit none
*     IMPORT/EXPORT
      integer,             intent(in) :: n
      complex,dimension(n),intent(in) :: q
      complex,dimension(n),intent(out):: r
*     LOCAL
      real     rr(n),ri(n)
*
      call matBvec(real(q) ,rr)
      call matBvec(aimag(q),ri)
      r = cmplx(rr,ri)
*
      END
****************************************************************
      SUBROUTINE precon(n,q)
      USE m_gsolve
*     Since there is (not yet) a preconditioning available
*     for complex matrices we will use a complete LU decom
*     from the NAG library.
*     r = A(q)
      implicit none
*     IMPORT/EXPORT
      integer  n
      complex  q(n)
*     LOCAL
      integer  id
      real     nrm,det
      real     qr(n),qi(n)
*
      id = 1
      qr = real(q)
      qi = aimag(q)
      nrm = sqrt(dot_product(qr,qr))
      IF (nrm.GT.1e-12) THEN
       call g_solve(qr,det,id)
      ELSE
       write(97,*) '|qr| lt eps'
      ENDIF
      nrm = sqrt(dot_product(qi,qi))
      IF (nrm.GT.1e-12) THEN
       call g_solve(qi,det,id)
      ELSE
       write(97,*) '|qi| lt eps'
      ENDIF
      q = cmplx(qr,qi)
*
      END
**************************************************************
      SUBROUTINE sort(s,v)
*     This routine sorts eigenvalues and vectors to their real parts.
      implicit none
      include 'bag.com'
*     IMPORT/EXPORT
      real    s(nf,2),v(ndim,nf)
*     LOCAL
      integer j,k,jmax
      real    s0(2),v0(ndim)
*
      DO k=1,nf-1
        s0(1) = s(k,1)
        s0(2) = s(k,2)
        jmax=k
        DO j=k+1,nf
          IF (real(s(j,1)).GT.real(s0(1))) THEN
            s0(1) = s(j,1)
            s0(2) = s(j,2)
            jmax=j
          ENDIF
        ENDDO
        v0(:) = v(:,jmax)
        s(jmax,1)=s(k,1)
        s(jmax,2)=s(k,2)
        s(k,1)=s0(1)
        s(k,2)=s0(2)
        v(:,jmax)=v(:,k)
        v(:,k)=v0(:)
      ENDDO
*
      END

