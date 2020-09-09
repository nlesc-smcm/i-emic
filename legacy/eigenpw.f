**********************************************************************
      SUBROUTINE backs(res,index)
      implicit none
      include 'bag.com'
*     IMPORT/EXPORT
      integer   index(nf)
      real      res(nf)
*     LOCAL     
      integer   LL
      real      er(ndim),ei(ndim)
      real      Awr(ndim),Awi(ndim),Bwr(ndim),Bwi(ndim)
*     FUNCTIONS
      real      l2nrm,linrm
*
      call matrix(u,0.0)
      LL = 1
 1000 IF (LL.LE.nf) THEN
*       (A-SIG*B*I)*W 
         IF (index(L).EQ.0) THEN
           call matAvec(w(:,LL),Awr)
           call matBvec(w(:,LL),Bwr)
           er = Awr - sig(LL,1)*Bwr
           res(LL)  = l2nrm(er,ndim)
           LL = LL+1
         ELSE
           call matAvec(w(:,LL),Awr)
           call matBvec(w(:,LL),Bwr)
           call matAvec(w(:,LL+1),Awi)
           call matBvec(w(:,LL+1),Bwi)
           er = Awr - ( sig(LL,1)*Bwr - sig(LL,2)*Bwi )
           ei = Awi - ( sig(LL,1)*Bwi + sig(LL,2)*Bwr )
           res(LL)  = l2nrm(er,ndim)
           res(LL+1)= l2nrm(ei,ndim)
           LL = LL+2
         ENDIF
         GOTO 1000
      ENDIF
      write(98,*) ' Re(sig)     Im(sig)     Err'
      write(98,999) (sig(LL,1),sig(LL,2),res(LL),LL=1,nf)
*
 999  format(3e12.4)
      END
*********************************************************
      SUBROUTINE linsolv(ui,uu,id)
      USE m_mat
      USE m_gsolve
      USE m_start
      implicit none
      include 'bag.com'
!      include 'mat.com'
!      include 'start.com'
*     IMPORT/EXPORT
      integer   id
      real      ui(ndim,nf),uu(ndim,nf)
*     LOCAL
      integer  L1,L2,ising,ninv
      real     coC(ndim)
*
      ising= 1
* LtR 05/10/01      ninv = 3
      ninv = 8
*LtR      tolabs = 1.0e-05
*LtR      tolred = 1.0e-04
*      tolabs = 1.0e-07
*      tolred = 1.0e-06
*ACdN SET tolabs/tolred in bgspars.F90
*
*     FILTERING BY INVERSE ITERATION: (A-B)u=Bu
*
      IF ((ising.EQ.1).AND.(id.EQ.0)) THEN
*     IF (ising.EQ.1) THEN
        call matrix(u,1.0)
        coC = coB
        DO L2=1,ninv
          DO L1=1,nf
            coB = coC
            call matBvec(ui(1,L1),uu(1,L1))
            call g_solve(uu(1,L1),det,id)
          ENDDO
          ui = uu
          write(97,*) 'inverse it: ',L2
        ENDDO
      ENDIF
*
*    END INVERSE ITERATION CYCLE
*
      call matrix(u,-1.0)
*
*     Generate right hand side: Multiply to obtain -(A+B)u
*
      DO L1=1,nf
        call matAvec(ui(1,L1),uu(1,L1))
      ENDDO
      uu = -uu
        
*
*     Solve (A-B)uu = -(A+B)uu
*
      call matrix(u,1.0)
      DO L1=1,nf
        call g_solve(uu(1,L1),det,id)
      ENDDO
      END
************************************************************
      SUBROUTINE eigen
      USE m_mat
      implicit none
      include 'bag.com'
!      include 'mat.com'
*     LOCAL
      integer  iz,index(nf),intger(nf)
      real     v(ndim,nf),uv(ndim,nf),z(1),ut(nf,ndim),
     +         h(nf,nf),g(nf,nf),b(nf,nf),p(nf,nf),
     +         rr(nf),ri(nf),d(nf)
      integer  i,k,L1,L2,j,kit,jeig,ifail,ipiv(nf),info,kmax,lmax
      real     res,eps,err(nf)
      integer id
      common /PRC/ id
*     FUNCTION
      integer  x02bhf,ib
      real     x02ajf
      real     l2nrm,linrm
	write(99,*) 'eigen:start'
*
      kmax = 8
*LtR 03/10/01      lmax = 5
      lmax = 10
      res=1.0e+00
      eps=1.0e-5
      iz =1
      ifail=0
      id =0
*
      DO kit=1,kmax
        IF (res.GT.eps) THEN
          DO L1=1,lmax
            index = 0
            DO i=1,nf
              call xnorm(index,i)
            ENDDO
            uv = w
            call linsolv(uv,w,id)
          ENDDO
	write(99,*) 'eigen: before reorientation'
*
*        REORIENTATION
*
          DO L1=1,nf
            DO i=1,ndim
              ut(L1,i)=uv(i,L1)
              v(i,L1) =w(i,L1)
            ENDDO
          ENDDO
          call f01ckf(G,UT,Uv,nf,nf,ndim,z,iz,1,ifail)
          call f01ckf(H,UT,V ,nf,nf,ndim,z,iz,1,ifail)
*
          call dgetrf(NF,NF,G,NF,IPIV,INFO)
          call dgetrs('N',NF,NF,G,NF,IPIV,H,NF,INFO)
          B = H
          ib = x02bhf()
          call f01atf(nf,ib,B,nf,k,l2,d)
          call f01akf(nf,k,l2,B,nf,intger)
          call f01apf(nf,k,l2,intger,B,nf,P,nf)
          call f02aqf(nf,k,l2,x02ajf(),B,nf,P,nf,rr,ri,intger,IFAIL)
          call f01auf(nf,k,l2,nf,d,P,nf)
*
          call sort(rr,ri,P,nf,index)
          call f01ckf(w,v,P,ndim,nf,nf,z,iz,1,ifail)
*
          jeig=1
	write(99,*) 'hallo'
 100      IF (jeig.LE.nf) THEN
            call xnorm(index,jeig)
            IF (index(jeig).EQ.1) THEN
              jeig=jeig+2
            ELSE
              jeig=jeig+1
            END IF
            GOTO 100
          END IF
          call xmap(rr,ri,sig,nf)
          call backs(err,index)
          res = l2nrm(err,2)
          write(99,999) kit,res
         ENDIF
      ENDDO 
      write(99,*) ' Re(sig)     Im(sig)     Err'
      write(99,998) (sig(k,1),sig(k,2),err(k),k=1,nf)
 999  format(' Iteration: ',i2,' Residue: ',e12.4)
 998  format(3e12.4)
      END
******************************************************************
      SUBROUTINE sort(rr,ri,p,nf,index)
*     This routine sorts (transformed) eigenvalues to their norms
*
      implicit none
*     IMPORT/EXPORT
      integer  nf,index(nf)
      real     rr(nf),ri(nf),p(nf,nf)
*     LOCAL
      integer  i,k,imax
      real     x0,y0,x,dum,eps
      eps=1.0e-06
      DO k=1,nf-1
         x0=rr(k)
         y0=ri(k)
         imax=k
         DO i=k+1,nf
            x=rr(i)*rr(i)+ri(i)*ri(i)
            IF (x.GT.(x0*x0+y0*y0)) THEN
               x0=rr(i)
               y0=ri(i)
               imax=i
            END IF
         ENDDO
         rr(imax)=rr(k)
         ri(imax)=ri(k)
         rr(k)=x0
         ri(k)=y0
         DO i=1,nf
            dum=p(i,k)
            p(i,k)=p(i,imax)
            p(i,imax)=dum
         ENDDO
      ENDDO
*
      DO i=1,nf
         IF (abs(ri(i)).GT.eps) THEN
            index(i)=1
         ELSE
            index(i)=0
         END IF
      ENDDO
      END
******************************************************************
      SUBROUTINE xnorm(index,jeig)
*     This routine normalizes the eigenvectors.
      implicit none
      include 'bag.com'
      integer  index(nf),jeig
*     LOCAL
      integer  i
      real     xmax,smax,e,f,s,emax,fmax,dum1,dum2
*     FUNCTIONS
      real     xrm
*
      xmax=0.0e+00
      smax=0.0e+00
      IF (index(jeig).EQ.0) THEN
         xmax=xrm(w(1,jeig),ndim,2)
      ELSE
         DO i=1,ndim
            e=w(i,jeig)
            f=w(i,jeig+1)
            s=abs(e*e+f*f)
            IF (s.GT.smax) THEN
               emax=e
               fmax=f
               smax=s
            END IF
         ENDDO
      END IF
      DO i=1,ndim
         IF (index(jeig).EQ.0) THEN
            dum1=w(i,jeig)/xmax
            w(i,jeig)=dum1
         ELSE
            dum1 = ( emax*w(i,jeig)+fmax*w(i,jeig+1) )/smax
            dum2 = (-fmax*w(i,jeig)+emax*w(i,jeig+1) )/smax
            w(i,jeig)  = dum1
            w(i,jeig+1)= dum2
         END IF
      ENDDO
      END
*******************************************************
      SUBROUTINE xmap(rr,ri,sig,nf)
*     This routine does the inverse Jennings transformation
*     sig=(r-1)/(r+1) to find the eigenvalues in the real domain
      implicit none
*     IMPORT/EXPORT
      integer  nf
      real     rr(nf),ri(nf),sig(nf,2)
*     LOCAL
      integer  i
      real     xnul,een,xl1p,xl12,xl22,xnrm

      een=1.0e+00
      DO i=1,nf
         xl1p=rr(i)+een
         xl12=rr(i)*rr(i)
         xl22=ri(i)*ri(i)
         xnrm=xl1p*xl1p+xl22
         sig(i,1)=(xl12+xl22-een)/xnrm
         sig(i,2)=2*ri(i)/xnrm
      ENDDO 
*
      END
*****************************************************************
      real FUNCTION xrm(u,ndim,iopt)
      implicit none
*     IMPORT/EXPORT
      integer  ndim,iopt
      real     u(ndim)
*     LOCAL
      integer  i
      real     x
*     FUNCTION
      real     f06ejf
*
*     IOPT=0 : MAX NORM
*     IOPT=1 : L1 NORM
*     IOPT=2 : L2 NORM
*
      x=0.0e+00
      IF (iopt.EQ.0) THEN
         DO i=1,ndim
            IF (abs(u(i)).GT.x) THEN
               x=abs(u(i))
            END IF
         ENDDO
         xrm=x
      END IF
      IF (iopt.EQ.1) THEN
         DO i=1,ndim
            x=x+abs(u(i))
         ENDDO 
         xrm=x
      END IF
      IF (IOPT.EQ.2) THEN
         xrm=f06ejf(ndim,u,1)
      END IF
      END
            
            
