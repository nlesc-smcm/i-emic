*********************************************
* NB Subr. thdepth: thermocline depth taken
* at m-11 instead of m-8!!
*********************************************
      SUBROUTINE monitor(itp,ifile)
      implicit none
      integer nmt
      parameter(nmt=2)
      include 'bag.com'
*     IMPORT/EXPORT
      integer  itp,ifile,i
*     LOCAL
      real,    parameter :: pi   =  3.14159265358979323846
      real mheat(nmt),msalt(nmt),period,growth
      real timesc
      real eh,ahe,psimax
      real dep, depk,buof
*     COMMON
      integer ist,iend,jst,jend
      common /integrate/ ist,iend,jst,jend
*
      eh = par(EK_H)
      call transport(u,eh,psimax,mheat,msalt,nmt,ahe,timesc)
      growth = (sig(1,1)/timesc) * (24.*3600.*365.)
      IF (abs(sig(1,2)).gt.1.0e-06) THEN
      period = (2*pi/sig(1,2))*timesc/(24.*3600.*365.)
      ELSE
      period = 0.0
      END IF
      call thdepth(u,dep)
      call buoforc(u,buof)
      write(ifile,999) itp,xl,det,ahe,growth,period,
     +       mheat(1),mheat(2),msalt(1),msalt(2),psimax,
     +       dep, buof
*
* Define region of integration
*
      ist  = 1
      iend = n
      jst  = 1
      jend = m
*
* LtR 22/10/01      call kinetic(u,psimax,xl)
*
* LtR 22/10/01      call potential(u,psimax,xl)
*
 999  format(i5,12e12.4)
      END
********************************************************************
      SUBROUTINE monitor_eigen
      implicit none
      include 'bag.com'
*     LOCAL
      real psimax,uacc(ndim)
      integer nvec,i,k
*     COMMON
      integer ist,iend,jst,jend
      common /integrate/ ist,iend,jst,jend
*
* Define region of integration
*
      ist  = 1
      iend = n
      jst  = 1
      jend = m
*
      call streamf(u,psimax)
*
* Determine position of most unstable eigenvector
*
      do i=1,nf
         if (abs(sig(i,1)).gt.1e-10) then
            nvec = i
            goto 776
         endif
      enddo
 776  continue
      write(99,*) 'nvec = ',nvec,' sig = ',sig(nvec,1)
*
      uacc = 0.0
      do i=1,ndim
        uacc(i) = w(i,nvec)
      enddo
*
      call pot_exch(u,uacc,psimax,xl)
*
      END
********************************************************************
      SUBROUTINE monitor2(un,ifile)
      implicit none
      integer nmt
      parameter(nmt=2)
      include 'bag.com'
*     IMPORT/EXPORT
      integer  ifile,un(ndim)
*     LOCAL
      real eh,ahe,psimax,i
      real dep, buof
      real depk,dr,drh,depk2,dr2,dr3,dr4,dr5,dr6
      real mheat(nmt),msalt(nmt),period,growth
      real timesc
*
      eh = par(EK_H)
*
      call transport(un,eh,psimax,mheat,msalt,nmt,ahe,timesc)
      call thdepth(un,dep)
      call buoforc(un,buof)
      call drho(un,depk,depk2,dr,drh,dr2,dr3,dr4,dr5,dr6)
*      write(ifile,996)
*      write(ifile,999) (par(i),i=1,npar)
      write(ifile,996)
      write(ifile,997) icp
*      write(ifile,992) xl,psimax,dep,depk,buof,dr,drh
      write(ifile,992) xl,psimax,dep,depk,buof,dr,drh,depk2,
     &                 dr2,dr3,dr4,dr5,dr6
*
 999  format('*',5e15.5)
 997  format('*',3x,'par(',i2,')',4x,'Psimax',5x,
     &       'Th depth', 4x, 'Thd Klinger',2x,
     &       'wB', 9x,'T_n-T_s' ,2x, 'T_n-T_s*dz*dz')
 996  format('*')
 992  format(13e12.4)
      END
********************************************************************
      SUBROUTINE transport(un,eh,psimax,mheat,msalt,nmt,ahe,timesc)
*     write solution
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      integer  nmt
      real     un(ndim)
      real     mheat(nmt),msalt(nmt),psimax
      real     eh,ahe,timesc
*     LOCAL
      integer i,j
      real   u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1)
      real   w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1)
      real   t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      real   transc,heatsc,saltsc
*
      timesc = r0dim/udim
      transc = r0dim*hdim*udim
      heatsc = rhodim*udim*cp0*hdim*r0dim  ! heat
      saltsc = udim*hdim*r0dim             ! salt
      ahe = eh * (2*omegadim*r0dim*r0dim)
      call streamf(un,psimax)
*  Compute overturning transport in Sv
      psimax = psimax * transc/1.0e+06
      call heatsalt(un,mheat,msalt,nmt)
      mheat = heatsc * mheat/1.0e+15
      msalt = saltsc * msalt/1.0e+06
      END
********************************************************************************
      SUBROUTINE streamf(un,psimax)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real   psimax,un(ndim)
*     LOCAL
      real  u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +      w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +      t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      integer    i,j,k
      real       vs(0:m,l),psim(0:m,0:l),cs
      real       dum
*
*     Calculation meridional overturning streamfunction
*
      call usol(un,u,v,w,p,t,s)
      vs = 0.0
      DO j=0,m
         DO k=1,l
            dum = 0.0
            DO i=1,n
               dum=v(i,j,k)*dx+dum
            ENDDO
            vs(j,k) = dum
         ENDDO
      ENDDO
      psim(:,0) = 0.0
      psim(0,:) = 0.0
      DO j=1,m
         cs = cos(yv(j))
         DO k=1,l
            psim(j,k)=cs*vs(j,k)*dz*dfzT(k)  + psim(j,k-1)
         ENDDO
      ENDDO
*      psimax = maxval(abs(psim))
*      psimax = maxval(psim)-minval(psim);
      psimax = 0
      dum = 0
      do j = 0, m
         do k = 0, l
            if( psim(j,k) > psimax ) psimax = psim(j,k)
            if( psim(j,k) < dum    ) dum    = psim(j,k)
         enddo
      enddo
      psimax = psimax - dum
* writeout separate file for psim
      call idl2out(psim,psimax,15)
      END
*****************************************************************
      SUBROUTINE maxminw(un,xl)
      implicit none
      include 'usr.com'
      real   un(ndim),xl
*     LOCAL
      real  u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +      w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +      t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      real wmax, wmin
*
      wmax=0.0
      wmin=0.0
      call usol(un,u,v,w,p,t,s)
      wmax=maxval(w)
      wmin=minval(w)
*
*   case 10N-74N: j=5 and m/2 correspond to 30N and 42N
*   case 54S-74N: j=21 and 24 correspond to 30N and 42N
*   both cases: 30N -> j=m-11; 42N -> m-8
*
*      write(88,992) xl,wmax,wmin,w(n/2,m-11,12),w(n/2,m-8,12),
*     +            w(n/4,m-11,12),w(n/4,m-8,12)
 992  format(13e12.4)
*
      END

*******************************************************************************
      SUBROUTINE heatsalt(un,mheat,msalt,nmt)
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      integer  nmt
      real     un(ndim)
      real     mheat(nmt),msalt(nmt)
*     LOCAL
      real  u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +      w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +      t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      integer i,j,k
      real   hfu(n,m,l),hfv(n,m,l),sfu(n,m,l),sfv(n,m,l)
* fluxes
      real   thz(m,l),thzm(n,l),thzvi(n),
     +       shz(m,l),shzm(n,l),shzvi(n)
* zonal transport
      real   thm(n,l),thmvi(m),thmn(m,l),
     +       shm(n,l),shmvi(m),shmn(m,l)
* meridional transport
      real   tdx,tdy,cs,phinv
*
      phinv = kappah/(r0dim*udim)
      call usol(un,u,v,w,p,t,s)
*
* compute dimensionless heat flux and salt flux vector
*
      tdx = 2*dx
      tdy = 2*dy
      hfv = 0.0
      hfu = 0.0
      sfv = 0.0
      sfu = 0.0
      DO k=1,l
       DO i=1,n
         DO j=1,m
            cs = cos(y(j))
            hfv(i,j,k) = 0.5*(v(i,j,k) + v(i,j-1,k)) * t(i,j,k) -
     +                 phinv*(t(i,j+1,k) - t(i,j-1,k))/tdy
            sfv(i,j,k) = 0.5*(v(i,j,k) + v(i,j-1,k)) * s(i,j,k) -
     +                 phinv*(s(i,j+1,k) - s(i,j-1,k))/tdy
            hfu(i,j,k) = 0.5*(u(i,j,k) + u(i-1,j,k)) * t(i,j,k) -
     +                 phinv*(t(i+1,j,k) - t(i-1,j,k))/(tdx*cs)
            sfu(i,j,k) = 0.5*(u(i,j,k) + u(i-1,j,k)) * s(i,j,k) -
     +                 phinv*(s(i+1,j,k) - s(i-1,j,k))/(tdx*cs)
         ENDDO
       ENDDO
      ENDDO
* zonal heat and salt transport
      DO k=1,l
       DO i=1,n
        thz(1,k) = hfu(i,1,k)
        shz(1,k) = sfu(i,1,k)
        DO j=2,m
           thz(j,k) = thz(j-1,k) + hfu(i,j,k)
           shz(j,k) = shz(j-1,k) + sfu(i,j,k)
        ENDDO
      	thzm(i,k) = thz(m,k)*dy
      	shzm(i,k) = shz(m,k)*dy
       ENDDO
      ENDDO
* meridional heat and salt transport
      DO k=1,l
       DO j=1,m
        cs = cos(y(j))
        thm(1,k) = cs*hfv(1,j,k)
        shm(1,k) = cs*sfv(1,j,k)
        DO i=2,n
           thm(i,k) = thm(i-1,k)+cs*hfv(i,j,k)
           shm(i,k) = shm(i-1,k)+cs*sfv(i,j,k)
        ENDDO
      	thmn(j,k) = thm(n,k)*dx
      	shmn(j,k) = shm(n,k)*dx
       ENDDO
      ENDDO
* vertical integration
      thmvi = 0.0
      shmvi = 0.0
      DO j=1,m
       DO k=1,l
          thmvi(j) = thmn(j,k) * dz *dfzT(k) + thmvi(j)
          shmvi(j) = shmn(j,k) * dz *dfzT(k) + shmvi(j)
       ENDDO
      ENDDO
      mheat(1) = thmvi(m/4)
      mheat(2) = thmvi(m/2)
      msalt(1) = shmvi(m/4)
      msalt(2) = shmvi(m/2)
      END
******************************************************************
      SUBROUTINE thdepth(un,dep)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      dep,un(ndim)
*     LOCAL
      real      u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +          w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +          t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      integer   i,j,k
      real      rho(0:n+1,0:m+1,0:l+1)
      real      lambda
      real      rhoe(n,m), depf(n,m), depm(m)
*
*     Calculation thermocline depth depf(i,j)
*
      call usol(un,u,v,w,p,t,s)
      lambda = par(LAMB)
      rho = lambda*s - t
*
      DO i=1,n
        DO j=1,m
          rhoe(i,j)=-0.9*(rho(i,j,l)-rho(i,j,1))+rho(i,j,l)
        ENDDO
      ENDDO
*
      DO i=1,n
        DO j=1,m
          depf(i,j)=0.0
          DO k=1,l
            IF ((rho(i,j,k+1).LE.rhoe(i,j)).AND.
     +            (rho(i,j,k).GT.rhoe(i,j))) THEN
              depf(i,j)=zmax-(dz *dfzW(k)*(rhoe(i,j)-
     +           rho(i,j,k))/(rho(i,j,k+1)-rho(i,j,k))
     +           + z(k))
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*
* take zonal average of thermocline depth field depf(i,j)
*
      depm=0.0
      DO j=1,m
         DO i=1,n
           depm(j) = depf(i,j)/n + depm(j)
         ENDDO
      ENDDO
*
* take as a measure the thermocline depth at 42N
* FOR BASIN FROM 10 TO 74 j = m-8 = 8 GIVES 42 N
* FOR BASIN FROM -22 TO +74 j = m-8 = 16 GIVES 42 N
* FOR BASIN FROM -54 TO +74 j = m-8= 24 GIVES 42 N
*
c$$$      dep=depm(m-11)
c$$$      write(24,*) "dep normal= ", dep
*
* For testing: print depth
*
*      DO i=1,n
*        DO j=1,m
*            write(24,*) i,j,depf(i,j)
*        ENDDO
*      ENDDO
*
*      write(24,*) "*******"
*        DO j=1,m
*            write(24,*) j,depm(j)
*        ENDDO
*
*      write(24,*) "*******"
*      DO k=1,l
*         write(24,*) z(k)
*      ENDDO
*
      END
******************************************************************
      SUBROUTINE drho(un,depk,depk2,dr,drh,dr2,dr3,dr4,dr5,dr6)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      dr,drh,depk,un(ndim),depk2
      real      dr2,dr3,dr4,dr5,dr6
*     LOCAL
      real      u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +          w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +          t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      real      rho(0:n+1,0:m+1,0:l+1)
      integer   i,j,k
      real      lambda,cs
      real      tzo(m,l),drl(l),zth
      real      tsum,tzsum,tzonk,tzon(l+1)
      integer   k0
*
      call usol(un,u,v,w,p,t,s)
      lambda = par(LAMB)
      rho = lambda*s - t
*
* Thermocline depth Klinger et al ('99)
*
      tzon=0.0
      DO k=1,l+1
          DO i=1,n
*            tzon(k)=tzon(k) + t(i,3,k)/n
            tzon(k)=tzon(k) + t(i,m/2,k)/n
          ENDDO
      ENDDO
      k0=0.0
      tzonk=-0.99*(tzon(l)-tzon(1))+tzon(l)
      DO k=1,l
        IF ((tzon(k+1).GT.tzonk).AND.(tzon(k).LT.tzonk)) k0=k
      ENDDO
      zth=z(k0+1) - (z(k0+1)-z(k0))*(tzon(k0+1)-tzonk)/
     +    (tzon(k0+1)-tzon(k0))
*
      tsum=0.0
      tzsum=0.0
      DO k=k0+2,l
         tsum=tsum + tzon(k)*dz *dfzT(k)
         tzsum=tzsum + tzon(k)*(-z(k))*dz *dfzT(k)
*         tzsum=tzsum + tzon(k)*(1.0+z(k))*dz *dfzT(k)
      ENDDO
      tsum = tsum +  tzon(k0+1)*(dz*dfzT(k)+(zw(k0)-zth))
      tzsum = tzsum + tzon(k0+1)*(-z(k0+1))*(dz*dfzT(k)+(zw(k0)-zth))
*      tzsum = tzsum + tzon(k0+1)*(1.0+z(k0+1))*(dz*dfzT(k)+(zw(k0)-zth))
*      DO k=k0+1,l
*         tsum=tsum + tzon(k)*dz*dfzT(k)
*         tzsum=tzsum + tzon(k)*(1.0+z(k))*dz*dfzT(k)
*      ENDDO
*
      depk=tzsum/tsum
*      write(24,*) "dep normal= ", dep
      write(24,*) "k0= ", k0
      write(24,*) "z(k0)=", z(k0)
      write(24,*) "zth= ", zth
      write(24,*) "z(k0+1)=", z(k0+1)
      write(24,*)  "z(k)="
      DO k=1,l
      write(24,*) z(k)
      ENDDO
      write(24,*) "depk= ", depk
      write(24,*) "tzsum= " ,tzsum
      write(24,*) "tsum= " ,tsum
      write(24,*) "tzon(k)                  1+z(k) "
      DO k=1,l
      write(24,*) tzon(k) , 1.0+z(k)
      ENDDO
      write(24,*) "tzonk=", tzonk
*
      tsum=0.0
      tzsum=0.0
      DO k=k0+2,l
         tsum=tsum + tzon(k)*dz *dfzT(k)
         tzsum=tzsum + tzon(k)*(1.0+z(k))*dz*dfzT(k)
      ENDDO
      tsum = tsum +  tzon(k0+1)*(dz*dfzT(k)+(zw(k0)-zth))
      tzsum = tzsum + tzon(k0+1)*(1.0+z(k0+1))*
     &        (dz*dfzT(k)+(zw(k0)-zth))
*
      depk2=tzsum/tsum
*
*     Other measures for north-south temp. diff:
*
*   Difference between near-surface temperature of two
*   latitudinal strips
*
      tzo=0.0
      DO k=1,l
        DO j=1,m
          DO i=1,n
            tzo(j,k)= t(i,j,k)/n+tzo(j,k)
          ENDDO
        ENDDO
        drl(k)=(tzo(2,k)-tzo(m-1,k))
*        drl(k)=-(tzo(m-2,k)+tzo(m-1,k)-tzo(2,k)-tzo(3,k))/2
      ENDDO
      dr=drl(l-1)
      dr2=(tzo(1,l-1)-tzo(m,l-1))
      dr3=(tzo(1,l)-tzo(m,l))
      dr4=(tzo(2,l)-tzo(m-1,l))
      dr5=(t(n/2,1,l)-t(n/2,m,l))
*
* For basin 54S-74N, drho 0N - 74N (Tmax-Tmin)
*
      dr6=(tzo(14,l)-tzo(m,l))
*
      write(24,*) 'dr= ', dr
      write(24,*) 'drl(k)='
*      DO k=1,l
*      write(24,*) drl(k)
*      ENDDO
*      write(24,*) 'tzo(2,10)=', tzo(2,10),'tzo(14,10)=',tzo(14,10)
*      write(24,*) 't(i,2,10)=         ' ,'t(i,14,10)='
*      DO i=1,n
*      write(24,*) t(i,2,10), t(i,14,10)
*      ENDDO
*
*    Double vertical integral of this quantity, to imitate
*    depth-integrated steric height
*
*
      tzonk=0.0
      tzonk=-0.9*(tzon(l)-tzon(1))+tzon(l)
      DO k=1,l
        IF ((tzon(k+1).GT.tzonk).AND.(tzon(k).LT.tzonk)) k0=k
      ENDDO
      zth=0.0
      zth=z(k0+1) - (z(k0+1)-z(k0))*(tzon(k0+1)-tzonk)/
     +    (tzon(k0+1)-tzon(k0))
*
* take double vertical integral over thermocline
*
      tsum=0.0
      tsum = tsum +  tzon(k0+1)*(dz*dfzT(k)+(zw(k0)-zth))
      drh=0.0
      DO k=k0+2,l
        drh=drl(k)*dz*dz*dfzT(k)*dfzT(k)+drh
      ENDDO
      drh = drh +drl(k0+1)*(dz*dfzT(k)+(zw(k0)-zth))*
     &      (dz*dfzT(k)+(zw(k0)-zth))
      write(24,*) 'k0 (0.9) =', k0
      write(24,*) 'drh=', drh
*
      END
******************************************************************
      SUBROUTINE buoforc(un,buof)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      buof,un(ndim)
*     LOCAL
      real      u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +          w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +          t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      real      rho(0:n+1,0:m+1,0:l+1)
      integer   i,j,k
      real      lambda, ra,cs
      real      bbb(n,m,l)
*
*     Calculation buoyancy forcing wB
*
      call usol(un,u,v,w,p,t,s)
      lambda = par(LAMB)
      rho = lambda*s - t
      ra = par(RAYL)
      bbb = 0.0
      buof = 0.0
*
      DO i=1,n
        DO j=1,m
          DO k=1,l
            bbb(i,j,k)=-(w(i,j,k-1)+w(i,j,k))*ra*rho(i,j,k)
          ENDDO
        ENDDO
      ENDDO
*
* Take volume integral
*
      DO i=1,n
        DO j=1,m
         cs = cos(y(j))
          DO k=1,l
             buof=bbb(i,j,k)*dx*cs*dy*dz*dfzT(k)+buof
          ENDDO
        ENDDO
      ENDDO
*
      END
******************************************************************
      SUBROUTINE termprep(un)
      implicit none
      include 'usr.com'
*     CONSTANT
      integer  ntermu, ntermv, ntermw, ntermp,ntermt,
     +         nterms,dimvect
      parameter(ntermu=4, ntermv=4, ntermw=2, ntermp=3,
     +          ntermt=6, nterms=6, dimvect=l*(ntermu+ntermv+
     +          ntermw+ntermp+ntermt+nterms))
*     IMPORT/EXPORT
      integer  ii,jj
      real     vecterm(dimvect)
      real     un(ndim)
*     LOCAL
      integer i,j,k
*
      i=n/4
      j=m/4
      call terms(un,vecterm,i,j)
      call idl4out(vecterm,i,j,31)

      END
**************************************************************
      SUBROUTINE idl4out(f,i,j,nout)
      implicit none
      include 'usr.com'
*     CONSTANT
      integer  ntermu, ntermv, ntermw, ntermp,ntermt,
     +         nterms,dimvect
      parameter(ntermu=4, ntermv=4, ntermw=2, ntermp=3,
     +          ntermt=6, nterms=6, dimvect=l*(ntermu+ntermv+
     +          ntermw+ntermp+ntermt+nterms))

*     IMPORT/EXPORT
      integer   i,j,nout
      real      f(dimvect)
*     LOCAL
      integer   p,k,p1,p2,p3,p4,p5
      real      maxf,eps
*
*   schrijf eerst i,j en l weg, dan xh en yh, dan zh en z
*   en dan de data
*
      eps=1.0e-10
      maxf =-10.
*
      p1=l*ntermu
      p2=l*(ntermu+ntermv)
      p3=l*(ntermu+ntermv+ntermw)
      p4=l*(ntermu+ntermv+ntermw+ntermp)
      p5=l*(ntermu+ntermv+ntermw+ntermp+ntermt)
      write(nout,999) i,j,l
      write(nout,997) x(i),y(j)
      DO p=1,p1
         maxf=max(abs(f(p)),maxf)
      ENDDO
      write(nout,998) maxf
      maxf =-10.
      DO p=p1+1,p2
         maxf=max(abs(f(p)),maxf)
      ENDDO
      write(nout,998) maxf
      maxf =-10.
      DO p=p2+1,p3
         maxf=max(abs(f(p)),maxf)
      ENDDO
      write(nout,998) maxf
      maxf =-10.
      DO p=p3+1,p4
         maxf=max(abs(f(p)),maxf)
      ENDDO
      write(nout,998) maxf
      maxf =-10.
      DO p=p4+1,p5
         maxf=max(abs(f(p)),maxf)
      ENDDO
      write(nout,998) maxf
      maxf =-10.
      DO p=p5+1,dimvect
         maxf=max(abs(f(p)),maxf)
      ENDDO
      write(nout,998) maxf
      DO k=1,L
         write(nout,998) zw(k)
      ENDDO
      DO k=1,L
         write(nout,998) z(k)
      ENDDO
      DO p=1,dimvect
         write(nout,998) f(p)
      ENDDO
 999  format(3i8)
 998  format(e18.10)
 997  format(2e18.10)
      END
******************************************************************
      SUBROUTINE terms(un,vecterm,ii,jj)
      implicit none
      include 'usr.com'
*     CONSTANT
      integer  ntermu, ntermv, ntermw, ntermp,ntermt,
     +         nterms,dimvect
      parameter(ntermu=4, ntermv=4, ntermw=2, ntermp=3,
     +          ntermt=6, nterms=6, dimvect=l*(ntermu+ntermv+
     +          ntermw+ntermp+ntermt+nterms))
*     IMPORT/EXPORT
      integer  ii,jj
      real     vecterm(dimvect)
      real     un(ndim)
*     EXTERNAL
      real      hs
*     LOCAL
      real      u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +          w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +          t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      real      rho(0:n+1,0:m+1,0:l+1)
      integer   i,j,k,l1,l2,l3,l4,l5,l6
      real      eh,ev,ra,lambda,ph,pv,pvc,rota
      real      sn(m),cs(m),tn(m),snv(m),csv(m),tnv(m)
      real      xcor(l),xpres(l),fruh1(l),fruh2(l),fruh3(l),
     +          fruh4(l),fruh5(l),fruh(l),fruv(l),
     +          ycor(l),ypres(l),frvh1(l),frvh2(l),frvh3(l),
     +          frvh4(l),frvh5(l),frvh(l),frvv(l),
     +          zpres(l),buow(l),dudx(l),dvdy(l),dwdz(l),
     +          xadt(l),yadt(l),zadt(l),frth1(l),frth2(l),
     +          frth(l),frtv(l),cadt(l),
     +          xads(l),yads(l),zads(l),frsh1(l),frsh2(l),
     +          frsh(l),frsv(l),cads(l)
      integer   ind
      real      h1,h2
*
      xcor=0.0
      xpres=0.0
      fruh1=0.0
      fruh2=0.0
      fruh3=0.0
      fruh4=0.0
      fruh5=0.0
      fruh=0.0
      fruv=0.0
      DO j=1,m
         sn(j)=sin(y(j))
         cs(j)=cos(y(j))
         tn(j)=tan(y(j))
         snv(j)=sin(yv(j))
         csv(j)=cos(yv(j))
         tnv(j)=tan(yv(j))
      ENDDO
*
      eh     = par(EK_H)
      ev     = par(EK_V)
      ra     = par(RAYL)
      lambda = par(LAMB)
      ph     = par(PE_H)
      pv     = par(PE_V)
      pvc    = par(P_VC)
      rota   = 0.0 !par(ROTN) ! ATvS-Mix
*
      i=ii
      j=jj
*
      call usol(un,u,v,w,p,t,s)
*
      rho = lambda*s - t
*
      DO k=1,l
         h1 = 1./(dfzT(k)*dfzW(k))
         h2 = 1./(dfzT(k)*dfzW(k-1))
          xcor(k)  = 0.25*rota*sn(j)*(v(i,j,k)+v(i+1,j,k)+v(i,j-1,k)+
     &               v(i+1,j-1,k))
          xpres(k) = -(p(i+1,j,k)-p(i,j,k))/(cs(j)*dx)
           fruh1(k) = (u(i+1,j,k)+u(i-1,j,k)-2*u(i,j,k))/
     &                (cs(j)*cs(j)*dx*dx)
           fruh2(k) =-u(i,j,k)/(cs(j)*cs(j))
           fruh3(k) =-tn(j)*(u(i,j+1,k)-u(i,j-1,k))/(2*dy)
           fruh4(k) = (u(i,j+1,k)+u(i,j-1,k)-2*u(i,j,k))/(dy*dy)
           fruh5(k) =-tn(j)*(v(i+1,j,k)-v(i,j,k)+v(i+1,j-1,k)-
     &                v(i,j-1,k))/(cs(j)*dx)
          fruh(k)  = eh*(fruh1(k)+fruh2(k)+fruh3(k)+fruh4(k)+fruh5(k))
*          fruv(k)  = ev*(u(i,j,k+1)+u(i,j,k-1)-2*u(i,j,k))/
*     &               (dz*dz*dfzT(k)*dfzT(k))
          fruv(k)  = ev*(h1*u(i,j,k+1)+h2*u(i,j,k-1)-(h1+h2)*u(i,j,k))/
     &               (dz*dz)
      ENDDO
*
      DO k=1,l
         L1 = ntermu*(k-1)
         vecterm(L1+1) = xcor(k)
         vecterm(L1+2) = xpres(k)
         vecterm(L1+3) = fruh(k)
         vecterm(L1+4) = fruv(k)
      ENDDO
*
      write(27,*) 'x-equation:'
      DO k=1,l
      write(27,919) xcor(k),xpres(k),fruh(k),fruv(k)
      ENDDO
      write(29,*) 'x-equation:'
      DO k=1,L
      write(29,909) xcor(k)+xpres(k)+fruh(k)+fruv(k)
      ENDDO
      write(28,*) 'x-equation:'
      DO ind=1,dimvect
      write(28,909) vecterm(ind)
      ENDDO
*
      DO k=1,l
         h1 = 1./(dfzT(k)*dfzW(k))
         h2 = 1./(dfzT(k)*dfzW(k-1))
          ycor(k)  = -0.25*rota*snv(j)*(u(i,j,k)+u(i,j+1,k)+
     &               u(i-1,j,k)+u(i-1,j+1,k))
          ypres(k) = -(p(i,j+1,k)-p(i,j,k))/dy
           frvh1(k) = (v(i+1,j,k)+v(i-1,j,k)-2*v(i,j,k))/
     &                (csv(j)*csv(j)*dx*dx)
           frvh2(k) =-v(i,j,k)/(csv(j)*csv(j))
           frvh3(k) =-tnv(j)*(v(i,j+1,k)-v(i,j-1,k))/(2*dy)
           frvh4(k) = (v(i,j+1,k)+v(i,j-1,k)-2*v(i,j,k))/(dy*dy)
*           frvh5(k) = tnv(j)*(u(i,j,k)-u(i-1,j,k)+u(i,j+1,k)-
*     &                u(i-1,j+1,k))/(csv(j)*dx)
           frvh5(k) = tn(j)*(u(i,j,k)-u(i-1,j,k)+u(i,j+1,k)-
     &                u(i-1,j+1,k))/(cs(j)*dx)
          frvh(k)  = eh*(frvh1(k)+frvh2(k)+frvh3(k)+frvh4(k)+frvh5(k))
*          frvv(k)  = ev*(v(i,j,k+1)+v(i,j,k-1)-2*v(i,j,k))/
*     &               (dz*dz*dfzT(k)*dfzT(k))
          frvv(k)  = ev*(h1*v(i,j,k+1)+h2*v(i,j,k-1)-(h1+h2)*v(i,j,k))/
     &               (dz*dz)
      ENDDO
*
      DO k=1,l
         L2 = ntermu*l+ntermv*(k-1)
         vecterm(L2+1) = ycor(k)
         vecterm(L2+2) = ypres(k)
         vecterm(L2+3) = frvh(k)
         vecterm(L2+4) = frvv(k)
      ENDDO
*
      write(27,*) 'y-equation:'
      DO k=1,l
      write(27,919) ycor(k),ypres(k),frvh(k),frvv(k)
      ENDDO
      write(29,*) 'y-equation:'
      DO k=1,L
      write(29,909) ycor(k)+ypres(k)+frvh(k)+frvv(k)
      ENDDO
*
      DO k=1,l
          zpres(k) = -(p(i,j,k+1)-p(i,j,k))/(dz*dfzW(k))
          buow(k)  = ra*(t(i,j,k+1)+t(i,j,k)-
     &               lambda*(s(i,j,k+1)+s(i,j,k)))/2.
      ENDDO
*
      DO k=1,l
         L3 = ntermu*l+ntermv*l+ntermw*(k-1)
         vecterm(L3+1) = zpres(k)
         vecterm(L3+2) = buow(k)
      ENDDO
*
      write(27,*) 'z-equation:'
      DO k=1,l
      write(27,919) zpres(k),buow(k)
      ENDDO
      write(29,*) 'z-equation:'
      DO k=1,L
      write(29,909) zpres(k)+buow(k)
      ENDDO
*
      DO k=1,l
          dudx(k) = (u(i,j,k)-u(i-1,j,k))/(cs(j)*dx)
          dvdy(k) = (v(i,j,k)*csv(j)-v(i,j-1,k)*csv(j-1))/(cs(j)*dy)
          dwdz(k) = (w(i,j,k)-w(i,j,k-1))/(dz*dfzT(k))
      ENDDO
*
      DO k=1,l
         L4 = ntermu*l+ntermv*l+ntermw*l+ntermp*(k-1)
         vecterm(L4+1) = dudx(k)
         vecterm(L4+2) = dvdy(k)
         vecterm(L4+3) = dwdz(k)
      ENDDO
*
      write(27,*) 'p-equation:'
      DO k=1,l
      write(27,919) dudx(k),dvdy(k),dwdz(k)
      ENDDO
      write(29,*) 'p-equation:'
      DO k=1,L
      write(29,909) dudx(k)+dvdy(k)+dwdz(k)
      ENDDO
*
      DO k=1,l
        h1 = 1.0/(dfzT(k)*dfzW(k))
        h2 = 1.0/(dfzT(k)*dfzW(k-1))
          xadt(k) =-(u(i,j,k)*(t(i+1,j,k)+t(i,j,k))/(2*dx*cs(j))-
     &              u(i-1,j,k)*(t(i,j,k)+t(i-1,j,k))/(2*dx*cs(j)))
          yadt(k) =-(v(i,j,k)*(t(i,j+1,k)+t(i,j,k))*csv(j)/(cs(j)*2*dy)-
     &              v(i,j-1,k)*(t(i,j,k)+t(i,j-1,k))*csv(j-1)/
     &              (cs(j)*2*dy))
*          zadt(k) =-(w(i,j,k)*(t(i,j,k+1)+t(i,j,k))/(2*dz)-
*     &              w(i,j,k-1)*(t(i,j,k)+t(i,j,k-1))/(2*dz))
          zadt(k) =-(w(i,j,k)*(t(i,j,k+1)+t(i,j,k))/(2*dz*dfzT(k))-
     &              w(i,j,k-1)*(t(i,j,k)+t(i,j,k-1))/(2*dz*dfzT(k)))
           frth1(k) = (t(i+1,j,k)+t(i-1,j,k)-2*t(i,j,k))/
     &                (cs(j)*cs(j)*dx*dx)
           frth2(k) = ((t(i,j+1,k)-t(i,j,k))*csv(j)-
     &                (t(i,j,k)-t(i,j-1,k))*csv(j-1))/
     &                (cs(j)*dy*dy)
          frth(k) = ph*(frth1(k)+frth2(k))
*          frtv(k) = pv*(t(i,j,k+1)+t(i,j,k-1)-2*t(i,j,k))/(dz*dz)
          frtv(k) = pv*(h1*t(i,j,k+1)+h2*t(i,j,k-1)-
     &             (h1+h2)*t(i,j,k))/(dz*dz)
*          cadt(k) = pvc*(hs(rho(i,j,k+1)-rho(i,j,k))*
*     &              (t(i,j,k+1)-t(i,j,k))-hs(rho(i,j,k)-rho(i,j,k-1))*
*     &              (t(i,j,k)-t(i,j,k-1)))/(dz*dz)
          cadt(k) = pvc*(h1*hs(rho(i,j,k+1)-rho(i,j,k))*t(i,j,k+1)
     &              + h2*hs(rho(i,j,k)-rho(i,j,k-1))*t(i,j,k-1)
     &              - (h1+h2)*( hs(rho(i,j,k+1)-rho(i,j,k))+
     &              hs(rho(i,j,k)-rho(i,j,k-1))*t(i,j,k) ))/(dz*dz)
      ENDDO
*
      DO k=1,l
         L5 = ntermu*l+ntermv*l+ntermw*l+ntermp*l+ntermt*(k-1)
         vecterm(L5+1) = xadt(k)
         vecterm(L5+2) = yadt(k)
         vecterm(L5+3) = zadt(k)
         vecterm(L5+4) = frth(k)
         vecterm(L5+5) = frtv(k)
         vecterm(L5+6) = cadt(k)
      ENDDO
*
      write(27,*) 't-equation:'
      DO k=1,l
      write(27,919) xadt(k),yadt(k),zadt(k),frth(k),frtv(k),cadt(k)
      ENDDO
      write(29,*) 't-equation:'
      DO k=1,L
      write(29,909) xadt(k)+yadt(k)+zadt(k)+frth(k)+frtv(k)+cadt(k)
      ENDDO
*
      DO k=1,l
        h1 = 1.0/(dfzT(k)*dfzW(k))
        h2 = 1.0/(dfzT(k)*dfzW(k-1))
          xads(k) =-(u(i,j,k)*(s(i+1,j,k)+s(i,j,k))/(2*dx*cs(j))-
     &              u(i-1,j,k)*(s(i,j,k)+s(i-1,j,k))/(2*dx*cs(j)))
          yads(k) =-(v(i,j,k)*(s(i,j+1,k)+s(i,j,k))*csv(j)/(cs(j)*2*dy)-
     &              v(i,j-1,k)*(s(i,j,k)+s(i,j-1,k))*csv(j-1)/
     &              (cs(j)*2*dy))
*          zads(k) =-(w(i,j,k)*(s(i,j,k+1)+s(i,j,k))/(2*dz)-
*     &              w(i,j,k-1)*(s(i,j,k)+s(i,j,k-1))/(2*dz))
          zads(k) =-(w(i,j,k)*(s(i,j,k+1)+s(i,j,k))/(2*dz*dfzT(k))-
     &              w(i,j,k-1)*(s(i,j,k)+s(i,j,k-1))/(2*dz*dfzT(k)))
           frsh1(k) = (s(i+1,j,k)+s(i-1,j,k)-2*s(i,j,k))/
     &                (cs(j)*cs(j)*dx*dx)
           frsh2(k) = ((s(i,j+1,k)-s(i,j,k))*csv(j)-
     &                (s(i,j,k)-s(i,j-1,k))*csv(j-1))/
     &                (cs(j)*dy*dy)
          frsh(k) = ph*(frsh1(k)+frsh2(k))
*          frsv(k) = pv*(s(i,j,k+1)+s(i,j,k-1)-2*s(i,j,k))/(dz*dz)
          frsv(k) = pv*(h1*s(i,j,k+1)+h2*s(i,j,k-1)-
     &             (h1+h2)*s(i,j,k))/(dz*dz)
*          cads(k) = pvc*(hs(rho(i,j,k+1)-rho(i,j,k))*
*     &              (s(i,j,k+1)-s(i,j,k))-
*     &              hs(rho(i,j,k)-rho(i,j,k-1))*
*     &              (s(i,j,k)-s(i,j,k-1)))/(dz*dz)
          cads(k) = pvc*(h1*hs(rho(i,j,k+1)-rho(i,j,k))*s(i,j,k+1)
     &              + h2*hs(rho(i,j,k)-rho(i,j,k-1))*s(i,j,k-1)
     &              - (h1+h2)*( hs(rho(i,j,k+1)-rho(i,j,k))+
     &              hs(rho(i,j,k)-rho(i,j,k-1))*s(i,j,k)) )/(dz*dz)
      ENDDO
*
      DO k=1,l
         L6 = ntermu*l+ntermv*l+ntermw*l+ntermp*l+
     +        ntermt*l+nterms*(k-1)
         vecterm(L6+1) = xads(k)
         vecterm(L6+2) = yads(k)
         vecterm(L6+3) = zads(k)
         vecterm(L6+4) = frsh(k)
         vecterm(L6+5) = frsv(k)
         vecterm(L6+6) = cads(k)
      ENDDO
*
      write(27,*) 's-equation:'
      DO k=1,l
      write(27,919) xads(k),yads(k),zads(k),frsh(k),frsv(k),cads(k)
      ENDDO
      write(29,*) 's-equation:'
      DO k=1,L
      write(29,909) xads(k)+yads(k)+zads(k)+frsh(k)+frsv(k)+cads(k)
      ENDDO
 909  format(e18.10)
 919  format(6e18.10)
*
      END
******************************************************************
      SUBROUTINE idl2out(f,psimax,nout)
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      integer  nout
      real     f(0:m,0:l),psimax
*     LOCAL
      integer  j,k
      real     yc(0:m),zc(0:l)
*
      yc(0) = ymin
      DO j=1,m
         yc(j) = yv(j)
      ENDDO
      zc(0) = zmin
      DO k=1,l
         zc(k) = zw(k)
      ENDDO
      write(nout,999) m+1,l+1
      write(nout,998) psimax
      DO j=0,m
         DO k=0,l
            write(nout,998) yc(j),zc(k),f(j,k)
         ENDDO
      ENDDO
 999  format(2i8)
 998  format(3e18.10)
      END

**********************************************************
*** This function returns the index of the nearest cell
*** center for the specified latitude in degrees
      integer function latitude_index( lat )
      implicit none
      include 'usr.com'
      real lat, dif, difn
      integer i, j
      lat = lat * pi / 180.
      dif = pi
      do j = 1, m
         difn = abs(lat - y(j))
         if( difn .lt. dif) then
            dif = difn
            i = j
         endif
      enddo
      latitude_index = i
      end


