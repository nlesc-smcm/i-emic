********************************************************************
* Postprocessing written for analysis interdecadal mode
* Assumes that real part of eigenvector is in column jevr,
* imaginary part in column jevi
********************************************************************
      SUBROUTINE wts(un,w,t,s)
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      real    un(ndim)
      real    w(0:n+1,0:m+1,0:l),t(0:n+1,0:m+1,0:l+1)
      real    s(0:n+1,0:m+1,0:l+1)
*     LOCAL
      real    u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1)
      real    p(0:n+1,0:m+1,0:l+1)
*
* Compute fields and return only w, t and s
*
      call usol(un,u,v,w,p,t,s)
      END
********************************************************************
      SUBROUTINE bsper(upt,vpt,wpt,ppt,tpt,spt,jevr,jevi,ewre,ewim)
      implicit none
      include 'bag.com'
*     IMPORT/EXPORT
      integer jevr,jevi
*     LOCAL
      real    ewre, ewim
      integer i,j,k
      real    ur(0:n,0:m+1,0:l+1),vr(0:n+1,0:m,0:l+1)
      real    wr(0:n+1,0:m+1,0:l),pr(0:n+1,0:m+1,0:l+1)
      real    tr(0:n+1,0:m+1,0:l+1),sr(0:n+1,0:m+1,0:l+1)
      real    ui(0:n,0:m+1,0:l+1),vi(0:n+1,0:m,0:l+1)
      real    wi(0:n+1,0:m+1,0:l),pii(0:n+1,0:m+1,0:l+1)
      real    ti(0:n+1,0:m+1,0:l+1),si(0:n+1,0:m+1,0:l+1)
      real    test,e1,e2,sigia
      real    ure(ndim),uim(ndim)
      real    cst,cs(m),cst2
      real    pi
      integer t,tend
      parameter(tend=100)
      real    upt(0:n,0:m+1,0:l+1,0:tend)
      real    vpt(0:n+1,0:m,0:l+1,0:tend)
      real    wpt(0:n+1,0:m+1,0:l,0:tend)
      real    tpt(0:n+1,0:m+1,0:l+1,0:tend)
      real    spt(0:n+1,0:m+1,0:l+1,0:tend)
      real    ppt(0:n+1,0:m+1,0:l+1,0:tend)
*
* Compute fields for eigenvalue jevr
*
      IF ((jevr.EQ.0).OR.(jevi.EQ.0)) stop 'in bsper: 
     &          no of eigenvalue 0'
      ewre=sig(jevr,1)
      ewim=sig(jevr,2)
      sigia=abs(ewim)
	write(99,*) 'jevi=',jevi
	write(99,*) 'jevr=',jevr
	write(99,*) 'sigma=',ewre ,'+',ewim,'i'
	write(99,*) 'abs(ewim)=',sigia
        write(99,*) 'sigi1' ,ewim
      test=sig(jevr,2)+sig(jevi,2)
        write(99,*) 'sigi2' ,ewim
      ure=w(:,jevr)
      uim=w(:,jevi)
      write(99,*) 'sigi3' ,ewim
      call usol(ure,ur,vr,wr,pr,tr,sr)
      write(99,*) 'sigi3a' ,ewim
      call usol(uim,ui,vi,wi,pii,ti,si)
      write(99,*) 'sigi4' ,ewim
*
* Compute perturbations u, v, w, T and S
*
      pi    =  3.14159265358979323846
*      cst=2*pi/sigia*tend
      write(99,*) 'sigi5' ,ewim
      cst2=2*pi/tend
*	write(99,*) 'after cst'
*	write(99,*) cst
      write(99,*) 'sigi6' ,ewim
      DO i=1,n
       DO j=1,m
        DO k=1,l
         DO t=0,tend
           upt(i,j,k,t)=ur(i,j,k)*cos(cst2*t)-
     &                 ui(i,j,k)*sin(cst2*t)
           vpt(i,j,k,t)=vr(i,j,k)*cos(cst2*t)-
     &                 vi(i,j,k)*sin(cst2*t)
           wpt(i,j,k,t)=wr(i,j,k)*cos(cst2*t)-
     &                 wi(i,j,k)*sin(cst2*t)
           tpt(i,j,k,t)=tr(i,j,k)*cos(cst2*t)-
     &                 ti(i,j,k)*sin(cst2*t)
           spt(i,j,k,t)=sr(i,j,k)*cos(cst2*t)-
     &                 si(i,j,k)*sin(cst2*t)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
	write(99,*) 'end of bsper, ewim=', ewim
      END
********************************************************************
      SUBROUTINE structpert_orig(un,jevr,jevi)
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      real    un(ndim)
      integer jevr,jevi
*     LOCAL
      real   sigmar,sigmai
      integer i,j,k
      real    cst,cs(m)
      integer t,tend
      parameter(tend=100)
      real  wbs(0:n+1,0:m+1,0:l)
      real  tbs(0:n+1,0:m+1,0:l+1),sbs(0:n+1,0:m+1,0:l+1)
      real  rbs(0:n+1,0:m+1,0:l+1)
      real  upt(0:n,0:m+1,0:l+1,0:tend),vpt(0:n+1,0:m,0:l+1,0:tend)
      real  wpt(0:n+1,0:m+1,0:l,0:tend)
      real  ppt(0:n+1,0:m+1,0:l+1,0:tend)
      real  tpt(0:n+1,0:m+1,0:l+1,0:tend)
      real  spt(0:n+1,0:m+1,0:l+1,0:tend)
      real  rpt(0:n+1,0:m+1,0:l+1,0:tend)
      real    upsi(0:n,m,l,0:tend),vpsi(n,0:m,l,0:tend)
      real    us(0:n,l,0:tend),psiz(0:m,0:l,0:tend)
      real    vs(0:m,l,0:tend),psim(0:m,0:l,0:tend)
      real    dum,cs2
      real maxpm(0:tend),minpm(0:tend),maxpz(0:tend),minpz(0:tend)
      real mpm(0:tend),mpz(0:tend)
      real lambda
      real  sigia
*
      lambda = par(LAMB)
*
* Compute basic state w,T,S
*
      call wts(un,wbs,tbs,sbs)
*
* Compute perturbations u', v', w', p', T' and S'
*
      call bsper(upt,vpt,wpt,ppt,tpt,spt,jevr,jevi,sigmar,sigmai)
      write(99,*) 'beginning of structpert, sigmai=',sigmai
      sigia=abs(sigmai)
      write(99,*) 'jevr,jevi=', jevr,jevi
      write(99,*) 'structpert: sigmai=',sigmai
      write(99,*) 'structpert: abs(sigmai)=',sigia
      IF (sigia.EQ.0) stop 'abs(sigmai)=0 in structpert'
	write(99,*) 'hallo00'
      rpt=lambda*spt-tpt
	write(99,*) 'hallo0'
      rbs=lambda*sbs-tbs
*******************************************************************
* Temperature, salinity and density perturbations
*******************************************************************
*
* Write T' at surface at 8 snapshots 
*
* LtR      cst=2*pi/(abs(sigmai)*tend)
      cst=2*pi/(sigia*tend)
      DO i=1,n
         DO j=1,m
           write(70,998) tpt(i,j,l,0),tpt(i,j,l,tend/16),
     &                 tpt(i,j,l,tend/8),tpt(i,j,l,3*tend/16),
     &                 tpt(i,j,l,tend/4),tpt(i,j,l,5*tend/16),
     &                 tpt(i,j,l,3*tend/8),tpt(i,j,l,7*tend/16)
        ENDDO
      ENDDO

*
* Write S' at surface at 8 snapshots 
*
      DO i=1,n
         DO j=1,m
           write(71,998) spt(i,j,l,0),spt(i,j,l,tend/16),
     &                 spt(i,j,l,tend/8),spt(i,j,l,3*tend/16),
     &                 spt(i,j,l,tend/4),spt(i,j,l,5*tend/16),
     &                 spt(i,j,l,3*tend/8),spt(i,j,l,7*tend/16)
        ENDDO
      ENDDO

*
* Write rho' at surface at 8 snapshots 
*
      DO i=1,n
         DO j=1,m
           write(72,998) rpt(i,j,l,0),rpt(i,j,l,tend/16),
     &                 rpt(i,j,l,tend/8),rpt(i,j,l,3*tend/16),
     &                 rpt(i,j,l,tend/4),rpt(i,j,l,5*tend/16),
     &                 rpt(i,j,l,3*tend/8),rpt(i,j,l,7*tend/16)
        ENDDO
      ENDDO
	write(99,*) 'hallo1'
*******************************************************************
* Vertical velocity perturbations
*******************************************************************
*
* Write w' at surface at 8 snapshots 
*
      DO i=1,n
         DO j=1,m
           write(73,998) 0.5*(wpt(i,j,l,0)+wpt(i,j,l-1,0)),
     &              0.5*(wpt(i,j,l,tend/16)+wpt(i,j,l-1,tend/16)),
     &              0.5*(wpt(i,j,l,tend/8)+wpt(i,j,l-1,tend/8)),
     &              0.5*(wpt(i,j,l,3*tend/16)+wpt(i,j,l-1,3*tend/16)),
     &              0.5*(wpt(i,j,l,tend/4)+wpt(i,j,l-1,tend/4)),
     &              0.5*(wpt(i,j,l,5*tend/16)+wpt(i,j,l-1,5*tend/16)),
     &              0.5*(wpt(i,j,l,3*tend/8)+wpt(i,j,l-1,3*tend/8)),
     &              0.5*(wpt(i,j,l,7*tend/16)+wpt(i,j,l-1,7*tend/16))
        ENDDO
      ENDDO
*******************************************************************
* Meridional and zonal overturning perturbations
*******************************************************************
*
* Compute Psi_M'
*
      vs = 0.0
      vpsi(:,0,:,:)=0.0
      DO t=0,tend
       DO i=1,n
         DO j=1,m
           DO k=1,l
            vpsi(i,j,k,t)=vpt(i,j,k,t)
           ENDDO
         ENDDO
       ENDDO
      ENDDO
      DO t=0,tend
       DO j=0,m
         DO k=1,l
            dum = 0.0
            DO i=1,n
               dum=vpsi(i,j,k,t)*dx+dum
            ENDDO
            vs(j,k,t) = dum
         ENDDO
       ENDDO
      ENDDO
      psim(:,0,:) = 0.0
      psim(0,:,:) = 0.0
      DO t=0,tend
       DO j=1,m
         cs2 = cos(yv(j))
         DO k=1,l
            psim(j,k,t)=cs2*vs(j,k,t)*dz*dfzT(k)  + psim(j,k-1,t)
         ENDDO
       ENDDO
      ENDDO
	write(99,*) 'hallo2'
*
* Compute Psi_Z' 
*
      us = 0.0
      upsi(0,:,:,:)=0.0
      DO t=0,tend
       DO i=1,n
         DO j=1,m
           DO k=1,l
            upsi(i,j,k,t)=upt(i,j,k,t)
           ENDDO
         ENDDO
       ENDDO
      ENDDO
      DO t=0,tend
       DO i=0,n
         DO k=1,l
            dum = 0.0
            DO j=1,m
               cs2 = cos(y(j))
* Bug in computation Psi_Z' repaired!!!! 15/10/01               
*               dum=upsi(i,j,k,t)*dy/cs2+dum
               dum=upsi(i,j,k,t)*dy+dum
            ENDDO
            us(i,k,t) = dum
         ENDDO
       ENDDO
      ENDDO
      psiz(:,0,:) = 0.0
      psiz(0,:,:) = 0.0
      DO t=0,tend
       DO i=1,n
         DO k=1,l
            psiz(i,k,t)=us(i,k,t)*dz*dfzT(k)  + psiz(i,k-1,t)
         ENDDO
       ENDDO
      ENDDO
	write(99,*) 'hallo3'
*
* Write Psi_M' at 8 snapshots during oscillation
*
        DO j=1,m
          DO k=1,l
             write(74,998) psim(j,k,0),psim(j,k,tend/16),
     &                      psim(j,k,tend/8),psim(j,k,3*tend/16),
     &                      psim(j,k,tend/4),psim(j,k,5*tend/16),
     &                      psim(j,k,3*tend/8),psim(j,k,7*tend/16)
          ENDDO
        ENDDO
	write(99,*) 'hallo3'
*
* Write Psi_Z' at 8 snapshots during oscillation
*
        DO i=1,n
          DO k=1,l
             write(75,998) psiz(i,k,0),psiz(i,k,tend/16),
     &                      psiz(i,k,tend/8),psiz(i,k,3*tend/16),
     &                      psiz(i,k,tend/4),psiz(i,k,5*tend/16),
     &                      psiz(i,k,3*tend/8),psiz(i,k,7*tend/16)
          ENDDO
        ENDDO
*
 995  format(i4,6(X,e16.8e3))
 998  format(8(X,e16.8e3))
*
      END
********************************************************************
      SUBROUTINE tempdiff(un,jevr,jevi)
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      real    un(ndim)
      integer jevr,jevi
*     LOCAL
      integer i,j,k
      real    cst,cs(m)
      integer t,tend
      parameter(tend=100)
      real  sigmar,sigmai
      real  wbs(0:n+1,0:m+1,0:l)
      real  tbs(0:n+1,0:m+1,0:l+1),sbs(0:n+1,0:m+1,0:l+1)
      real  upt(0:n,0:m+1,0:l+1,0:tend),vpt(0:n+1,0:m,0:l+1,0:tend)
      real  wpt(0:n+1,0:m+1,0:l,0:tend)
      real  ppt(0:n+1,0:m+1,0:l+1,0:tend)
      real  tpt(0:n+1,0:m+1,0:l+1,0:tend)
      real  spt(0:n+1,0:m+1,0:l+1,0:tend)
      real    tavpt(n,m,0:tend)
      real    tnspt(1:n,0:tend), tewpt(1:m,0:tend)
      real    tnsxavpt(0:tend),tewyavpt(0:tend)
      real    tnspts(1:n,0:tend), tewpts(1:m,0:tend)
      real    tnsxavpts(0:tend),tewyavpts(0:tend)
      real    savpt(n,m,0:tend)
      real    snspt(1:n,0:tend), sewpt(1:m,0:tend)
      real    snsxavpt(0:tend),sewyavpt(0:tend)
      real    snspts(1:n,0:tend), sewpts(1:m,0:tend)
      real    snsxavpts(0:tend),sewyavpts(0:tend)
      real    rewyavpt(0:tend),rewyavpts(0:tend)
      real    rnsxavpt(0:tend),rnsxavpts(0:tend)
      real    lambda
      real    sigia
      
      lambda=par(LAMB)
*
* Compute basic state w,T,S
*
      call wts(un,wbs,tbs,sbs)
*
* Compute perturbations u', v', w', p', T' and S'
*
      call bsper(upt,vpt,wpt,ppt,tpt,spt,jevr,jevi,sigmar,sigmai)
      sigia=abs(sigmai)
      write(99,*) 'tempdiff: sigmai=',sigmai
      write(99,*) 'tempdiff: abs(sigmai)=',sigia
      IF (sigia.EQ.0) stop 'abs(sigmai)=0 in tempdiff'
*
* Integrate T' over upper layers
*
      tavpt=0.0
      DO t=0,tend
      DO i=1,n
         DO j=1,m
           DO k=6,l
             tavpt(i,j,t)=tpt(i,j,k,t)*dz*dfzT(k) + tavpt(i,j,t)
           ENDDO
         ENDDO
      ENDDO
      ENDDO
*
* Integrate S' over upper layers
*
      savpt=0.0
      DO t=0,tend
      DO i=1,n
         DO j=1,m
           DO k=6,l
             savpt(i,j,t)=spt(i,j,k,t)*dz*dfzT(k) + savpt(i,j,t)
           ENDDO
         ENDDO
      ENDDO
      ENDDO
****************************************************************
* East-west differences
****************************************************************
*
* East-west temperature difference, T_E - T_W
*
      DO t=0,tend
         DO j=1,m
             tewpt(j,t)=tavpt(n,j,t)-tavpt(1,j,t)
             tewpts(j,t)=tpt(n,j,l,t)-tpt(1,j,l,t)
         ENDDO
      ENDDO
*
* Average over y
*
      tewyavpt=0.0
      DO t=0,tend
         DO j=1,m
             tewyavpt(t)=tewpt(j,t)/m + tewyavpt(t)
             tewyavpts(t)=tewpts(j,t)/m + tewyavpts(t)
         ENDDO
      ENDDO
*
* East-west salinity difference, S_E - S_W
*
      DO t=0,tend
         DO j=1,m
             sewpt(j,t)=savpt(n,j,t)-savpt(1,j,t)
             sewpts(j,t)=spt(n,j,l,t)-spt(1,j,l,t)
         ENDDO
      ENDDO
*
* Average over y
*
      sewyavpt=0.0
      DO t=0,tend
         DO j=1,m
             sewyavpt(t)=sewpt(j,t)/m + sewyavpt(t)
             sewyavpts(t)=sewpts(j,t)/m + sewyavpts(t)
         ENDDO
      ENDDO
*
* East-west density difference, averaged over y
*
       rewyavpt=lambda*sewyavpt-tewyavpt
       rewyavpts=lambda*sewyavpts-tewyavpts
****************************************************************
* North-south differences
****************************************************************
*
* North-south temperature difference, T_n - T_s 
*
      DO t=0,tend
         DO i=1,n
             tnspt(i,t)=tavpt(i,m,t)-tavpt(i,1,t)
             tnspts(i,t)=tpt(i,m,l,t)-tpt(i,1,l,t)
         ENDDO
      ENDDO
*
* Average over x
*
      tnsxavpt=0.0
      DO t=0,tend
         DO i=1,n
             tnsxavpt(t)=tnspt(i,t)/n + tnsxavpt(t)
             tnsxavpts(t)=tnspts(i,t)/n + tnsxavpts(t)
         ENDDO
      ENDDO
*
      DO t=0,tend
         write(82,995) t,tewyavpt(t), tnsxavpt(t),
     &                 tewyavpts(t), tnsxavpts(t)
      ENDDO
*
*
* North-south salinity difference, S_n - S_s 
*
      DO t=0,tend
         DO i=1,n
             snspt(i,t)=savpt(i,m,t)-savpt(i,1,t)
             snspts(i,t)=spt(i,m,l,t)-spt(i,1,l,t)
         ENDDO
      ENDDO
*
* Average over x
*
      snsxavpt=0.0
      DO t=0,tend
         DO i=1,n
             snsxavpt(t)=snspt(i,t)/n + snsxavpt(t)
             snsxavpts(t)=snspts(i,t)/n + snsxavpts(t)
         ENDDO
      ENDDO
*
      DO t=0,tend
         write(83,995) t,sewyavpt(t), snsxavpt(t),
     &                 sewyavpts(t), snsxavpts(t)
      ENDDO
*
*
* North-south density difference, averaged over x
*
       rnsxavpt=lambda*snsxavpt-tnsxavpt
       rnsxavpts=lambda*snsxavpts-tnsxavpts
*
      DO t=0,tend
         write(84,995) t,rewyavpt(t), rnsxavpt(t),
     &                 rewyavpts(t), rnsxavpts(t)
      ENDDO
*
 995  format(i4,4(X,e16.8e3))
      END
********************************************************************
      SUBROUTINE wtpert(un,jevr,jevi)
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      real    un(ndim)
      integer jevr,jevi
*     LOCAL
      integer i,j,k
      real    cst,cs(m)
      real   sigmar,sigmai
      integer t,tend
      parameter(tend=100)
      real  wbs(0:n+1,0:m+1,0:l)
      real  tbs(0:n+1,0:m+1,0:l+1),sbs(0:n+1,0:m+1,0:l+1)
      real  upt(0:n,0:m+1,0:l+1,0:tend),vpt(0:n+1,0:m,0:l+1,0:tend)
      real  wpt(0:n+1,0:m+1,0:l,0:tend)
      real  ppt(0:n+1,0:m+1,0:l+1,0:tend)
      real  tpt(0:n+1,0:m+1,0:l+1,0:tend)
      real  spt(0:n+1,0:m+1,0:l+1,0:tend)
      real  rho(0:n+1,0:m+1,0:l+1)
      real  convt(n,m,l,0:tend),convs(n,m,l,0:tend)
      real  vdift(n,m,l,0:tend),vdifs(n,m,l,0:tend)
      real  iconv(3,tend), ivdif(3,tend)
      real  pv,pvc,ra,lambda,h1,h2,dum3,dum4
      real  rdz2i
      real  wbtp(n,m,l,0:tend),wptb(n,m,l,0:tend)
      real  wbsp(n,m,l,0:tend),wpsb(n,m,l,0:tend)
      real  intwbt(0:tend),intwtb(0:tend)
      real  intwbs(0:tend),intwsb(0:tend)
      real  intwbr(0:tend),intwrb(0:tend)
      real  dum1n(m,tend),dum2n(m,tend)
      real  dum1(0:m,tend),dum2(0:m,tend)
      real  iwbtp1(tend),iwptb1(tend),iwptb2(tend),iwbtp2(tend)
      real  iwbtp3(tend),iwptb3(tend),iwbtp4(tend),iwptb4(tend)
      real   sigia
*     EXTERNAL FUNCTIONS
      real      integral,hs
*      
      lambda = par(LAMB)
*
* Compute basic state w,T,S 
*
      call wts(un,wbs,tbs,sbs)
*
* Compute perturbations u', v', w', p', T' and S'
*
      call bsper(upt,vpt,wpt,ppt,tpt,spt,jevr,jevi,sigmar,sigmai)
      sigia=abs(sigmai)
      write(99,*) 'wtpert: sigmai=',sigmai
      write(99,*) 'wtpert: abs(sigmai)=',sigia
      IF (sigia.EQ.0) stop 'abs(sigmai)=0 in wtpert'
*
* Compute wT', w'T, wS' and w'S at T-points
*
      DO i=1,n
       DO j=1,m
        DO k=1,l
         DO t=0,tend
          wbtp(i,j,k,t)=0.5*(wbs(i,j,k)+wbs(i,j,k-1))*tpt(i,j,k,t)
          wptb(i,j,k,t)=0.5*(wpt(i,j,k,t)+wpt(i,j,k-1,t))*tbs(i,j,k)
          wbsp(i,j,k,t)=0.5*(wbs(i,j,k)+wbs(i,j,k-1))*spt(i,j,k,t)
          wpsb(i,j,k,t)=0.5*(wpt(i,j,k,t)+wpt(i,j,k-1,t))*sbs(i,j,k)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
*
* Compute volume integral over basin 
*
      cs = cos(y)
      intwbt = 0.0
      intwtb = 0.0
      intwbs = 0.0
      intwsb = 0.0
*
      DO t=1,tend
        DO i=1,n
          DO j=1,m
            DO k=1,l
           intwbt(t)=wbtp(i,j,k,t)*dx*cs(j)*dy*dz*dfzT(k)+intwbt(t)
           intwtb(t)=wptb(i,j,k,t)*dx*cs(j)*dy*dz*dfzT(k)+intwtb(t)
           intwbs(t)=wbsp(i,j,k,t)*dx*cs(j)*dy*dz*dfzT(k)+intwbs(t)
           intwsb(t)=wpsb(i,j,k,t)*dx*cs(j)*dy*dz*dfzT(k)+intwsb(t)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*
      DO t=1,tend
        intwbr(t)=lambda*intwbs(t)-intwbt(t)
        intwrb(t)=lambda*intwsb(t)-intwtb(t)
      ENDDO
*
* Write <w'T>, <wT'>, <w'S>, <wS'>, <w'rho> and <wrho'> 
* as a function of time
*
      cst=2*pi/(sigia*tend)
      write(76,997) tend
      write(76,998) sigmar,sigmai
      DO t=1,tend
        write(76,999) t,cst*t,intwbt(t),intwtb(t),intwbs(t),
     &                intwsb(t),intwbr(t),intwrb(t)   
      ENDDO

*
* Compute <z(P_VT'_z)_z> for comparison with <w'T> and <wT'>
*
      lambda = par(LAMB)
      ra     = par(RAYL)
      rho=lambda*sbs-tbs
*
      vdift = 0.0
      vdifs = 0.0
      ivdif = 0.0
*
      pv    = par(PE_V)
      rdz2i = (1.0/dz)**2
      DO t=0,tend
      DO k=1,l
        h1 = rdz2i/(dfzT(k)*dfzW(k))
        h2 = rdz2i/(dfzT(k)*dfzW(k-1))
        DO i=1,n
          DO j=1,m
            vdift(i,j,k,t) = (  h1*(tpt(i,j,k+1,t)-tpt(i,j,k,t))
     +                      + h2*(tpt(i,j,k-1,t)-tpt(i,j,k,t)) )*z(k)
            vdifs(i,j,k,t) = (  h1*(spt(i,j,k+1,t)-spt(i,j,k,t))
     +                      + h2*(spt(i,j,k-1,t)-spt(i,j,k,t)) )*z(k)
          ENDDO
        ENDDO
      ENDDO
      ENDDO
*
      DO t=1,tend
      ivdif(1,t) =   pv*integral(vdift(:,:,:,t))
      ivdif(2,t) =   pv*integral(vdifs(:,:,:,t))
      ivdif(3,t) = - ra*(ivdif(1,t)-lambda*ivdif(2,t))
      ENDDO
*
      convt = 0.0
      convs = 0.0
      iconv = 0.0
*
*      pvc   = par(P_VC)
*      rdz2i = (1.0/dz)**2
*      DO t=0,tend
**     DO k=2,l-1
*      DO k=1,l
*        h1 = 1.0/(dfzT(k)*dfzW(k))
*        h2 = 1.0/(dfzT(k)*dfzW(k-1))
*        DO i=1,n
*          DO j=1,m
*            dum3 = hs(rho(i,j,k+1)-rho(i,j,k))
*            dum4 = hs(rho(i,j,k)-rho(i,j,k-1))
*            convt(i,j,k,t)= ( h1*dum3*(tpt(i,j,k+1,t)-tpt(i,j,k,t))
*     +                    + h2*dum4*(tpt(i,j,k-1,t)-tpt(i,j,k,t)) )*z(k)
*            convs(i,j,k,t)= ( h1*dum3*(spt(i,j,k+1,t)-spt(i,j,k,t))
*     +                    + h2*dum4*(spt(i,j,k-1,t)-spt(i,j,k,t)) )*z(k)
*          ENDDO
*        ENDDO
*      ENDDO
*      ENDDO
*
*      DO t=1,tend
*      iconv(1,t) =   pvc*integral(convt(:,:,:,t))
*      iconv(2,t) =   pvc*integral(convs(:,:,:,t))
*      iconv(3,t) = - ra*(iconv(1,t)-lambda*iconv(2,t))
*      ENDDO
*
      write(77,997) tend
      write(77,998) sigmar,sigmai
      DO t=1,tend
        write(77,999) t,intwbt(t),intwtb(t),iconv(1,t)+ivdif(1,t),
     &                iconv(2,t)+ivdif(2,t),iconv(3,t)+ivdif(3,t)    
      ENDDO
*      
 997  format(i6)
 998  format(8(X,e16.8e3))
 999  format(i3,8(X,e16.8e3))
 996  format(2(i4),6(X,e16.8e3))
      END
********************************************************************
      real FUNCTION integral(fld)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      fld(n,m,l)
*     LOCAL
      integer   i,j,k
      real      cs(m)
*
* Take volume integral
*
      cs = cos(y)
      integral = 0.0
*
      DO i=1,n
        DO j=1,m
          DO k=1,l
            integral=fld(i,j,k)*dx*cs(j)*dy*dz*dfzT(k)+integral
          ENDDO
        ENDDO
      ENDDO
*
      END
********************************************************************
      SUBROUTINE structpert(un,jevr,jevi)
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      real    un(ndim)
      integer jevr,jevi
*     LOCAL
      real   sigmar,sigmai
      integer i,j,k
      real    cst,cs(m)
      integer t,tend
      parameter(tend=100)
      real  wbs(0:n+1,0:m+1,0:l)
      real  tbs(0:n+1,0:m+1,0:l+1),sbs(0:n+1,0:m+1,0:l+1)
      real  rbs(0:n+1,0:m+1,0:l+1)
      real  upt(0:n,0:m+1,0:l+1,0:tend),vpt(0:n+1,0:m,0:l+1,0:tend)
      real  wpt(0:n+1,0:m+1,0:l,0:tend)
      real  ppt(0:n+1,0:m+1,0:l+1,0:tend)
      real  tpt(0:n+1,0:m+1,0:l+1,0:tend)
      real  spt(0:n+1,0:m+1,0:l+1,0:tend)
      real  rpt(0:n+1,0:m+1,0:l+1,0:tend)
      real    upsi(0:n,m,l,0:tend),vpsi(n,0:m,l,0:tend)
      real    us(0:n,l,0:tend),psiz(0:m,0:l,0:tend)
      real    vs(0:m,l,0:tend),psim(0:m,0:l,0:tend)
      real    dum,cs2
      real maxpm(0:tend),minpm(0:tend),maxpz(0:tend),minpz(0:tend)
      real mpm(0:tend),mpz(0:tend)
      real lambda
      real  sigia
*
      lambda = par(LAMB)
*
* Compute basic state w,T,S
*
      call wts(un,wbs,tbs,sbs)
*
* Compute perturbations u', v', w', p', T' and S'
*
      call bsper(upt,vpt,wpt,ppt,tpt,spt,jevr,jevi,sigmar,sigmai)
      write(99,*) 'beginning of structpert, sigmai=',sigmai
      sigia=abs(sigmai)
      write(99,*) 'jevr,jevi=', jevr,jevi
      write(99,*) 'structpert: sigmai=',sigmai
      write(99,*) 'structpert: abs(sigmai)=',sigia
      IF (sigia.EQ.0) stop 'abs(sigmai)=0 in structpert'
	write(99,*) 'hallo00'
      rpt=lambda*spt-tpt
	write(99,*) 'hallo0'
      rbs=lambda*sbs-tbs
*******************************************************************
* Temperature, salinity and density perturbations
*******************************************************************
*
* Write T' at surface at 8 snapshots 
*
* LtR      cst=2*pi/(abs(sigmai)*tend)
      cst=2*pi/(sigia*tend)
      DO i=1,n
         DO j=1,m
           write(70,998) tpt(i,j,l,0),tpt(i,j,l,tend/16),
     &                 tpt(i,j,l,tend/8),tpt(i,j,l,3*tend/16),
     &                 tpt(i,j,l,tend/4),tpt(i,j,l,5*tend/16),
     &                 tpt(i,j,l,3*tend/8),tpt(i,j,l,7*tend/16)
        ENDDO
      ENDDO

*
* Write S' at surface at 8 snapshots 
*
      DO i=1,n
         DO j=1,m
           write(71,998) spt(i,j,l,0),spt(i,j,l,tend/16),
     &                 spt(i,j,l,tend/8),spt(i,j,l,3*tend/16),
     &                 spt(i,j,l,tend/4),spt(i,j,l,5*tend/16),
     &                 spt(i,j,l,3*tend/8),spt(i,j,l,7*tend/16)
        ENDDO
      ENDDO

*
* Write rho' at surface at 8 snapshots 
*
      DO i=1,n
         DO j=1,m
           write(72,998) rpt(i,j,l,0),rpt(i,j,l,tend/16),
     &                 rpt(i,j,l,tend/8),rpt(i,j,l,3*tend/16),
     &                 rpt(i,j,l,tend/4),rpt(i,j,l,5*tend/16),
     &                 rpt(i,j,l,3*tend/8),rpt(i,j,l,7*tend/16)
        ENDDO
      ENDDO
	write(99,*) 'hallo1'
*******************************************************************
* Vertical velocity perturbations
*******************************************************************
*
* Write w' at surface at 8 snapshots 
*
      DO i=1,n
         DO j=1,m
           write(73,998) 0.5*(wpt(i,j,l,0)+wpt(i,j,l-1,0)),
     &              0.5*(wpt(i,j,l,tend/16)+wpt(i,j,l-1,tend/16)),
     &              0.5*(wpt(i,j,l,tend/8)+wpt(i,j,l-1,tend/8)),
     &              0.5*(wpt(i,j,l,3*tend/16)+wpt(i,j,l-1,3*tend/16)),
     &              0.5*(wpt(i,j,l,tend/4)+wpt(i,j,l-1,tend/4)),
     &              0.5*(wpt(i,j,l,5*tend/16)+wpt(i,j,l-1,5*tend/16)),
     &              0.5*(wpt(i,j,l,3*tend/8)+wpt(i,j,l-1,3*tend/8)),
     &              0.5*(wpt(i,j,l,7*tend/16)+wpt(i,j,l-1,7*tend/16))
        ENDDO
      ENDDO
*******************************************************************
* Meridional and zonal overturning perturbations
*******************************************************************
*
* Compute Psi_M'
*
      vs = 0.0
      vpsi(:,0,:,:)=0.0
      DO t=0,tend
       DO i=1,n
         DO j=1,m
           DO k=1,l
            vpsi(i,j,k,t)=vpt(i,j,k,t)
           ENDDO
         ENDDO
       ENDDO
      ENDDO
      DO t=0,tend
       DO j=0,m
         DO k=1,l
            dum = 0.0
            DO i=1,n
               dum=vpsi(i,j,k,t)*dx+dum
            ENDDO
            vs(j,k,t) = dum
         ENDDO
       ENDDO
      ENDDO
      psim(:,0,:) = 0.0
      psim(0,:,:) = 0.0
      DO t=0,tend
       DO j=1,m
         cs2 = cos(yv(j))
         DO k=1,l
            psim(j,k,t)=cs2*vs(j,k,t)*dz*dfzT(k)  + psim(j,k-1,t)
         ENDDO
       ENDDO
      ENDDO
	write(99,*) 'hallo2'
*
* Compute Psi_Z' 
*
      us = 0.0
      upsi(0,:,:,:)=0.0
      DO t=0,tend
       DO i=1,n
         DO j=1,m
           DO k=1,l
            upsi(i,j,k,t)=upt(i,j,k,t)
           ENDDO
         ENDDO
       ENDDO
      ENDDO
      DO t=0,tend
       DO i=0,n
         DO k=1,l
            dum = 0.0
            DO j=1,m
               cs2 = cos(y(j))
* Bug in computation Psi_Z' repaired!!!! 15/10/01               
*               dum=upsi(i,j,k,t)*dy/cs2+dum
               dum=upsi(i,j,k,t)*dy+dum
            ENDDO
            us(i,k,t) = dum
         ENDDO
       ENDDO
      ENDDO
      psiz(:,0,:) = 0.0
      psiz(0,:,:) = 0.0
      DO t=0,tend
       DO i=1,n
         DO k=1,l
            psiz(i,k,t)=us(i,k,t)*dz*dfzT(k)  + psiz(i,k-1,t)
         ENDDO
       ENDDO
      ENDDO
	write(99,*) 'hallo3'
*
* Write Psi_M' at 8 snapshots during oscillation
*
        DO j=1,m
          DO k=1,l
             write(74,998) psim(j,k,0),psim(j,k,tend/16),
     &                      psim(j,k,tend/8),psim(j,k,3*tend/16),
     &                      psim(j,k,tend/4),psim(j,k,5*tend/16),
     &                      psim(j,k,3*tend/8),psim(j,k,7*tend/16)
          ENDDO
        ENDDO
	write(99,*) 'hallo3'
*
* Write Psi_Z' at 8 snapshots during oscillation
*
        DO i=1,n
          DO k=1,l
             write(75,998) psiz(i,k,0),psiz(i,k,tend/16),
     &                      psiz(i,k,tend/8),psiz(i,k,3*tend/16),
     &                      psiz(i,k,tend/4),psiz(i,k,5*tend/16),
     &                      psiz(i,k,3*tend/8),psiz(i,k,7*tend/16)
          ENDDO
        ENDDO
*
 995  format(i4,6(X,e16.8e3))
 998  format(8(X,e16.8e3))
*
      END
