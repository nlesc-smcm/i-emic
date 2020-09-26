******************************************************************
      SUBROUTINE potential(un,psimax,xl)
*
*     This routine reconstructs the advection-diffusion equations
*     and calculates the basin-integrated sources and sinks of
*     potential energy.
*
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      psimax,xl,un(ndim)
*     LOCAL
      real      ihdif(3),ivdif(3),idadv(3),ivadv(3),iconv(3)
      real      isurf(3),isum(3)
      real      u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +          w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +          t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      real      rho(0:n+1,0:m+1,0:l+1)
      integer   i,j,k
      real      hdift(n,m,l),hdifs(n,m,l),vdift(n,m,l),vdifs(n,m,l)
      real      vadvt(n,m,l),vadvs(n,m,l),dadvt(n,m,l),dadvs(n,m,l)
      real      convt(n,m,l),convs(n,m,l),surft(n,m,l),surfs(n,m,l)
      real      lambda,ra,ph,pv,pvc,h1,h2,dum1,dum2,gamma,etabi,bi
      real      tandyi(m),cosdx2i(m),csv(0:m),rdy2i(m),rdz2i
      real      cosdxi(m),rdyi(m),rdzi(l)
      real      advut(0:n,0:m+1,0:l+1),advus(0:n,0:m+1,0:l+1)
      real      advvt(0:n+1,0:m,0:l+1),advvs(0:n+1,0:m,0:l+1)
      real      advwt(0:n+1,0:m+1,0:l),advws(0:n+1,0:m+1,0:l)
      integer   ist,iend,jst,jend
      integer   west,south,center,north,east,bottom,top
      real      te(0:n+1,0:m+1,0:l+1),se(0:n+1,0:m+1,0:l+1)
      real      tw(0:n+1,0:m+1,0:l+1),sw(0:n+1,0:m+1,0:l+1)
      real      tn(0:n+1,0:m+1,0:l+1),sn(0:n+1,0:m+1,0:l+1)
      real      ts(0:n+1,0:m+1,0:l+1),sso(0:n+1,0:m+1,0:l+1)
      real      tto(0:n+1,0:m+1,0:l+1),st(0:n+1,0:m+1,0:l+1)
      real      tb(0:n+1,0:m+1,0:l+1),sb(0:n+1,0:m+1,0:l+1)
      common /integrate/ ist,iend,jst,jend
*     EXTERNAL FUNCTIONS
      real      int_vol,int_sur,hs
*
      ra     = par(RAYL)
      lambda = par(LAMB)
      ph     = par(PE_H)
      pv     = par(PE_V)
      pvc    = par(P_VC)
      gamma  = par(SALT)*(1-SRES+SRES*par(BIOT))
      etabi  = par(TEMP)*(1-TRES+TRES*par(BIOT))
      bi     = par(BIOT)
*
      call usol(un,u,v,w,p,t,s)
      rho    = lambda*s - t
*
* Initialize stencil points
*
      do k = 1, l
      do j = 1, m
      do i = 1, n
         west   = landm(i-1,j  ,k  )   ! 1
         south  = landm(i  ,j-1,k  )   ! 3
         center = landm(i  ,j  ,k  )   ! 4
         north  = landm(i  ,j+1,k  )   ! 5
         east   = landm(i+1,j  ,k  )   ! 7
         bottom = landm(i  ,j  ,k-1)   ! 8
         top    = landm(i  ,j  ,k+1)   ! 9
         if (center == OCEAN) then
            te(i,j,k) = t(i+1,j  ,k  )
            tw(i,j,k) = t(i-1,j  ,k  )
            tn(i,j,k) = t(i  ,j+1,k  )
            ts(i,j,k) = t(i  ,j-1,k  )
            tto(i,j,k) = t(i  ,j  ,k+1)
            tb(i,j,k) = t(i  ,j  ,k-1)
            se(i,j,k) = s(i+1,j  ,k  )
            sw(i,j,k) = s(i-1,j  ,k  )
            sn(i,j,k) = s(i  ,j+1,k  )
            sso(i,j,k) = s(i  ,j-1,k  )
            st(i,j,k) = s(i  ,j  ,k+1)
            sb(i,j,k) = s(i  ,j  ,k-1)
            if (west == LAND) then
               tw(i,j,k) = t(i,j,k)
               sw(i,j,k) = s(i,j,k)
            elseif (west == PERIO) then
               te(0,j,k) = t(1,j,k)
               se(0,j,k) = s(1,j,k)
            endif
            if (east == LAND) then
               te(i,j,k) = t(i,j,k)
               se(i,j,k) = s(i,j,k)
            endif
            if (south == LAND) then
               ts(i,j,k) = t(i,j,k)
               sso(i,j,k) = s(i,j,k)
            endif
            if (north == LAND) then
               tn(i,j,k) = t(i,j,k)
               sn(i,j,k) = s(i,j,k)
            endif
            if (bottom == LAND) then
               tb(i,j,k) = t(i,j,k)
               sb(i,j,k) = s(i,j,k)
            endif
*           if (top == LAND) then
            if (k.eq.l) then
               tto(i,j,k) = t(i,j,k)
               st(i,j,k) = s(i,j,k)
            endif
         endif
      enddo
      enddo
      enddo
*
*  A)  The terms in the advection-diffusion balance are reconstructed
*
*  1)  Horizontal diffusion
*
      hdift  = 0.0
      hdifs  = 0.0
      ihdif  = 0.0
*
      csv      = cos(yv)
      rdy2i    = 1.0/(cos(y)*dy**2)
      tandyi   = tan(y)/(2*dy)
      cosdx2i  = (1.0/(cos(y)*dx))**2
      DO k=1,l
        DO j=1,m
          DO i=1,n
            hdift(i,j,k)= cosdx2i(j)*( te(i,j,k)+tw(i,j,k)-2*t(i,j,k))
     +                  +   rdy2i(j)*((tn(i,j,k)-t(i,j,k))*csv(j)
     +                               +(ts(i,j,k)-t(i,j,k))*csv(j-1))
*    +                  -  tandyi(j)*( tn(i,j,k)-ts(i,j,k)) )
            hdifs(i,j,k)= cosdx2i(j)*( se(i,j,k)+sw(i,j,k)-2*s(i,j,k))
     +                  +   rdy2i(j)*((sn(i,j,k)-s(i,j,k))*csv(j)
     +                               +(sso(i,j,k)-s(i,j,k))*csv(j-1))
*    +                  -  tandyi(j)*( sn(i,j,k)-sso(i,j,k)) )
          ENDDO
        ENDDO
      ENDDO
*
      ihdif(1)  =   ph*int_vol(hdift,ist,iend,jst,jend,0)
      ihdif(2)  =   ph*int_vol(hdifs,ist,iend,jst,jend,0)
      ihdif(3)  = - ra*(ihdif(1)-lambda*ihdif(2))  ! < div_h . grad_h rho >
*
*
*  2) Vertical diffusion
*
      vdift  = 0.0
      vdifs  = 0.0
      ivdif  = 0.0
*
      rdz2i = (1.0/dz)**2
      DO k=1,l
        h1 = rdz2i/(dfzT(k)*dfzW(k))
        h2 = rdz2i/(dfzT(k)*dfzW(k-1))
        DO j=1,m
          DO i=1,n
            vdift(i,j,k)  = (  h1*(tto(i,j,k)-t(i,j,k))
     +                      + h2*(tb(i,j,k)-t(i,j,k)) )
            vdifs(i,j,k)  = (  h1*(st(i,j,k)-s(i,j,k))
     +                      + h2*(sb(i,j,k)-s(i,j,k)) )
          ENDDO
        ENDDO
      ENDDO
*
      ivdif(1)  =   pv*int_vol(vdift,ist,iend,jst,jend,0)
      ivdif(2)  =   pv*int_vol(vdifs,ist,iend,jst,jend,0)
      ivdif(3)  = - ra*(ivdif(1)-lambda*ivdif(2)) ! < d^2/dz^2 rho >
*
*  3) Convection (Enhanced vertical diffusion: old and new procedure)
*
      convt = 0.0
      convs = 0.0
      iconv = 0.0
*
      rdz2i = (1.0/dz)**2
      DO k=1,l
        h1 = rdz2i/(dfzT(k)*dfzW(k))
        h2 = rdz2i/(dfzT(k)*dfzW(k-1))
        DO j=1,m
          DO i=1,n
            dum1 = pvc * hs(rho(i,j,k+1)-rho(i,j,k))
     +           + pv  * (pv_adj(i,j,k  )-1.)
            dum2 = pvc * hs(rho(i,j,k)-rho(i,j,k-1))
     +           + pv  * (pv_adj(i,j,k-1)-1.)
            convt(i,j,k)= ( h1*dum1*(tto(i,j,k)-t(i,j,k))
     +                    + h2*dum2*(tb(i,j,k)-t(i,j,k)) )
            convs(i,j,k)= ( h1*dum1*(st(i,j,k)-s(i,j,k))
     +                    + h2*dum2*(sb(i,j,k)-s(i,j,k)) )
          ENDDO
        ENDDO
      ENDDO
*
      iconv(1) =   int_vol(convt,ist,iend,jst,jend,0)
      iconv(2) =   int_vol(convs,ist,iend,jst,jend,0)
      iconv(3) = - ra*(iconv(1)-lambda*iconv(2))
*
*
*  4) Advection
*
* Calculate advection on u,v en w-points first
*
      advut = 0.0
      advus = 0.0
      advvt = 0.0
      advvs = 0.0
      advwt = 0.0
      advws = 0.0
*
      DO k=1,l
        DO j=1,m
          DO i=0,n
            advut(i,j,k) = 0.5*u(i,j,k)*(te(i,j,k)+t(i,j,k))
            advus(i,j,k) = 0.5*u(i,j,k)*(se(i,j,k)+s(i,j,k))
            advvt(i,j,k) = 0.5*v(i,j,k)*(tn(i,j,k)+t(i,j,k))*csv(j)
            advvs(i,j,k) = 0.5*v(i,j,k)*(sn(i,j,k)+s(i,j,k))*csv(j)
            advwt(i,j,k) = 0.5*w(i,j,k)*(tto(i,j,k)+t(i,j,k))
            advws(i,j,k) = 0.5*w(i,j,k)*(st(i,j,k)+s(i,j,k))
          ENDDO
        ENDDO
      ENDDO
*
*  4a) Divergence of vertical advection
*
      rdzi    = 1.0/(dfzT*dz)
      DO k=1,l
        DO j=1,m
          DO i=1,n
            vadvt(i,j,k)= - rdzi(k)*(advwt(i,j,k)-advwt(i,j,k-1)) 
            vadvs(i,j,k)= - rdzi(k)*(advws(i,j,k)-advws(i,j,k-1))
          ENDDO
        ENDDO
      ENDDO
*
      ivadv(1) = int_vol(vadvt,ist,iend,jst,jend,0)
      ivadv(2) = int_vol(vadvs,ist,iend,jst,jend,0)
      ivadv(3) = - ra*(ivadv(1)-lambda*ivadv(2))   ! - < d/dz w rho  >
*        
*  4b) Divergence of horizontal advection
*
      dadvt = 0.0
      dadvs = 0.0
      idadv = 0.0
*
      rdyi    = 1.0/(cos(y)*dy)
      cosdxi  = 1.0/(cos(y)*dx)
      rdzi    = 1.0/(dfzT*dz)
      DO k=1,l
        DO j=1,m
          DO i=1,n
            dadvt(i,j,k)= - cosdxi(j)*(advut(i,j,k)-advut(i-1,j,k))
     +                    - rdyi(j)  *(advvt(i,j,k)-advvt(i,j-1,k)) 
            dadvs(i,j,k)= - cosdxi(j)*(advus(i,j,k)-advus(i-1,j,k))
     +                    - rdyi(j)  *(advvs(i,j,k)-advvs(i,j-1,k)) 
          ENDDO
        ENDDO
      ENDDO
*
      idadv(1) = int_vol(dadvt,ist,iend,jst,jend,0)
      idadv(2) = int_vol(dadvs,ist,iend,jst,jend,0)
      idadv(3) = - ra*(idadv(1)-lambda*idadv(2))  ! - < nabla_h u rho >
*
*  5) Surface forcing
*
      surft = 0.0
      surfs = 0.0
      isurf = 0.0
*
      k=l
      DO j=1,m
        DO i=1,n
          surft(i,j,k)  = - bi*TRES*t(i,j,l) + etabi*ft(i,j)
          surfs(i,j,k)  = - bi*SRES*s(i,j,l) + gamma*fs(i,j)
        ENDDO
      ENDDO
*
      isurf(1) =   int_vol(surft,ist,iend,jst,jend,0)
      isurf(2) =   int_vol(surfs,ist,iend,jst,jend,0)
      isurf(3) = - ra*(isurf(1)-lambda*isurf(2))  ! < F_rho g(z) >
*
*   Write integrals or pointwise balances to file:
*
*      call write_ts_int(ihdif,ivdif,iconv,ivadv,idadv,isurf,xl,psimax,83)
*
*      call write_ts(hdift,vdift,convt,vadvt,dadvt,surft,'T')
*      call write_ts(hdifs,vdifs,convs,vadvs,dadvs,surfs,'S')
*
*
*  B) Multiply terms of advection-diffusion equations with z to obtain
*     potential energy production/dissipation terms
*
      DO k=1,l
        DO j=1,m
          DO i=1,n
            hdift(i,j,k) = z(k)*hdift(i,j,k)
            hdifs(i,j,k) = z(k)*hdifs(i,j,k)
            vdift(i,j,k) = z(k)*vdift(i,j,k)
            vdifs(i,j,k) = z(k)*vdifs(i,j,k)
            vadvt(i,j,k) = z(k)*vadvt(i,j,k)
            vadvs(i,j,k) = z(k)*vadvs(i,j,k)
            dadvt(i,j,k) = z(k)*dadvt(i,j,k)
            dadvs(i,j,k) = z(k)*dadvs(i,j,k)
            convs(i,j,k) = z(k)*convs(i,j,k)
            convt(i,j,k) = z(k)*convt(i,j,k)
            surft(i,j,k) = z(l)*surft(i,j,k)
            surfs(i,j,k) = z(l)*surfs(i,j,k)
          ENDDO
        ENDDO
      ENDDO
*
*  4c) Vertical advection: Explicit calculation of buoyancy production
*      <w rho>. In integrated sense it should equal -<z d/dz (w rho)>.
*      Not correct for local balances, since: 
*       -z d/dz (w rho) = - d/dz (z w rho) + w rho 
*
      if (.false.) then
      DO k=1,l
        DO j=1,m
          DO i=1,n
            h1 = dfzW(k)/(2*dfzT(k))
            h2 = dfzW(k-1)/(2*dfzT(k))
            vadvt(i,j,k)= h2*advwt(i,j,k-1)+h1*advwt(i,j,k)
            vadvs(i,j,k)= h2*advws(i,j,k-1)+h1*advws(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      endif
*
      ihdif(1) =   ph*int_vol(hdift,ist,iend,jst,jend,0)
      ihdif(2) =   ph*int_vol(hdifs,ist,iend,jst,jend,0)
      ihdif(3) = - ra*(ihdif(1)-lambda*ihdif(2))      ! < z nabla_h^2 rho >
*
      ivdif(1) =   pv*int_vol(vdift,ist,iend,jst,jend,0)
      ivdif(2) =   pv*int_vol(vdifs,ist,iend,jst,jend,0)
      ivdif(3) = - ra*(ivdif(1)-lambda*ivdif(2))      ! < z d^2/dz^2 rho >
*
      iconv(1) =      int_vol(convt,ist,iend,jst,jend,0)
      iconv(2) =      int_vol(convs,ist,iend,jst,jend,0)
      iconv(3) = - ra*(iconv(1)-lambda*iconv(2))
*
      ivadv(1) =      int_vol(vadvt,ist,iend,jst,jend,0)
      ivadv(2) =      int_vol(vadvs,ist,iend,jst,jend,0)
      ivadv(3) = - ra*(ivadv(1)-lambda*ivadv(2))      ! -<z d/dz (w rho)  >
*
      idadv(1) =      int_vol(dadvt,ist,iend,jst,jend,0)
      idadv(2) =      int_vol(dadvs,ist,iend,jst,jend,0)
      idadv(3) = - ra*(idadv(1)-lambda*idadv(2))      ! -< z nabla_h (u rho) >
*
      isurf(1) =      int_vol(surft,ist,iend,jst,jend,0)
      isurf(2) =      int_vol(surfs,ist,iend,jst,jend,0)
      isurf(3) = - ra*(isurf(1)-lambda*isurf(2))      ! < z F_rho g(z) >
*
      call write_ts_int(ihdif,ivdif,iconv,ivadv,idadv,isurf,xl,
     +                    psimax,81)
*
*      call write_ts(hdift,vdift,convt,vadvt,dadvt,surft,'T')
*      call write_ts(hdifs,vdifs,convs,vadvs,dadvs,surfs,'S')
*
      END
******************************************************************
      SUBROUTINE write_ts_int(ihdif,ivdif,iconv,ivadv,
     +                                 idadv,isurf,xl,psimax,lu)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      ihdif(3),ivdif(3),iconv(3),ivadv(3),idadv(3),isurf(3)
      real      xl,psimax
      integer   lu
*     LOCAL
      real      isum(3)
      logical   pc
*
      isum = ihdif+ivdif+ivadv+idadv+iconv+isurf
* 
* pc = true : plot format (single line)
* pc = false: colums
*
      pc = .true.
      if (pc) then
        write(lu,981)xl,ihdif,ivdif,iconv,ivadv,idadv,
     +                  isurf,isum,psimax
      else
        write(lu,981)xl
        write(lu,982)'hdif ',ihdif
        write(lu,982)'vdif ',ivdif
        write(lu,982)'conv ',iconv
        write(lu,982)'vadv ',ivadv
        write(lu,982)'hadv ',idadv
        write(lu,982)'surf ',isurf
        write(lu,*)'-----------------------------------------+'
        write(lu,982)'sum  ',isum
        write(lu,*)             
      endif
*
 981  format(23e12.4)
 982  format(a5,3e12.4)
*
      return
      end
******************************************************************
      SUBROUTINE write_ts(hdif,vdif,conv,vadv,dadv,surf,TS)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      hdif(n,m,l),vdif(n,m,l),conv(n,m,l),vadv(n,m,l),
     +          dadv(n,m,l),surf(n,m,l)
      character*1 TS
*     LOCAL
      real      sum(n,m,l)
      real      ph,pv,pvc
      real      maxsum
      integer   i,j,k,imax,jmax,kmax
*
      ph     = par(PE_H)
      pv     = par(PE_V)
      pvc    = par(P_VC)
*
      sum = ph*hdif+pv*vdif+pvc*conv+vadv+dadv+surf
*
      maxsum = 0.0
      imax = 0
      jmax = 0
      kmax = 0
*
      do k=1,l
        do j=1,m
          do i=1,n
          if (landm(i,j,k).eq.OCEAN) then
            write(96,996) i,j,k,sum(i,j,k),ph*hdif(i,j,k),
     +                      pv*vdif(i,j,k),pvc*conv(i,j,k),
     +                      vadv(i,j,k),dadv(i,j,k),
     +                      surf(i,j,k)
            if ( abs(sum(i,j,k)).gt.maxsum ) then 
              maxsum = abs(sum(i,j,k))
              imax = i
              jmax = j
              kmax = k
            endif
          endif
          enddo
        enddo
        write(96,*) 
      enddo
      write(96,997) TS,'-equation largest residual : ',
     +                   maxsum, ' at ',imax,jmax,kmax
      write(96,*) 
*
 996  format(3i3,7e12.4)
 997  format(a1,a,1e12.4,a,3i3)
*
      return
      end
******************************************************************
      SUBROUTINE kinetic(un,psimax,xl)
*
* This routine reconstructs the momentum equations and computes 
* the source/sink terms of kinetic energy. 
* 
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      psimax,xl,un(ndim)
*     LOCAL
      real      ibuf,idpa,iwsw,ifrv,ifrh,icor,inln,iugp,isum
      real      u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +          w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +          t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      real      uun(0:n,0:m+1,0:l+1),uus(0:n,0:m+1,0:l+1)
      real      uub(0:n,0:m+1,0:l+1),uut(0:n,0:m+1,0:l+1)
      real      vve(0:n+1,0:m,0:l+1),vvw(0:n+1,0:m,0:l+1)
      real      vvb(0:n+1,0:m,0:l+1),vvt(0:n+1,0:m,0:l+1)
      real      vune(0:n,0:m+1,0:l+1),vuse(0:n+1,0:m,0:l+1)
      real      vunw(0:n,0:m+1,0:l+1),vusw(0:n+1,0:m,0:l+1)
      real      uvne(0:n,0:m+1,0:l+1),uvse(0:n,0:m+1,0:l+1)
      real      uvnw(0:n,0:m+1,0:l+1),uvsw(0:n,0:m+1,0:l+1)
      real      rho(0:n+1,0:m+1,0:l+1)
      integer   i,j,k
      real      lambda,ra,sigma,eh,ev,rota,epsr
      real      cs,c1,c2,c3,ci,h1,h2,p5
      real      rdy2i,rdz2i,fact,usum1,usum2,vsum
      real      rdyi,cosdxi(m),cos2dxi(m),cos2dyi(m),cos2dyiv(0:m)
      real      cos2dxiv(0:m)
      real      cos2i(m),tandyi(m),cosdx2i(m),sincos2i(m),tdzi(l)
      real      cos2iv(0:m),tandyiv(0:m),cosdx2iv(0:m),sincos2iv(0:m)
      real      buf(n,m,l),dpa(n,m,l),wsw(n,m,l),frv(n,m,l),frh(n,m,l)
      real      cor(n,m,l),nln(n,m,l),ugp(n,m,l)
      real      hvu(n,m,l,np),hvv(n,m,l,np)
      real      frhu(0:n,0:m+1,0:l+1),frhv(0:n+1,0:m,0:l+1)
      real      frvu(0:n,0:m+1,0:l+1),frvv(0:n+1,0:m,0:l+1)
      real      coru(0:n,0:m+1,0:l+1),corv(0:n+1,0:m,0:l+1)
      real      ugpu(0:n,0:m+1,0:l+1),ugpv(0:n+1,0:m,0:l+1)
      real      wswu(0:n,0:m+1,0:l+1),wswv(0:n+1,0:m,0:l+1)
      real      nlnu(0:n,0:m+1,0:l+1),nlnv(0:n+1,0:m,0:l+1)
      real      uux(0:n,0:m+1,0:l+1),uvx(0:n+1,0:m,0:l+1)
      real      vuy(0:n,0:m+1,0:l+1),vvy(0:n+1,0:m,0:l+1)
      real      wuz(0:n,0:m+1,0:l+1),wvz(0:n+1,0:m,0:l+1)
      real      uvt(0:n,0:m+1,0:l+1),u2t(0:n+1,0:m,0:l+1)
      integer east, west, north, south, nowe, soea, top, bottom, center
*     COMMON
      integer ist,iend,jst,jend
      common /integrate/ ist,iend,jst,jend
*     EXTERNAL FUNCTIONS
      real      int_vol,int_sur
*
      ra     = par(RAYL)
      epsr   = par(ROSB)
      lambda = par(LAMB)
      ev     = par(EK_V)
      eh     = par(EK_H)
      sigma  = par(WIND)*par(AL_T)/(dz*dfzT(l))
      rota   = 0.0 !par(ROTN) ! ATvS-Mix
*
      call usol(un,u,v,w,p,t,s)
      rho    = lambda*s - t
*
* Initialize stencil points and uv and vu
*
      vune = 0.0
      vuse = 0.0
      vunw = 0.0
      vusw = 0.0
      uvne = 0.0
      uvse = 0.0
      uvnw = 0.0
      uvsw = 0.0
      uun  = 0.0
      uus  = 0.0
      uub  = 0.0
      uut  = 0.0
      vve  = 0.0
      vvw  = 0.0
      vvb  = 0.0
      vvt  = 0.0
      do k = 1, l
      do j = 1, m
      do i = 1, n
         west   = landm(i-1,j  ,k  )   ! 1
         nowe   = landm(i-1,j+1,k  )   ! 2
         south  = landm(i  ,j-1,k  )   ! 3
         center = landm(i  ,j  ,k  )   ! 4
         north  = landm(i  ,j+1,k  )   ! 5
         soea   = landm(i+1,j-1,k  )   ! 6
         east   = landm(i+1,j  ,k  )   ! 7
         bottom = landm(i  ,j  ,k-1)   ! 8
         top    = landm(i  ,j  ,k+1)   ! 9
         if (center == OCEAN) then
            vune(i,j,k) = v(i+1,j  ,k)
            vuse(i,j,k) = v(i+1,j-1,k)
            vunw(i,j,k) = v(i  ,j  ,k)
            vusw(i,j,k) = v(i  ,j-1,k)
            uvne(i,j,k) = u(i  ,j+1,k)
            uvse(i,j,k) = u(i  ,j  ,k)
            uvnw(i,j,k) = u(i-1,j+1,k)
            uvsw(i,j,k) = u(i-1,j  ,k)
            uun(i,j,k)  = u(i  ,j+1,k)
            uus(i,j,k)  = u(i  ,j-1,k)
            uub(i,j,k)  = u(i  ,j,k-1)
            uut(i,j,k)  = u(i  ,j,k+1)
            vve(i,j,k)  = v(i+1,j  ,k)
            vvw(i,j,k)  = v(i-1,j  ,k)
            vvb(i,j,k)  = v(i  ,j,k-1)
            vvt(i,j,k)  = v(i  ,j,k+1)
            if (west == LAND) then
               vvw(i,j,k)  = SLIP * v(i,j,k)
               uvnw(i,j,k) = 0.0
               uvsw(i,j,k) = 0.0
            endif
            if (east == LAND) then
               vve(i,j,k)  = SLIP * v(i,j,k)
               vune(i,j,k) = 0.0
               vuse(i,j,k) = 0.0
               vunw(i,j,k) = 0.0
               vusw(i,j,k) = 0.0
               uvne(i,j,k) = 0.0
               uvse(i,j,k) = 0.0
            endif
            if (south == LAND) then
               uus(i,j,k)   = SLIP * u(i,j,k)
               vuse(i,j,k) = 0.0
               vusw(i,j,k) = 0.0
            endif
            if (north == LAND) then
               uun(i,j,k)   = SLIP * u(i,j,k)
               uvne(i,j,k) = 0.0
               uvse(i,j,k) = 0.0
               uvnw(i,j,k) = 0.0
               uvsw(i,j,k) = 0.0
               vune(i,j,k) = 0.0
               vunw(i,j,k) = 0.0
            endif
            if (bottom == LAND) then
               uub(i,j,k) = SLIP * u(i,j,k)
               vvb(i,j,k) = SLIP * v(i,j,k)
            endif
*           if (top == LAND) then
            if (k.eq.l) then
               uut(i,j,k) = u(i,j,k)
               vvt(i,j,k) = v(i,j,k)
            endif
         endif
      enddo
      enddo
      enddo
*
*  1)  Pressure Work P
*
*  1a) - u_h . grad_h p  ( = - div . (u p) + w dp/dz )
*
      ugp   = 0.0
      iugp  = 0
      ugpu  = 0.0
      ugpv  = 0.0
*
      rdyi   = 1.0/dy
      cosdxi = 1.0/(cos(y)*dx)
      DO k=1,l
        DO j=1,m
          DO i=1,n
            ugpu(i,j,k) = - cosdxi(j)*(p(i+1,j,k)-p(i,j,k))
            ugpv(i,j,k) = - rdyi*(p(i,j+1,k)-p(i,j,k))
          ENDDO
          ugpu(0,j,k) = ugpu(n,j,k)
        ENDDO
      ENDDO
*
      call u_in_fu(u,v,ugpu,ugpv,ugp)
*
      iugp = int_vol(ugp,ist,iend,jst,jend,0)
*
*  1b)  buoyancy forcing wB:  w dp/dz = - ra w rho
*
      buf  = 0.0
      ibuf = 0.0
*
      DO k=1,l
        DO j=1,m
          DO i=1,n
            h1 = dfzW(k)/(2*dfzT(k))
            h2 = dfzW(k-1)/(2*dfzT(k))
            buf(i,j,k)=-h2*w(i,j,k-1)*ra*0.5*(rho(i,j,k)+rho(i,j,k-1))
     +                 -h1*w(i,j,k)  *ra*0.5*(rho(i,j,k)+rho(i,j,k+1))
          ENDDO
        ENDDO
      ENDDO
*
      ibuf = int_vol(buf,ist,iend,jst,jend,0)
*
*  1c) divergence of pressure-work: - div . (u p)
*
      dpa  = 0.0
      idpa = 0.0
*
      DO k=1,l
        DO j=1,m
          cs = cos(y(j))
          c1 = cos(yv(j))
          c2 = cos(yv(j-1))
          DO i=1,n
            dpa(i,j,k)=-((p(i+1,j,k)+p(i,j,k))*u(i,j,k) -
     +                   (p(i,j,k)+p(i-1,j,k))*u(i-1,j,k))/(2*dx*cs)
     +                 -((p(i,j+1,k)+p(i,j,k))*v(i,j,k)*c1 - 
     +                   (p(i,j,k)+p(i,j-1,k))*v(i,j-1,k)*c2)/(2*dy*cs) 
     +                 -((p(i,j,k+1)+p(i,j,k))*w(i,j,k)-
     +                  (p(i,j,k)+p(i,j,k-1))*w(i,j,k-1))/(2*dz*dfzT(k))
          ENDDO
        ENDDO
      ENDDO
*
      idpa = int_vol(dpa,ist,iend,jst,jend,0)
*
*
*  2)  Wind Stress Work T
*
      wsw  = 0.0
      iwsw = 0.0
      wswu  = 0.0
      wswv  = 0.0
*
      DO j=1,m
        DO i=1,n-1
          wswu(i,j,l) =  0.5*(tx(i,j)+tx(i+1,j))
        ENDDO
      ENDDO
      DO j=1,m-1
        DO i=1,n
          wswv(i,j,l) =  0.5*(ty(i,j)+ty(i,j+1))
        ENDDO
      ENDDO
*
      call u_in_fu(u,v,wswu,wswv,wsw)
*
      iwsw = sigma * int_vol(wsw,ist,iend,jst,jend,0)
*
*
*  3)  Dissipation by Vertical Friction D_V
*
      frv   = 0.0
      ifrv  = 0.0
      frvu  = 0.0
      frvv  = 0.0
*
      rdz2i = (1.0/dz)**2
      DO k=1,l
        h1 = rdz2i/(dfzT(k)*dfzW(k))
        h2 = rdz2i/(dfzT(k)*dfzW(k-1))
        DO j=1,m
          DO i=1,n
            frvu(i,j,k) = h2*uub(i,j,k)+h1*uut(i,j,k)
     +                    -(h2+h1)*u(i,j,k)
            frvv(i,j,k) = h2*vvb(i,j,k)+h1*vvt(i,j,k)
     +                     -(h2+h1)*v(i,j,k)
          ENDDO
          frvu(0,j,k) = frvu(n,j,k)
        ENDDO
      ENDDO
*
      call u_in_fu(u,v,frvu,frvv,frv)
*
      ifrv = ev * int_vol(frv,ist,iend,jst,jend,0)
*
*
*  4)  Dissipaton by Horizontal Friction D_H
*
      frh   = 0.0
      ifrh  = 0.0
      frhu  = 0.0
      frhv  = 0.0
*
*     call hovisc(1,hvu)
*     call hovisc(2,hvv)
*
      rdy2i    = (1.0/dy)**2
      cos2i    = 1.0/cos(y)**2
      tandyi   = tan(y)/(2*dy)
      cosdx2i  = (1.0/(cos(y)*dx))**2
      sincos2i = sin(y)/(dx*cos(y)**2)
      DO k=1,l
        DO j=1,m
*
          DO i=1,n-1
            frhu(i,j,k)=cosdx2i(j)*(u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))
     +                + rdy2i     *(uun(i,j,k)-2*u(i,j,k)+uus(i,j,k))
     +                -  tandyi(j)*(uun(i,j,k)-uus(i,j,k))
     +                -   cos2i(j)* u(i,j,k)
     +                -sincos2i(j)*(vune(i,j,k)+vuse(i,j,k)
     +                             -vunw(i,j,k)-vusw(i,j,k))
*            frhu(i,j,k) = frhu(i,j,k)*hvu(i,j,k,1)
          ENDDO
          if (periodic) then
            i=n
            frhu(n,j,k)=cosdx2i(j)*(u(1,j,k)-2*u(n,j,k)+u(n-1,j,k))
     +              + rdy2i     *(uun(n,j,k)-2*u(n,j,k)+uus(n,j,k))
     +              -  tandyi(j)*(uun(n,j,k)-uus(n,j,k))
     +              -   cos2i(j)* u(n,j,k)
     +              -sincos2i(j)*(vune(n,j,k)+vuse(n,j,k)
     +                           -vunw(n,j,k)-vusw(n,j,k))
            frhu(0,j,k)=frhu(n,j,k)
*            frhu(n,j,k) = frhu(n,j,k)*hvu(n,j,k,1)
*            frhu(0,j,k) = frhu(0,j,k)*hvu(0,j,k,1)
          endif
*
        ENDDO
      ENDDO
*
      cos2iv    = 1.0/cos(yv)**2
      tandyiv   = tan(yv)/(2*dy)
      cosdx2iv  = (1.0/(cos(yv)*dx))**2
      sincos2iv = sin(yv)/(dx*cos(yv)**2)
      DO k=1,l
        DO j=1,m-1
          DO i=1,n
            frhv(i,j,k) = cosdx2iv(j)*(vve(i,j,k)-2*v(i,j,k)
     +                              +vvw(i,j,k))
     +                  +   rdy2i    *(v(i,j+1,k)-2*v(i,j,k)+v(i,j-1,k))
     +                  -  tandyiv(j)*(v(i,j+1,k)-v(i,j-1,k))
     +                  -   cos2iv(j)* v(i,j,k)
     +                  + sincos2iv(j)*( uvne(i,j,k)+uvse(i,j,k)
     +                                  -uvnw(i,j,k)-uvsw(i,j,k) )
*            frhv(i,j,k) = frhv(i,j,k)*hvv(i,j,k,1)
          ENDDO
        ENDDO
      ENDDO
*
      call u_in_fu(u,v,frhu,frhv,frh)
*
      ifrh = eh * int_vol(frh,ist,iend,jst,jend,0)
*
*
*  5)  Coriolis effects C
*
      coru  = 0.0
      corv  = 0.0
      cor   = 0.0
      icor  = 0.0
*
*     Extra minus-sign for transfer to right-hand side
*
      DO k=1,l
        DO j=1,m
          fact = 0.25*sin(y(j))
          DO i=1,n
            coru(i,j,k) = fact*(vune(i,j,k)+vuse(i,j,k)
     +                         +vunw(i,j,k)+vusw(i,j,k))
          ENDDO
          coru(0,j,k) = coru(n,j,k)
        ENDDO
        DO j=1,m-1
          fact = 0.25*sin(yv(j))
          DO i=1,n
            corv(i,j,k) = - fact*(uvne(i,j,k)+uvse(i,j,k)
     +                           +uvnw(i,j,k)+uvsw(i,j,k))
          ENDDO
        ENDDO
      ENDDO
*
      call u_in_fu(u,v,coru,corv,cor)
*
      icor = rota * int_vol(cor,ist,iend,jst,jend,0)
*
*
*  6)  Non-linear terms
*
      nln  = 0.0
      inln = 0.0
      uux  = 0.0
      uvx  = 0.0
      vuy  = 0.0
      vvy  = 0.0
      wuz  = 0.0
      wvz  = 0.0
      uvt  = 0.0
      u2t  = 0.0
*
      tdzi     = 1.0/(4*dz*dfzT)
      cos2dxi  = 1.0/(cos(y)*4*dx)
      cos2dyi  = 1.0/(cos(y)*4*dy)
      DO k=1,l
        DO j=1,m
          c1 = cos(yv(j))
          c2 = cos(yv(j-1))
          c3 = tan(y(j))/4.
          DO i=1,n-1
            uux(i,j,k) = cos2dxi(j)*(u(i+1,j,k)*u(i+1,j,k)
     +                              -u(i-1,j,k)*u(i-1,j,k))
          ENDDO
          DO i=1,n
            vuy(i,j,k) = cos2dyi(j)*( c1*(v(i,j  ,k)+v(i+1,j  ,k))
     +                                  *(uun(i,j,k)+u(i  ,j  ,k))
     +                              - c2*(v(i,j-1,k)+v(i+1,j-1,k))
     +                                  *(u(i,j  ,k)+uus(i,j,k)) )
            wuz(i,j,k) = tdzi(k)*( (w(i,j,k  ) + w(i+1,j,k  ))
     +                            *(uut(i,j,k) + u(i  ,j,k  ))
     +                           - (w(i,j,k-1) + w(i+1,j,k-1))
     +                            *(u(i,j,k  ) + uub(i,j,k  )) )
            vsum       =  v(i,j,k)+v(i+1,j,k)+v(i,j-1,k)+v(i+1,j-1,k)
            uvt(i,j,k) =  c3*vsum*u(i,j,k)
          ENDDO
          vuy(0,j,k) = vuy(n,j,k)
          wuz(0,j,k) = wuz(n,j,k)
          uvt(0,j,k) = uvt(n,j,k)
        ENDDO
      ENDDO
*
      tdzi     = 1.0/(4*dz*dfzT)
      cos2dxiv = 1.0/(cos(yv)*4*dx)
      cos2dyiv = 1.0/(cos(yv)*4*dy)
      DO k=1,l
        DO j=1,m-1
          c1 = cos(yv(j+1))
          c2 = cos(yv(j-1))
          c3 = tan(yv(j))/16.
          DO i=1,n
            vvy(i,j,k) = cos2dyiv(j)*(c1*v(i,j+1,k)*v(i,j+1,k)
     +                               -c2*v(i,j-1,k)*v(i,j-1,k))
          ENDDO
        ENDDO
        DO j=1,m
          c3 = tan(yv(j))/16.
          DO i=1,n
            uvx(i,j,k) = cos2dxiv(j)*((u(i  ,j+1,k)+u(i  ,j,k))
     +                              *(vve(i,j  ,k)+v(i  ,j,k))
     +                             - (u(i-1,j+1,k)+u(i-1,j,k))
     +                              *(v(i  ,j  ,k)+vvw(i,j,k)) )
            wvz(i,j,k) = tdzi(k)*( (w(i,j,k  ) + w(i,j+1,k  ))
     +                            *(vvt(i,j,k) + v(i,j  ,k  ))
     +                           - (w(i,j,k-1) + w(i,j+1,k-1))
     +                            *(v(i,j,k  ) + vvb(i,j  ,k)) )
            usum1 = u(i-1,j+1,k)+u(i,j+1,k) +u(i-1,j,k) +u(i,j,k)
            usum2 = uvnw(i,j,k) +uvne(i,j,k)+uvsw(i,j,k)+uvse(i,j,k)
            u2t(i,j,k) = c3*usum1*usum2
          ENDDO
        ENDDO
      ENDDO
*
*  Extra minus-sign for transfer to right-hand side
*
      nlnu = - (uux + vuy + wuz - uvt)
      nlnv = - (uvx + vvy + wvz + u2t)

      call u_in_fu(u,v,nlnu,nlnv,nln)
*
      inln = epsr * int_vol(nln,ist,iend,jst,jend,0)
*
*
* Sum of sources/sinks should be zero
*
      isum = iugp+iwsw+ifrv+ifrh+icor+inln
*
      write(80,980)xl,iugp,iwsw,ifrv,ifrh,icor,inln,isum,psimax
*
*      call write_momu(ugpu,wswu,frvu,frhu,coru,nlnu)
*      call write_momv(ugpv,wswv,frvv,frhv,corv,nlnv)

 980  format(10e12.4)
      END
******************************************************************
      SUBROUTINE u_in_fu(u,v,fu,fv,ufu)
      implicit none
      include 'usr.com'
*
* Inner product of (u,v) and terms from momentum equation (fu,fv)
*
      real      u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1)
      real      fu(0:n,0:m+1,0:l+1),fv(0:n+1,0:m,0:l+1)
      real      ufu(n,m,l)
      real      c1,c2,ci
      integer   i,j,k
*
      ufu = 0.0
      DO k=1,l
        DO j=1,m
          c1 = cos(yv(j))
          c2 = cos(yv(j-1))
          ci = 1./(c1+c2)
          DO i=1,n
            ufu(i,j,k)  = 0.5 * (    u(i-1,j,k)*fu(i-1,j,k) 
     +                  +            u(i  ,j,k)*fu(i  ,j,k) )
     +                  + ci  * ( c2*v(i,j-1,k)*fv(i,j-1,k) 
     +                  +         c1*v(i,j  ,k)*fv(i,j  ,k) )
          ENDDO
        ENDDO
      ENDDO
*
      return
      end
******************************************************************
      SUBROUTINE write_momu(ugpu,wswu,frvu,frhu,coru,nlnu)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      ugpu(0:n,0:m+1,0:l+1),wswu(0:n,0:m+1,0:l+1)
      real      frvu(0:n,0:m+1,0:l+1),frhu(0:n,0:m+1,0:l+1)
      real      coru(0:n,0:m+1,0:l+1),nlnu(0:n,0:m+1,0:l+1)
*     LOCAL
      real      momu(0:n,0:m+1,0:l+1)
      real      sigma,eh,ev,rota,epsr
      real      maxmom
      integer   i,j,k,imax,jmax,kmax
*
      epsr   = par(ROSB)
      ev     = par(EK_V)
      eh     = par(EK_H)
      sigma  = par(WIND)*par(AL_T)/(dz*dfzT(l))
      rota   = 0.0 ! par(ROTN) ! ATvS-Mix
      momu = ugpu+sigma*wswu+ev*frvu+eh*frhu+rota*coru+epsr*nlnu
*
      maxmom = 0.0
      imax = 0
      jmax = 0
      kmax = 0
*
      do k=1,l
        do j=1,m
          do i=1,n
          if (landm(i,j,k).eq.OCEAN .and. landm(i+1,j,k).ne.LAND) then
            write(96,996) i,j,k,momu(i,j,k),ugpu(i,j,k),
     +                      sigma*wswu(i,j,k),
     +                      ev*frvu(i,j,k),eh*frhu(i,j,k),
     +                      rota*coru(i,j,k),epsr*nlnu(i,j,k)
            if ( abs(momu(i,j,k)).gt.maxmom ) then 
              maxmom = abs(momu(i,j,k))
              imax = i
              jmax = j
              kmax = k
            endif
          endif
          enddo
        enddo
        write(96,*) 
      enddo
      write(96,997) 'zonal momentum equation largest residual : ',
     +                   maxmom, ' at ',imax,jmax,kmax
      write(96,*) 
*
 996  format(3i3,7e12.4)
 997  format(a,1e12.4,a,3i3)
*
      return
      end
******************************************************************
      SUBROUTINE write_momv(ugpv,wswv,frvv,frhv,corv,nlnv)
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      ugpv(0:n+1,0:m,0:l+1),wswv(0:n+1,0:m,0:l+1)
      real      frvv(0:n+1,0:m,0:l+1),frhv(0:n+1,0:m,0:l+1)
      real      corv(0:n+1,0:m,0:l+1),nlnv(0:n+1,0:m,0:l+1)
*     LOCAL
      real      momv(0:n+1,0:m,0:l+1)
      real      sigma,eh,ev,rota,epsr
      real      maxmom
      integer   i,j,k,imax,jmax,kmax
*
      epsr   = par(ROSB)
      ev     = par(EK_V)
      eh     = par(EK_H)
      sigma  = par(WIND)*par(AL_T)/(dz*dfzT(l))
      rota   = 0.0 !par(ROTN) ! ATvS-Mix
*
      momv = ugpv+sigma*wswv+ev*frvv+eh*frhv+rota*corv+epsr*nlnv
*
      maxmom = 0.0
      imax = 0
      jmax = 0
      kmax = 0
      do k=1,l
        do j=1,m
          do i=1,n
          if (landm(i,j,k).eq.OCEAN .and. landm(i,j+1,k).ne.LAND) then
            write(96,996) i,j,k,momv(i,j,k),ugpv(i,j,k),
     +                    sigma*wswv(i,j,k),
     +                    ev*frvv(i,j,k),eh*frhv(i,j,k),
     +                    rota*corv(i,j,k),epsr*nlnv(i,j,k)
            if ( abs(momv(i,j,k)).gt.maxmom ) then 
              maxmom = abs(momv(i,j,k))
              imax = i
              jmax = j
              kmax = k
            endif
          endif
          enddo
        enddo
      enddo
      write(96,997) 'meridional momentum equation largest residual : ',
     +                   maxmom, ' at ',imax,jmax,kmax
      write(96,*) 
*
 996  format(3i3,7e12.4)
 997  format(a,1e12.4,a,3i3) 
*
      return
      end
******************************************************************
      SUBROUTINE pot_exch(un,uw,psimax,xl)
*
* This routine calculates interaction terms between the basic state un
* and a disturbance uw
*
      implicit none
      include 'usr.com'
*     INPUT/OUTPUT
      real      un(ndim),uw(ndim),psimax,xl
*     LOCAL
      real      ihdif(3),ivdif(3),idadv(3,3),iconv(3),isurf(3),ipsum(3)
      real      ubar(0:n,0:m+1,0:l+1),vbar(0:n+1,0:m,0:l+1),
     +          wbar(0:n+1,0:m+1,0:l),pbar(0:n+1,0:m+1,0:l+1),
     +          tbar(0:n+1,0:m+1,0:l+1),sbar(0:n+1,0:m+1,0:l+1)
      real      u(0:n,0:m+1,0:l+1),v(0:n+1,0:m,0:l+1),
     +          w(0:n+1,0:m+1,0:l),p(0:n+1,0:m+1,0:l+1),
     +          t(0:n+1,0:m+1,0:l+1),s(0:n+1,0:m+1,0:l+1)
      real      rbar(0:n+1,0:m+1,0:l+1),r(0:n+1,0:m+1,0:l+1)
      integer   i,j,k
      real      hdift(n,m,l),hdifs(n,m,l),vdift(n,m,l),vdifs(n,m,l)
      real      dadvtx(n,m,l),dadvsx(n,m,l)
      real      dadvtz(n,m,l),dadvsz(n,m,l),convt(n,m,l),convs(n,m,l)
      real      surft(n,m),surfs(n,m),dadvty(n,m,l),dadvsy(n,m,l)
      real      lambda,ra,ph,pv,pvc,h1,h2,dum1,dum2,gamma,etat,bi
      real      tandyi(m),cosdx2i(m),csv(0:m),rdy2i(m),rdz2i
      real      cosdxi(m),rdyi(m),rdzi(l)
      real      advut(0:n,0:m+1,0:l+1),advus(0:n,0:m+1,0:l+1)
      real      advvt(0:n+1,0:m,0:l+1),advvs(0:n+1,0:m,0:l+1)
      real      advwt(0:n+1,0:m+1,0:l),advws(0:n+1,0:m+1,0:l)
      real      termt,terms
      integer west,north,east,south,bottom,top,center
      real te(0:n+1,0:m+1,0:l+1),se(0:n+1,0:m+1,0:l+1)
      real tw(0:n+1,0:m+1,0:l+1),sw(0:n+1,0:m+1,0:l+1)
      real tn(0:n+1,0:m+1,0:l+1),sn(0:n+1,0:m+1,0:l+1)
      real ts(0:n+1,0:m+1,0:l+1),sso(0:n+1,0:m+1,0:l+1)
      real tb(0:n+1,0:m+1,0:l+1),sb(0:n+1,0:m+1,0:l+1)
      real tto(0:n+1,0:m+1,0:l+1),st(0:n+1,0:m+1,0:l+1)
      real tbare(0:n+1,0:m+1,0:l+1),sbare(0:n+1,0:m+1,0:l+1)
      real tbarn(0:n+1,0:m+1,0:l+1),sbarn(0:n+1,0:m+1,0:l+1)
      real tbart(0:n+1,0:m+1,0:l+1),sbart(0:n+1,0:m+1,0:l+1)
      real rbarb(0:n+1,0:m+1,0:l+1),rbart(0:n+1,0:m+1,0:l+1)  
*     COMMON
      integer ist,iend,jst,jend
      common /integrate/ ist,iend,jst,jend
*     EXTERNAL FUNCTIONS
      real      int_vol,int_sur,hs
*
      ra     = par(RAYL)
      lambda = par(LAMB)
      ph     = par(PE_H)
      pv     = par(PE_V)
      pvc    = par(P_VC)
      bi     = par(BIOT)
*
      call usol(un,ubar,vbar,wbar,pbar,tbar,sbar)
      call usol(uw,u,v,w,p,t,s)
      rbar   = lambda*sbar - tbar
      r      = lambda*s    - t
*
* Initialize stencil points
*
      do k = 1, l
      do j = 1, m
      do i = 1, n
         west   = landm(i-1,j  ,k  )   ! 1
         south  = landm(i  ,j-1,k  )   ! 3
         center = landm(i  ,j  ,k  )   ! 4
         north  = landm(i  ,j+1,k  )   ! 5
         east   = landm(i+1,j  ,k  )   ! 7
         bottom = landm(i  ,j  ,k-1)   ! 8
         top    = landm(i  ,j  ,k+1)   ! 9
         if (center == OCEAN) then
            te(i,j,k) = t(i+1,j  ,k  )
            se(i,j,k) = s(i+1,j  ,k  )
            tw(i,j,k) = t(i-1,j  ,k  )
            sw(i,j,k) = s(i-1,j  ,k  )
            tn(i,j,k) = t(i  ,j+1,k  )
            sn(i,j,k) = s(i  ,j+1,k  )
            ts(i,j,k) = t(i  ,j-1,k  )
            sso(i,j,k)= s(i  ,j-1,k  )
            tto(i,j,k)= t(i  ,j  ,k+1)
            st(i,j,k) = s(i  ,j  ,k+1)
            tb(i,j,k) = t(i  ,j  ,k-1)
            sb(i,j,k) = s(i  ,j  ,k-1)
            tbare(i,j,k) = tbar(i+1,j  ,k  )
            sbare(i,j,k) = sbar(i+1,j  ,k  )
            tbarn(i,j,k) = tbar(i  ,j+1,k  )
            sbarn(i,j,k) = sbar(i  ,j+1,k  )
            tbart(i,j,k) = tbar(i  ,j  ,k+1)
            sbart(i,j,k) = sbar(i  ,j  ,k+1)
            rbart(i,j,k) = rbar(i  ,j  ,k+1)
            rbarb(i,j,k) = rbar(i  ,j  ,k-1)
            if (west == LAND) then
               tw(i,j,k)    = t(i,j,k)
               sw(i,j,k)    = s(i,j,k)
            endif
            if (east == LAND) then
               te(i,j,k)    = t(i,j,k)
               se(i,j,k)    = s(i,j,k)
               tbare(i,j,k) = tbar(i,j,k)
               sbare(i,j,k) = sbar(i,j,k)
            endif
            if (south == LAND) then
               ts(i,j,k)    = t(i,j,k)
               sso(i,j,k)   = s(i,j,k)
            endif
            if (north == LAND) then
               tn(i,j,k)    = t(i,j,k)
               sn(i,j,k)    = s(i,j,k)
               tbarn(i,j,k) = tbar(i,j,k)
               sbarn(i,j,k) = sbar(i,j,k)
            endif
            if (bottom == LAND) then
               tb(i,j,k)    = t(i,j,k)
               sb(i,j,k)    = s(i,j,k)
               rbarb(i,j,k) = rbar(i,j,k)
            endif
*           if (top == LAND) then
            if (k.eq.l) then
               tto(i,j,k)   = t(i,j,k)
               st(i,j,k)    = s(i,j,k)
               tbart(i,j,k) = tbar(i,j,k)
               sbart(i,j,k) = sbar(i,j,k)
               rbart(i,j,k) = rbar(i,j,k)
            endif
         endif
      enddo
      enddo
      enddo
*
*  1)  Horizontal diffusion:  < rho' div_h . grad_h rho' >
*
      hdift = 0.0
      hdifs = 0.0
      ihdif = 0.0
*
      csv      = cos(yv)
      rdy2i    = 1.0/(cos(y)*dy**2)
      tandyi   = tan(y)/(2*dy)
      cosdx2i  = (1.0/(cos(y)*dx))**2
      DO k=1,l
        DO j=1,m
          DO i=1,n
            hdift(i,j,k)= r(i,j,k)*( 
     +                    cosdx2i(j)*( te(i,j,k)+tw(i,j,k)-2*t(i,j,k))
     +                  +   rdy2i(j)*((tn(i,j,k)-t(i,j,k))*csv(j)
     +                               +(ts(i,j,k)-t(i,j,k))*csv(j-1)) )
*    +                  -  tandyi(j)*( tn(i,j,k)-ts(i,j,k)) )
            hdifs(i,j,k)= r(i,j,k)*( 
     +                    cosdx2i(j)*( se(i,j,k)+sw(i,j,k)-2*s(i,j,k))
     +                  +   rdy2i(j)*((sn(i,j,k)-s(i,j,k))*csv(j)
     +                               +(sso(i,j,k)-s(i,j,k))*csv(j-1)) )
*    +                  -  tandyi(j)*( sn(i,j,k)-sso(i,j,k)) )
          ENDDO
        ENDDO
      ENDDO
*
      ihdif(1) = - ph*int_vol(hdift,ist,iend,jst,jend,0)
      ihdif(2) =   ph*int_vol(hdifs,ist,iend,jst,jend,0)*lambda
      ihdif(3) =   ihdif(1) + ihdif(2)
*
*
*  2) Vertical diffusion:  < rho' d^2/dz^2 rho' >
*
      vdift = 0.0
      vdifs = 0.0
      ivdif = 0.0
  
      rdz2i = (1.0/dz)**2
      DO k=1,l
        h1 = rdz2i/(dfzT(k)*dfzW(k))
        h2 = rdz2i/(dfzT(k)*dfzW(k-1))
        DO j=1,m
          DO i=1,n
            vdift(i,j,k) = (  h1*(tto(i,j,k)-t(i,j,k))
     +                      + h2*(tb(i,j,k)-t(i,j,k)) )*r(i,j,k)
            vdifs(i,j,k) = (  h1*(st(i,j,k)-s(i,j,k))
     +                      + h2*(sb(i,j,k)-s(i,j,k)) )*r(i,j,k)
          ENDDO
        ENDDO
      ENDDO
  
      ivdif(1) = - pv*int_vol(vdift,ist,iend,jst,jend,0)
      ivdif(2) =   pv*int_vol(vdifs,ist,iend,jst,jend,0)*lambda
      ivdif(3) =   ivdif(1) + ivdif(2)
*
*  3) Convection (Enhanced vertical diffusion: old and new procedure)
*
      convt = 0.0
      convs = 0.0
      iconv = 0.0
*
      rdz2i = (1.0/dz)**2
      DO k=1,l
        h1 = 1.0/(dfzT(k)*dfzW(k))
        h2 = 1.0/(dfzT(k)*dfzW(k-1))
        DO j=1,m
          DO i=1,n
            dum1 = pvc * hs(rbart(i,j,k)-rbar(i,j,k)) 
     +           + pv  * (pv_adj(i,j,k)-1.)
            dum2 = pvc * hs(rbar(i,j,k)-rbarb(i,j,k))
     +           + pv  * (pv_adj(i,j,k-1)-1.)
            convt(i,j,k)= ( h1*dum1*(tto(i,j,k)-t(i,j,k))
     +                    + h2*dum2*(tb(i,j,k)-t(i,j,k)) )*r(i,j,k)
            convs(i,j,k)= ( h1*dum1*(st(i,j,k)-s(i,j,k))
     +                    + h2*dum2*(sb(i,j,k)-s(i,j,k)) )*r(i,j,k)
          ENDDO
        ENDDO
      ENDDO
*
      iconv(1) = - int_vol(convt,ist,iend,jst,jend,0)
      iconv(2) =   int_vol(convs,ist,iend,jst,jend,0)*lambda
      iconv(3) =   iconv(1) + iconv(2)
*
*
*  4) Advection: - < rho' div . (u rho) >
*
* Calculate advection on u,v en w-points first
*
      advut = 0.0
      advus = 0.0
      advvt = 0.0
      advvs = 0.0
      advwt = 0.0
      advws = 0.0
*
      DO k=1,l
        DO j=1,m
          DO i=1,n
            advut(i,j,k) = u(i,j,k)*(tbare(i,j,k)-tbar(i,j,k))
            advus(i,j,k) = u(i,j,k)*(sbare(i,j,k)-sbar(i,j,k))
            advvt(i,j,k) = v(i,j,k)*(tbarn(i,j,k)-tbar(i,j,k))*csv(j)
            advvs(i,j,k) = v(i,j,k)*(sbarn(i,j,k)-sbar(i,j,k))*csv(j)
            advwt(i,j,k) = w(i,j,k)*(tbart(i,j,k)-tbar(i,j,k))
            advws(i,j,k) = w(i,j,k)*(sbart(i,j,k)-sbar(i,j,k))
          ENDDO
        ENDDO
      ENDDO
*
      dadvtx = 0.0
      dadvty = 0.0
      dadvtz = 0.0
      dadvsx = 0.0
      dadvsy = 0.0
      dadvsz = 0.0
      idadv = 0.0
*
      rdyi    = 1.0/(2*cos(y)*dy)
      cosdxi  = 1.0/(2*cos(y)*dx)
      rdzi    = 1.0/(2*dfzT*dz)
      DO k=1,l
        DO j=1,m
          DO i=1,n
            dadvtx(i,j,k)= - r(i,j,k)*
     +         (advut(i,j,k)+advut(i-1,j,k))*cosdxi(j)
            dadvty(i,j,k)= - r(i,j,k)*
     +         (advvt(i,j,k)+advvt(i,j-1,k))*rdyi(j)
            dadvtz(i,j,k)= - r(i,j,k)*
     +         (advwt(i,j,k)+advwt(i,j,k-1))*rdzi(k)
*
            dadvsx(i,j,k)= - r(i,j,k)*
     +         (advus(i,j,k)+advus(i-1,j,k))*cosdxi(j)
            dadvsy(i,j,k)= - r(i,j,k)*
     +         (advvs(i,j,k)+advvs(i,j-1,k))*rdyi(j)
            dadvsz(i,j,k)= - r(i,j,k)*
     +         (advws(i,j,k)+advws(i,j,k-1))*rdzi(k)
          ENDDO
        ENDDO
      ENDDO
*
      idadv(1,1) = - int_vol(dadvtx,ist,iend,jst,jend,0)
      idadv(2,1) =   int_vol(dadvsx,ist,iend,jst,jend,0)*lambda
      idadv(3,1) =   idadv(1,1) + idadv(2,1)
      idadv(1,2) = - int_vol(dadvty,ist,iend,jst,jend,0)
      idadv(2,2) =   int_vol(dadvsy,ist,iend,jst,jend,0)*lambda
      idadv(3,2) =   idadv(1,2) + idadv(2,2)
      idadv(1,3) = - int_vol(dadvtz,ist,iend,jst,jend,0)
      idadv(2,3) =   int_vol(dadvsz,ist,iend,jst,jend,0)*lambda
      idadv(3,3) =   idadv(1,3) + idadv(2,3)
*
*
*  5) Surface forcing:  < bi rho' g(z) >
*
      surft = 0.0
      surfs = 0.0
      isurf = 0.0
*
      k=l
      DO j=1,m
        DO i=1,n
          surft(i,j) = - bi*TRES*r(i,j,l)*t(i,j,l)
          surfs(i,j) = - bi*SRES*r(i,j,l)*s(i,j,l)
        ENDDO
      ENDDO
*
      isurf(1) = - int_sur(surft,ist,iend,jst,jend,0)*dz*dfzT(l)
      isurf(2) =   int_sur(surfs,ist,iend,jst,jend,0)*dz*dfzT(l)*lambda
      isurf(3) =   isurf(1) + isurf(2)
*
      do i=1,3
        ipsum(i) = ihdif(i)+ivdif(i)+idadv(i,1)+idadv(i,2)+idadv(i,3)
     +                     +iconv(i)+isurf(i)
      enddo
*
      write(82,982) xl,ihdif,ivdif,idadv,iconv,isurf,ipsum,psimax
*
 982  format(26e12.4)
*
      END
