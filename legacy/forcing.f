*****************************************************************************
      SUBROUTINE forcing
*     shape of wind-forcing and buoyancy forcing 
      implicit none
      include 'usr.com'
      include 'atm.com'

      real    sigma, etabi, gamma
      real    salfun
      integer i, j, k, row
      real    wfun, hfun, qfun, qfun2(n,m), fsint,tauz(m)
      real    hom
      integer find_row2,ite,its

      if (iout.eq.0) then 
      write(99,*) '===========FORCING============'
      endif 
      if(SRES == 0) then  ! flux forcing 
        if (ifw.eq.1) then  ! only if idealized 
           do i=1,n
              do j=1,m
                 qfun2(i,j) = qfun(x(i),y(j))
              enddo
           enddo
        else ! flux will be read (-> emip) and can be perturbed here 
           call flux_perturbation(qfun2)
        endif
        call qint(qfun2,fsint) ! integral over domain is computed
        if (iout.eq.0) then 
        write(99,*) 'forcing with freshwater flux', fsint
        endif 
      else ! restoring salinity, profile is defined below 
         fsint = 0.        
         if (iout.eq.0) then 
         write(99,*) 'forcing with restoring salinity profile'
         endif
      endif
*
      hom = par(COMB) 
      sigma = hom*par(WIND)*par(AL_T)
      etabi = hom*par(TEMP)*(1- TRES + TRES*par(BIOT))
      gamma = hom*par(SALT)*(1- SRES + SRES*par(BIOT))
*
*  determine wind forcing field 
*
      SELECT CASE(iza) 
      CASE(0) ! from data 
      do j = 1, m
         do i = 1, n
            tx(i,j) = taux(i,j)
            ty(i,j) = tauy(i,j)
         enddo
      enddo
      CASE(1) ! zonally averaged
      DO j=1,m
	 tauz(j) = 0.0
         DO i=1,n
            tauz(j) = tauz(j) + taux(i,j)
         ENDDO
         tauz(j) = tauz(j)/n
      ENDDO
      do j = 1, m
         do i = 1, n
            tx(i,j) = tauz(j) 
            ty(i,j) = 0.0 
         enddo
      enddo
      CASE(2) ! idealized profile
      do j = 1, m
         do i = 1, n
            tx(i,j) = wfun(xu(i),yv(j),1) 
            ty(i,j) = wfun(xu(i),yv(j),2)                 
         enddo
      enddo
      END SELECT
      if (iout.eq.0) then 
      do i=1,n
        do j=1,m
           write(88,10) x(i)*180/pi,y(j)*180/pi,tx(i,j),ty(i,j)
        enddo 
      enddo
*
      write(99,*) 'you have chosen wind forcing option ',iza
      write(99,*) 'wind stress field written to unit 88'
      endif 
 10   format(1x,4(g12.5,1x))
* 
      Frc = 0.0
      do j = 1,m-1
         do i = 1, n                          
            Frc(find_row2(i,j,l,UU)) = sigma * tx(i,j) 
         enddo
      enddo   
      do j = 1, m-1                            
         do i = 1, n                         
            Frc(find_row2(i,j,l,VV)) = sigma * ty(i,j)
         enddo
      enddo  
c
      ite = 1
      its = 1
      do j = 1, m
         do i = 1, n
* top 
            if (la > 0 ) then ! coupling 
               ft(i,j) = par(COMB)*par(SUNP) * (suna(j) - amua)
               row = find_row2(i,j,l+1,TT)
               Frc(row) = ft(i,j)
               row = find_row2(i,j,l,TT)
               Frc(row) = par(COMB)* par(SUNP) * suno(j)
            else ! no coupling 
               row = find_row2(i,j,l,TT)
               Frc(row) = etabi * (ite*hfun(x(i),y(j)) + (1-ite)*tatm(i,j))
            endif
            if (SRES == 0) then  ! flux forcing 
                 if (ifw.eq.0) then ! flux forcing from data 
                 fs(i,j) = (1 - par(SPER))*(par(FPER)*(qfun2(i,j) - fsint) + 
     +                gamma*emip(i,j) )
                else ! idealized 
                 fs(i,j) = gamma*(qfun2(i,j) - fsint)
                endif
            else ! restoring forcing 
                fs(i,j) = gamma*(its*salfun(x(i),y(j))+(1-its)*emip(i,j)) 
            endif  
            row = find_row2(i,j,l,SS)
            Frc(row) = fs(i,j)
         enddo
      enddo
* 
      if (iout.eq.0) then 
      do i=1,n
        do j=1,m
           write(89,10) x(i)*180/pi,y(j)*180/pi,tatm(i,j),emip(i,j)
        enddo 
      enddo
*
      write(99,*) 'heat forcing and salinity forcing written to unit 89'
      write(99,*) '===========FORCING============'
      endif 

      iout = 1
      end
*****************************************************************************
      SUBROUTINE windfit
      USE m_itplbv
      implicit none
      include 'usr.com'
*     CONSTANT
      integer nx,ny
      parameter(nx=145,ny=72)
      integer,parameter::lwrk=4*n+nx+4,liwrk=n+nx
*     LOCAL
      integer i,j,ifail,px,py,iwrk(liwrk)
      real    wrk2(lwrk)
      real    xx(nx),yy(ny),ff1(nx*ny),ff2(nx*ny)
      real    lambda(nx+4),mu(ny+4),cspl(nx*ny),wrk( (nx+6)*(ny+6) )
      real    dumx(n*m),dumy(n*m)
      real    xh(n),yh(m)
      real    xi(n*m), yi(n*m)
      
      open(10,file=topdir//'wind/trtau.dat',action='read')
      read(10,*) 
      do i=1,nx
         read(10,*) xx(i)
      enddo
      do j=1,ny
         read(10,*) yy(j)
      enddo
      xx = xx*pi/180.
      yy = yy*pi/180.
      do i=1,nx
         do j=1,ny
            read(10,*) ff1(ny*(i-1)+j),ff2(ny*(i-1)+j)
         enddo
      enddo
      close(10)
      
      do i=1,n
        xh(i) = xu(i)
      enddo
      do j=1,m
        yh(j) = yv(j)
      enddo
      ifail = 0
      
! ACN alternative interpolation if nag-library is not available:
      do i = 1,n
         do j = 1,m
            xi(m*(i-1)+j) = xh(i)
            yi(m*(i-1)+j) = yh(j)
         enddo
      enddo
      call itplbv(99,ny,nx,yy,xx,ff1,n*m,yi,xi,dumx)
      call itplbv(99,ny,nx,yy,xx,ff2,n*m,yi,xi,dumy)
      
! ACN original interpolation:
!      call e01daf(nx,ny,xx,yy,ff1,px,py,lambda,mu,cspl,wrk,ifail)
!      call e02dff(n,m,px,py,xh,yh,lambda,mu,cspl,dumx,wrk2,lwrk,
!     &            iwrk,liwrk,ifail)

!      call e01daf(nx,ny,xx,yy,ff2,px,py,lambda,mu,cspl,wrk,ifail)
!      call e02dff(n,m,px,py,xh,yh,lambda,mu,cspl,dumy,wrk2,lwrk,
!     &            iwrk,liwrk,ifail)
!
      do i=1,n
         do j=1,m
            taux(i,j) = dumx(m*(i-1) + j)
            tauy(i,j) = dumy(m*(i-1) + j)
         enddo
      enddo
      tmax = max(maxval(taux),maxval(tauy))
!     call write_forcing('taux',taux,35)
!     call write_forcing('tauy',tauy,36)

      write(99,*) 'fit of windstress done, tmax = ', tmax

      end
*****************************************************************************
      real FUNCTION wfun(xx,yy,v1)
      implicit none
      include 'usr.com'
      integer v1
      real    xx,yy,asym
*      asym = par(CMPR) 
      select case(v1)
      case(1) ! x
*     F. Bryans (1987) analytical profile:
      wfun = 0.2 - 0.8*sin(6*abs(yy)) - 0.5*(1 - tanh(10*abs(yy)))
     +      -0.5*(1 - tanh(10*(pi/2 - abs(yy))))
*      wfun =-((1 - asym)*cos(2*pi * (yy - ymin) / (ymax-ymin))
*     +       + asym*cos(pi * (yy - ymin) / (ymax-ymin)) ) 
      case(2) ! y
         wfun = 0.0
      endselect
      tmax = 1.0
      end
*****************************************************************************
      real FUNCTION hfun(xx,yy)
      implicit none
      include 'usr.com'
      real    xx,yy
      if(ymin.gt.0.0) then  ! Only Northern hemisphere
        hfun = cos(pi*(yy-ymin)/(ymax-ymin))
      else
        hfun = cos(pi*yy/ymax)
      end if
      end
*****************************************************************************
      real FUNCTION salfun(xx,yy)
      implicit none
      include 'usr.com'
      real    xx,yy
      if(ymin.gt.0.0) then  ! Only Northern hemisphere
        salfun = cos(pi*(yy-ymin)/(ymax-ymin))
      else
        salfun = cos(pi*yy/ymax)
      end if
      end
*****************************************************************************
      real FUNCTION qfun(xx,yy)
      implicit none
      include 'usr.com'
      real    xx, yy,pps,qf1,qf2,pc
      pps = 2.0
*     qfun = cos(pps * pi * ((yy -ymin)/(ymax-ymin) -0.5)) 
*    + par(COMB) * sin(pps * pi * ((yy -ymin)/(ymax-ymin) -0.5))
*     qf1 = -(-cos (3*yy) +
*    +  2.38*exp(-(yy/0.2094395)**2))
*     qf2 = cos(pps * pi * ((yy -ymin)/(ymax-ymin) -0.5))  
*     pc = 0.0
*     pc = par(COMB)
*     qfun = qf2*(1 - pc) + pc*qf1
      qfun = cos(pi * (yy -ymin)/(ymax-ymin))  
      end
*****************************************************************************
      SUBROUTINE flux_perturbation(qfun2)
      implicit none
      include 'usr.com'
      real    yy,yc,qfun2(n,m)
      integer i, j, mask(0:n+1,0:m+1)
      character*23 filename

      if (rd_mask) then 
        filename = '_mkmask/mask.glo_atl2'
        open(unit=50,file=filename,err=995)
        do j = m+1, 0, -1
           read(50,'(98(i1,x))') (mask(i,j),i=0,n+1)
        enddo
        close(50)
      else 
        mask = 1
      endif 

      do j=1,m
        do i = 1,n
          qfun2(i,j) = real(1-mask(i,j))
        enddo
      enddo
*         IF (qfun2(i,j).gt.0) THEN
*            write(22,*) i,j,x(i)*180/pi,y(j)*180/pi
*         ENDIF
      return

 995  stop 'qfun2: open mask file'

      end
*****************************************************************************

      SUBROUTINE qint(qfun2,fsint)
*     compute average of surface salt flux perturbation
      implicit none
      include 'usr.com'
*     LOCAL
      integer i,j
      real qfun2(n,m),sint,fsint

      fsint = 0.0
      sint = 0.0
      do i=1,n
        do j=1,m
          fsint = qfun2(i,j) * cos(y(j)) * (1-landm(i,j,l)) + fsint
          sint = cos(y(j)) * (1-landm(i,j,l)) + sint
        enddo
      enddo
      fsint = fsint/sint

      end
***********************************************************************
      SUBROUTINE tempfit
      USE m_itplbv
      implicit none
      include 'usr.com'
*     CONSTANT
      integer nx,ny
      parameter(nx=144,ny=73)
      integer,parameter::lwrk=4*n+nx+4,liwrk=n+nx
*     LOCAL
      integer i,j,ifail,px,py,iwrk(liwrk)
      real    wrk2(lwrk)
      real    xx(nx),yy(ny),t1(nx*ny)
      real    lambda(nx+4),mu(ny+4),cspl(nx*ny),wrk( (nx+6)*(ny+6) )
      real    dumx(n*m), tatmmax
      real    xi(n*m), yi(n*m)
*
      open(10,file='/kaos_usr1/data/heat/nceptatm.dat')
c      open(10,file='/home/arie/data/heat/nceptatm.dat')
      read(10,*)
      DO i=1,nx
       read(10,*) xx(i)
      ENDDO
      DO j=1,ny
       read(10,*) yy(j)
      ENDDO
      xx = xx*pi/180.
      yy = yy*pi/180.
      DO i=1,nx
       DO j=1,ny
        read(10,*) t1(ny*(i-1)+j)
       ENDDO
      ENDDO
      close(10)
*
      ifail = 0
! ACN alternative interpolation if nag-library is not available:
      do i = 1,n
         do j = 1,m
            xi(m*(i-1)+j) = x(i)
            yi(m*(i-1)+j) = y(j)
         enddo
      enddo
      call itplbv(99,ny,nx,yy,xx,t1,n*m,yi,xi,dumx)
! ACN original interpolation:
!      call e01daf(nx,ny,xx,yy,t1,px,py,lambda,mu,cspl,wrk,ifail)
!      call e02dff(n,m,px,py,x,y,lambda,mu,cspl,dumx,wrk2,lwrk,
!     &            iwrk,liwrk,ifail)
!    
      DO i=1,n
       DO j=1,m
        tatm(i,j) = dumx(m*(i-1) + j)
       ENDDO
      ENDDO

      write(90,999) n,m
      write(90,998) xmin,xmax,ymin,ymax
      DO i=1,n
       DO j=1,m
         write(90,998) tatm(i,j)
       ENDDO
      ENDDO

      tatm = tatm - 273.15
      tatmmax = maxval(tatm)
      tatm = tatm / tatmmax
*
      write(99,*) 'fit of temperature field done, tatmmax  = ', tatmmax

 999  format(2i4)
 998  format(5e16.8,i5,e14.6)     

      END

***********************************************************************
      SUBROUTINE fwfluxfit
      implicit none
      include 'usr.com'
      integer, parameter :: nx = 180
      integer, parameter :: ny =  90
      integer, parameter :: nt =  12
      real*4 dat(nx, ny, 0:nt)
      real xx(0:nx), yy(0:ny)
      integer i, j, it, ii, jj
      real xi, yj, emipmax

      open(1,file='/kaos_usr1/data/salt/decdata.r4',form='unformatted',
c      open(1,file='/home/arie/data/salt/decdata.r4',form='unformatted',
     * status='unknown',access='direct',RECL=NX*NY)
      do IT=1,NT
         read(1,REC=IT)((dat(i,j,IT),i=1,nx),j=1,ny)
         dat(:,:,0) = dat(:,:,0) + dat(:,:,IT)
      enddo
      dat(:,:,0) = dat(:,:,0)/nt
      close(1)

      do i = 0, nx
         xx(i) = 2 * i - 1
      enddo
      do j = 0, ny
         yy(j) = 2 * j - 91
      enddo

      where(dat==-999.0) dat = 0. ! what else ?
      do j = 1, m
      do i = 1, n
         xi = 180.*x(i)/pi
         if(xi <= xx(0)) xi = xi + 360.
         if(xi > xx(nx)) xi = xi - 360.
         ii = floor((xi+1)/2.)
         xi = (xi - 2 * ii + 1)/2.
         yj = 180.*y(j)/pi
         jj = floor((yj+91)/2.)
         yj = (yj - 2 * jj + 91)/2.
         
         emip(i,j) = (1-xi)*(1-yj)*dat(ii  ,jj  ,0)
     >             + (  xi)*(1-yj)*dat(ii+1,jj  ,0)
     >             + (1-xi)*(  yj)*dat(ii  ,jj+1,0)
     >             + (  xi)*(  yj)*dat(ii+1,jj+1,0)
      enddo
      enddo

      emipmax = maxval(emip)
      emip = emip / emipmax
*
      write(99,*) 'fit of fresh water flux done, emipmax  = ', emipmax

      end
******************************************************************
      SUBROUTINE diagnose_restoring_sflux(un,filename)
      implicit none
      include 'usr.com'
      real, intent(in)   :: un(ndim)
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      integer i, j, k
      character*(*) filename
      real    sflux(n,m,l)
      real    amp,missing,p15i

      call usol(un,u,v,w,p,t,s)

      missing = -95.9999
      sflux = missing
      p15i = 1.0/par(SALT)
c
c Construct internal and surface salt flux implied by relaxation 
c
      do k = 1,l
*       amp = par(INTRS)
        amp = 1.0
        if (k.eq.l) amp = par(BIOT)
          do j = 1,m
            do i = 1,n
              if (landm(i,j,k).eq.OCEAN) then
                sflux(i,j,k) = amp*( fslev(i,j,k) - p15i*s(i,j,k) )
              endif
            enddo
         enddo
      enddo
c
c Write fluxes to file
c
      write(99,*) 'write diagnosed salt source field to: ',filename
      call write_internal_forcing(filename,sflux,41)

      return
      END

******************************************************************
      SUBROUTINE diagnose_restoring_tflux(un,filename)
      implicit none
      include 'usr.com'
      real, intent(in)   :: un(ndim)
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      integer i, j, k
      character*(*) filename
      real    tflux(n,m,l)
      real    amp,missing,p17i

      call usol(un,u,v,w,p,t,s)

      missing = -95.9999
      tflux = missing
      p17i = 1.0/par(TEMP)
c
c Construct internal and surface heat flux implied by relaxation 
c
      do k = 1,l
*       amp = par(INTRT)
        amp = 1.0
        if (k.eq.l) amp = par(BIOT)
        do j = 1,m
          do i = 1,n
            if (landm(i,j,k).eq.OCEAN) then
              tflux(i,j,k) = amp*( ftlev(i,j,k) - p17i*t(i,j,k) )
            endif
          enddo
        enddo
      enddo
c
c Write fluxes to file
c
      write(99,*) 'write diagnosed heat source field to: ',filename
      call write_internal_forcing(filename,tflux,41)

      return
      END
***********************************************************************
      SUBROUTINE qint3(fl)
*
* Subtracts net flux from flux field.
* Note that integral is according to dz(l)/dz(k) scaling of flux.
*
      implicit none
      include 'usr.com'
      real fl(n,m,l), sint, svol, int
      integer i,j,k

      call ssint(fl,1,l,sint,svol)

      int = sint/svol

      fl  = fl - int

      call ssint(fl,1,l,sint,svol)

      end
***********************************************************************

      SUBROUTINE ssint(fl,k0,k1,sint,svol)
      implicit none
      include 'usr.com'
      real fl(n,m,l), sint, svol
      integer i, j, k, k0, k1

      sint = 0.0
      svol = 0.0
      do k = k0,k1
        do j = 1,m
          do i = 1,n
            if (landm(i,j,k).eq.OCEAN) then
               sint = sint + cos(y(j))*fl(i,j,k)
               svol = svol + cos(y(j))
            endif
          enddo
        enddo
      enddo

      sint = sint*dx*dy*dfzT(l)*dz
      svol = svol*dx*dy*dfzT(l)*dz

      write(99,999) k0,k1,sint
 999  format(' Net flux from layers ',i3,' to ',i3,' :',e12.4)

      end
