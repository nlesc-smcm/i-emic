#include "fdefs.h"

!****************************************************************************
SUBROUTINE forcing
  !     Shape of wind-forcing and buoyancy forcing
  use m_usr
  use m_atm
  implicit none

  real    sigma, etabi, gamma
  real    wfun, temfun, salfun
  real    temcor, salcor, adapted_salcor, spertcor
  integer i, j, k, row
  real    qfun2(n,m), fsint
  integer find_row2
  real :: max_internal_forcing

  if (iout.eq.0) then
     write(f99,*) '===========FORCING============'
  endif

  Frc   = 0.0
  ! ------------------------------------------------------------------
  ! Determine wind forcing
  ! ------------------------------------------------------------------

  sigma = par(COMB)*par(WIND)*par(AL_T)
  
  if (iza.eq.2) then        ! idealized wind forcing
     do j=1,m
        do i=1,n
           taux(i,j) = wfun(xu(i),yv(j),1)
           tauy(i,j) = wfun(xu(i),yv(j),2)
        enddo
     enddo
  endif
  
  do j = 1,m-1
     do i = 1, n
        Frc(find_row2(i,j,l,UU)) = sigma * taux(i,j)
        Frc(find_row2(i,j,l,VV)) = sigma * tauy(i,j)
     enddo
  enddo

  ! ------------------------------------------------------------------
  ! Determine atmospheric temperature forcing
  ! ------------------------------------------------------------------

  etabi = par(COMB)*par(TEMP)*(1 - TRES + TRES*par(BIOT))

  ! idealized temperature forcing
  if ((ite.eq.1).and.(coupled_atm.eq.0)) then
     do j=1,m
        do i=1,n
           tatm(i,j) = temfun(x(i),y(j))
        enddo
     enddo
  endif

  ! correct for nonzero flux
  if ((TRES.eq.0).and.(coupled_atm.eq.0)) then
     call qint(tatm,temcor)
  else
     temcor = 0.0
  end if

  do j=1,m
     do i=1,n
        if (la > 0) then ! coupling to atmosphere
           Frc(find_row2(i,j,l+1,TT)) = &
                par(COMB) * par(SUNP) * (suna(j) - amua)
           Frc(find_row2(i,j,l,TT)) = &
                par(COMB) * par(SUNP) * suno(j) * (1 - landm(i,j,l))
           
        else if (coupled_atm.eq.1) then ! coupled externally
           Frc(find_row2(i,j,l,TT)) = &
                ( par(COMB) * par(SUNP) * suno(j)  &
                +  Ooa * tatm(i,j) &
                +  par(COMB)*par(TEMP)*lvsc * qdim * eta * qatm(i,j)   & ! latent heat 
                -  par(COMB)*par(TEMP)*lvsc * eo0 &                      ! latent heat
                ) * (1 - landm(i,j,l))

        else  ! ocean-only
           Frc(find_row2(i,j,l,TT)) = etabi * ( tatm(i,j) - temcor )
        endif
     enddo
  enddo

  ! write(*,*) ' tatm(n/2,m/2)=', tatm(n/2,m/2), &
  !     ' qatm(n/2,m/2)=',  qatm(n/2,m/2)

  ! ------------------------------------------------------------------
  ! Determine salinity forcing
  ! ------------------------------------------------------------------

  if (coupled_atm.eq.1) then
     gamma = par(COMB) * par(SALT) 
  else 
     gamma = par(COMB)*par(SALT)*(1 - SRES + SRES*par(BIOT))
  endif

  if (its.eq.1) then  ! idealized salinity forcing
     do j=1,m
        do i=1,n
           emip(i,j) = salfun(x(i),y(j))
        enddo
     enddo
     call qint(emip,  salcor)
  else
     salcor = 0.0;
  endif

  if (SRES.eq.0.and.coupled_atm.eq.0) then   ! correct for nonzero flux

     call qint(adapted_emip, adapted_salcor)
     call qint(spert, spertcor)

     !write(*,*) 'salcor=',salcor, ' adapted_salcor=', adapted_salcor, 'spertcor=', spertcor
     !write(*,*) 'emip(10,10)', emip(10,10)
     !write(*,*) 'adapted_emip(10,10)', adapted_emip(10,10)
     !write(*,*) 'spert(10,10)', spert(10,10)
     
  else
     adapted_salcor = 0.0
     spertcor = 0.0
  end if
  
  do j=1,m
     do i=1,n
        ! nus*qdim*(E-P) without the sst dependency, which is taken care of in usrc.F90
        if (coupled_atm.eq.1) then
           Frc(find_row2(i,j,l,SS)) = gamma * nus * qdim * &
                ( -eta * qatm(i,j) - pfield(i,j) ) &
                * (1 - landm(i,j,l))
        else
           Frc(find_row2(i,j,l,SS)) = gamma * (1 - par(HMTP)) * ( emip(i,j) - salcor ) + &
                gamma * par(HMTP) * ( adapted_emip(i,j) - adapted_salcor ) + &
                par(SPER) * (1 - SRES + SRES*par(BIOT)) * ( spert(i,j) - spertcor )
        end if
     enddo
  enddo

  ! ------------------------------------------------------------------
  ! Determine forcing in the z-direction
  ! ------------------------------------------------------------------

  max_internal_forcing=0.0
  !OPEN(50,FILE='frc.txt',STATUS='unknown')
  !REWIND(50)

  do k=1,l-1
     do j=1,m
        do i=1,n
           Frc(find_row2(i,j,k,WW)) = -par(COMB)*(1-landm(i,j,k))*par(RAYL)*  &
                (par(LAMB)*(internal_salt(i,j,k) + internal_salt(i,j,k+1))/2. &
                -(internal_temp(i,j,k)+internal_temp(i,j,k+1))/2.)
           max_internal_forcing = max(max_internal_forcing, abs(Frc(find_row2(i,j,k,WW))))
           !WRITE (50,*) Frc(find_row2(i,j,k,WW)),internal_temp(i,j,k)
        enddo
     enddo
  enddo

  !CLOSE(50)

  write(f99,*) 'par(RAYL): ',par(RAYL)
  write(f99,*) 'par(COMB): ',par(COMB)
  write(f99,*) 'par(LAMB): ',par(LAMB)
  write(f99,*) 'maximum of internal forcing on w: ',max_internal_forcing

end subroutine forcing

!****************************************************************************
SUBROUTINE windfit
  USE m_itplbv
  use m_global
  implicit none

  !     CONSTANT
  integer nx,ny
  parameter(nx=145,ny=72)
  ! parameter(nx=360,ny=180)
  integer ::lwrk,liwrk
  !     LOCAL
  integer i,j,ifail,px,py,iwrk(n+nx)
  real    wrk2(4*n+nx+4)
  real    xx(nx),yy(ny),ff1(nx*ny),ff2(nx*ny)
  real    lambda(nx+4),mu(ny+4),cspl(nx*ny),wrk( (nx+6)*(ny+6) )
  real    dumx(n*m),dumy(n*m)
  real    xh(n),yh(m)
  real    xi(n*m), yi(n*m)
  real :: tmax

  lwrk=4*n+nx+4
  liwrk=n+nx

  open(10,file=topdir//'wind/trtau.dat',action='read')
  ! open(10,file=topdir//'cesm/wind_38Ma.txt',action='read')
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

#ifndef HAVE_NAG
  ! ACN alternative interpolation if nag-library is not available:
  do i = 1,n
     do j = 1,m
        xi(m*(i-1)+j) = xh(i)
        yi(m*(i-1)+j) = yh(j)
     enddo
  enddo
  call itplbv(f99,ny,nx,yy,xx,ff1,n*m,yi,xi,dumx)
  call itplbv(f99,ny,nx,yy,xx,ff2,n*m,yi,xi,dumy)
#else
  ! ACN original interpolation:
  call e01daf(nx,ny,xx,yy,ff1,px,py,lambda,mu,cspl,wrk,ifail)
  call e02dff(n,m,px,py,xh,yh,lambda,mu,cspl,dumx,wrk2,lwrk,&
       &            iwrk,liwrk,ifail)

  call e01daf(nx,ny,xx,yy,ff2,px,py,lambda,mu,cspl,wrk,ifail)
  call e02dff(n,m,px,py,xh,yh,lambda,mu,cspl,dumy,wrk2,lwrk,&
       &            iwrk,liwrk,ifail)
#endif
  do i=1,n
     do j=1,m
        taux(i,j) = dumx(m*(i-1) + j)
        tauy(i,j) = dumy(m*(i-1) + j)
     enddo
  enddo
  tmax = max(maxval(taux),maxval(tauy))
  !      call write_forcing('taux',taux,35)
  !      call write_forcing('tauy',tauy,36)

  write(f99,*) 'fit of windstress done, tmax = ', tmax

end subroutine windfit

!****************************************************************************
! put wind field for a certain month (1-12) into arrays taux, tauy in
! module m_global. This is used for the seasonal cycle problem.
!****************************************************************************

SUBROUTINE windfit_monthly(month)
  USE m_itplbv
  use m_global
  implicit none

  integer, intent(in) :: month

  !     CONSTANT
  integer nx,ny
  parameter(nx=144,ny=73)
  integer ::lwrk,liwrk
  !     LOCAL
  integer i,j,k,ifail,px,py,iwrk(n+nx)
  real    wrk2(4*n+nx+4)
  real    xx(nx),yy(ny),ff1(nx*ny),ff2(nx*ny)
  real    lambda(nx+4),mu(ny+4),cspl(nx*ny),wrk( (nx+6)*(ny+6) )
  real    dumx(n*m),dumy(n*m)
  real    xh(n),yh(m)
  real    xi(n*m), yi(n*m)
  real :: delx,dely
  real :: tmax

  lwrk=4*n+nx+4
  liwrk=n+nx

  delx=360.0/nx
  dely=180.0/(ny-1)

  do i=1,nx
     xx(i) = dble(i-1)*delx
  enddo
  do j=1,ny
     yy(j) = -90.0 + dble(j-1)*dely
  enddo
  xx = xx*pi/180.
  yy = yy*pi/180.

  open(10,file=topdir//'wind/trenberth/motaux.txt',action='read')
  open(11,file=topdir//'wind/trenberth/motauy.txt',action='read')
  ! note: we should probably skip over the unwanted months instead
  ! of reading everything every time
  do k=1,month
     do j=1,ny
        read(10,*) ff1(j:j+(nx-1)*ny:ny)
        read(11,*) ff2(j:j+(nx-1)*ny:ny)
     end do
  end do
  close(10)
  close(11)

  ! rest of sub is like original 'windfit'

  do i=1,n
     xh(i) = xu(i)
  enddo
  do j=1,m
     yh(j) = yv(j)
  enddo
  ifail = 0

#ifndef HAVE_NAG
  ! ACN alternative interpolation if nag-library is not available:
  do i = 1,n
     do j = 1,m
        xi(m*(i-1)+j) = xh(i)
        yi(m*(i-1)+j) = yh(j)
     enddo
  enddo
  call itplbv(f99,ny,nx,yy,xx,ff1,n*m,yi,xi,dumx)
  call itplbv(f99,ny,nx,yy,xx,ff2,n*m,yi,xi,dumy)
#else
  ! ACN original interpolation:
  call e01daf(nx,ny,xx,yy,ff1,px,py,lambda,mu,cspl,wrk,ifail)
  call e02dff(n,m,px,py,xh,yh,lambda,mu,cspl,dumx,wrk2,lwrk,&
       &            iwrk,liwrk,ifail)

  call e01daf(nx,ny,xx,yy,ff2,px,py,lambda,mu,cspl,wrk,ifail)
  call e02dff(n,m,px,py,xh,yh,lambda,mu,cspl,dumy,wrk2,lwrk,&
       &            iwrk,liwrk,ifail)
#endif
  do i=1,n
     do j=1,m
        taux(i,j) = dumx(m*(i-1) + j)
        tauy(i,j) = dumy(m*(i-1) + j)
     enddo
  enddo
  tmax = max(maxval(taux),maxval(tauy))
  !      call write_forcing('taux',taux,35)
  !      call write_forcing('tauy',tauy,36)

  write(f99,*) 'fit of monthly windstress done, tmax = ', tmax

end subroutine windfit_monthly

!****************************************************************************
SUBROUTINE read_spertm
  use m_global
  implicit none
  integer dum(0:n+1,0:m+1)
  integer i, j
  integer status

  open(unit=42,file='spertm_name.txt',status='old',err=995)
  read(unit=42,fmt='(A100)',iostat=status,end=10) spertmaskfile

10 continue
  close(42)

  write(*,*) '===========SalinityPert============================================'
  write(*,*) 'Salinity pert. mask is read in from file mkmask/'//trim(spertmaskfile)
  write(*,*) '===========SalinityPert============================================'


  open(unit=50,file=topdir//'mkmask/'//trim(spertmaskfile),status='old',err=995)
  do j = m+1, 0, -1
     read(50,'(362i1)') (dum(i,j),i=0,n+1)
  enddo
  
  close(50)
  
  do j=m,1,-1
     do i = 1,n
        spert(i,j) = real(1-dum(i,j))
     enddo
     write(*,'(100i1)') (1-dum(i,j), i=1,n)
  enddo

  return

995 write(*,*) 'WARNING: failed to read salinity perturbation mask from file'
  write(*,*) '         specified in spertm_name.txt. Either the name was not'
  write(*,*) '         specified by the C++ caller, or the file'
  write(*,*) '         (mkmask/'//trim(spertmaskfile)//') does not exist.'

end subroutine read_spertm

!****************************************************************************
real FUNCTION wfun(xx,yy,v1)
  use m_par
  use m_global ! we need the global value of ymin and ymax here!
  implicit none
  integer v1
  real    xx,yy,asym
  real :: tmax
  select case(v1)
  case(1)  ! x
     !     F. Bryans (1987) analytical profile:
     wfun = 0.2 - 0.8*sin(6*abs(yy)) - 0.5*(1 - tanh(10*abs(yy)))&
          -0.5*(1 - tanh(10*(pi/2 - abs(yy))))
  case(2)  ! y
     wfun = 0.0
  endselect
  tmax = 1.0
end FUNCTION wfun

!****************************************************************************
real FUNCTION temfun(xx,yy)
  use m_par
  use m_global ! we need the global value of ymin and ymax here!
  implicit none
  real    xx,yy
  if (ymin.ge.0.0) then ! Northern hemisphere
     temfun = cos(pi*(yy-ymin)/(ymax-ymin))
  else
     temfun = cos(pi*yy/ymax) + par(CMPR)*sin(pi*yy/ymax)
  end if
end FUNCTION temfun

!****************************************************************************
real FUNCTION salfun(xx,yy)
  use m_par
  use m_global ! we need the global value of ymin and ymax here!
  implicit none
  real    xx,yy
  real    subtr
  if (ymin.ge.0.0) then ! Northern hemisphere
     salfun = cos(pi*(yy-ymin)/(ymax-ymin))
     ! salfun=0.0
  else
     salfun = cos(pi*yy/ymax)/cos(yy)

  end if
end FUNCTION salfun

!**********************************************************************
SUBROUTINE qint(field,cor)
  ! Compute correction for nonzero flux
  use m_usr
  implicit none
  integer i,j
  real field(n,m),cor,sint
  ! This is an evil breach of concept, we call a C++ function from F90
  ! to compute the global integral:
  external thcm_forcing_integral
  call thcm_forcing_integral(field,y,landm,cor)
end SUBROUTINE qint

!**********************************************************************
SUBROUTINE tempfit
  USE m_itplbv
  use m_usr
  implicit none
  !     CONSTANT
  integer nx,ny
  parameter(nx=144,ny=73)
  integer ::lwrk,liwrk
  !     LOCAL
  integer i,j,ifail,px,py,iwrk(n+nx)
  real    wrk2(4*n+nx+4)
  real    xx(nx),yy(ny),t1(nx*ny)
  real    lambda(nx+4),mu(ny+4),cspl(nx*ny),wrk( (nx+6)*(ny+6) )
  real    dumx(n*m), tatmmax
  real    xi(n*m), yi(n*m)
  !
  lwrk=4*n+nx+4
  liwrk=n+nx
  open(10,file=topdir//'heat/nceptatm.dat')
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
  !
  ifail = 0
#ifndef HAVE_NAG
  ! ACN alternative interpolation if nag-library is not available:
  do i = 1,n
     do j = 1,m
        xi(m*(i-1)+j) = x(i)
        yi(m*(i-1)+j) = y(j)
     enddo
  enddo
  call itplbv(f99,ny,nx,yy,xx,t1,n*m,yi,xi,dumx)
#else
  ! ACN original interpolation:
  call e01daf(nx,ny,xx,yy,t1,px,py,lambda,mu,cspl,wrk,ifail)
  call e02dff(n,m,px,py,x,y,lambda,mu,cspl,dumx,wrk2,lwrk,&
       &            iwrk,liwrk,ifail)

#endif
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
  !
  write(f99,*) 'fit of temperature field done, tatmmax  = ', tatmmax

999 format(2i4)
998 format(5e16.8,i5,e14.6)

END SUBROUTINE tempfit

!**********************************************************************
SUBROUTINE fwfluxfit
  use m_usr
  implicit none
  integer, parameter :: nx = 180
  integer, parameter :: ny =  90
  integer, parameter :: nt =  12
  real*4 dat(nx, ny, 0:nt)
  real xx(0:nx), yy(0:ny)
  integer i, j, it, ii, jj
  real xi, yj, emipmax

  open(1,file=topdir//'salt/decdata.r4',form='unformatted',&
       status='unknown',access='direct',RECL=NX*NY)
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

        emip(i,j) = (1-xi)*(1-yj)*dat(ii  ,jj  ,0)&
             + (  xi)*(1-yj)*dat(ii+1,jj  ,0)&
             + (1-xi)*(  yj)*dat(ii  ,jj+1,0)&
             + (  xi)*(  yj)*dat(ii+1,jj+1,0)
     enddo
  enddo

  emipmax = maxval(emip)
  emip = emip / emipmax
  !
  write(f99,*) 'fit of fresh water flux done, emipmax  = ', emipmax

end SUBROUTINE fwfluxfit

!*****************************************************************
SUBROUTINE diagnose_restoring_sflux(un,filename)
  use m_usr
  implicit none
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
  !
  ! Construct internal and surface salt flux implied by relaxation
  !
  do k = 1,l
     !       amp = par(INTRS)
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
  !
  ! Write fluxes to file
  !
  write(f99,*) 'write diagnosed salt source field to: ',filename
  call write_internal_forcing(filename,sflux,41)

  return
END SUBROUTINE diagnose_restoring_sflux

!*****************************************************************
SUBROUTINE diagnose_restoring_tflux(un,filename)
  use m_usr
  implicit none
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
  !
  ! Construct internal and surface heat flux implied by relaxation
  !
  do k = 1,l
     !       amp = par(INTRT)
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
  !
  ! Write fluxes to file
  !
  write(f99,*) 'write diagnosed heat source field to: ',filename
  call write_internal_forcing(filename,tflux,41)

  return
END SUBROUTINE diagnose_restoring_tflux

!**********************************************************************
SUBROUTINE qint3(fl)
  !
  ! Subtracts net flux from flux field.
  ! Note that integral is according to dz(l)/dz(k) scaling of flux.
  !
  use m_usr
  implicit none
  real fl(n,m,l), sint, svol, int
  integer i,j,k

  call ssint(fl,1,l,sint,svol)

  int = sint/svol

  fl  = fl - int

  call ssint(fl,1,l,sint,svol)

end SUBROUTINE qint3

!**********************************************************************
SUBROUTINE ssint(fl,k0,k1,sint,svol)
  use m_usr
  implicit none
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

  write(f99,999) k0,k1,sint
999 format(' Net flux from layers ',i3,' to ',i3,' :',e12.4)

end SUBROUTINE ssint
