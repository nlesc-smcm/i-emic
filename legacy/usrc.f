***************************************************************************
      SUBROUTINE init
      implicit none
      include 'usr.com'
*
      call grid
      call topofit
      call windfit
      call mixe
      call vmix_init	! ATvS-Mix
      if ((SRES.eq.0).and.(ifw.ne.1)) call read_forcing(emip,15) 
      if (SRES.eq.1) call levitus_sal
      if (TRES.eq.1) call levitus_sst
      call atmos_coef
      call forcing
*
      write(99,*) 'init done'

      end
*****************************************************************************
      SUBROUTINE mixe
      implicit none
      include 'usr.com'
      integer i,j,k
      do i=1,n
        do j=1,m
	   emix(i,j,0) = 0.0
	    do k=1,l-1
           emix(i,j,k)= - sin(pi*(x(i)-xmin)/(xmax-xmin))* 
     +    (cos((pi/2)*(y(j)-ymin)/(ymax-ymin)) + 0.2)
     +      *zw(k)
           enddo
           emix(i,j,l) = 0.0
        enddo
      enddo
*
      end
******************************************************************************
      SUBROUTINE matrix(un, sig)
      USE m_mat
*     construct the jacobian A and the 'mass' matrix B
*     Put B in coB,  A-sig*B in coA
      implicit none
      include 'usr.com'
      include 'mix.com'
      include 'atm.com'
      real un(ndim), sig, time0, time1
      integer i,j,k,row,ii,jj,kk,find_row2
      integer ix, iy, iz,ie

!     clean old arrays:

      coB  = 0
      Al   = 0
      begA = 0
      coA  = 0
      jcoA = 0
      call fillcolB
      call lin
      call nlin_jac(un)

! ATvS-Mix ---------------------------------------------------------------------              
      if (vmix_flag.ge.1) then
        call cpu_time(time0)
        if (vmix_out.gt.0) write (99,'(a26)')     'MIX| matrix...            '
        if (vmix_out.gt.0) write (99,'(a16,i10)') 'MIX|     fix:   ', vmix_fix
        if ((vmix_fix.eq.0).and.(vmix_flag.ge.2)) call vmix_control(un)
        if (vmix_out.gt.0) write (99,'(a16,i10)') 'MIX|     temp:  ', vmix_temp
        if (vmix_out.gt.0) write (99,'(a16,i10)') 'MIX|     salt:  ', vmix_salt
        if ((vmix_temp.eq.1).or.(vmix_salt.eq.1)) call vmix_jac(un)
        call cpu_time(time1)
        vmix_time=vmix_time+time1-time0
        write (99,'(a26, f10.3)') 'MIX|        ...matrix done', time1-time0
      endif
! --------------------------------------------------------------------- ATvS-Mix
      do k = 1, l+la
        do j = 1, m
          do i = 1, n
            do ii = 1,nun
               row = find_row2(i,j,k,ii)
               Al(i,j,k,5,ii,ii) = Al(i,j,k,5,ii,ii) - sig*coB(row)
            end do
          enddo
        enddo
      enddo

      if (frs) call frs_matrix(un,sig)
      call boundaries
      call assemble
!     call scale_mat_qz
      end
*****************************************************************************
      SUBROUTINE rhs(un,B)
*     construct the right hand side B
      USE m_mat
      implicit none
      include 'usr.com'
      include 'mix.com'
      include 'res.com'
      real    un(ndim),B(ndim), mix(ndim) ! ATvS-Mix
      real    Au(ndim), time0, time1
      integer i,j,k,k1,row,find_row2, mode

      mix=0.0
      Al   = 0
      begA = 0
      coA  = 0
      jcoA = 0
      call lin
      call nlin_rhs(un)
      call forcing
      if (frs) call frs_rhs(un)
      call boundaries
      call assemble
      call matAvec(un,Au)

! ATvS-Mix ---------------------------------------------------------------------              
      if (vmix_flag.ge.1) then
!        mode=vmix_fix
        call cpu_time(time0)
        if (vmix_out.gt.0) write (99,'(a26)')     'MIX| rhs...               '
        if (vmix_out.gt.0) write (99,'(a16,i10)') 'MIX|     fix:   ', vmix_fix
        if ((vmix_fix.eq.0).and.(vmix_flag.ge.2)) call vmix_control(un)
        if (vmix_out.gt.0) write (99,'(a16,i10)') 'MIX|     temp:  ', vmix_temp
        if (vmix_out.gt.0) write (99,'(a16,i10)') 'MIX|     salt:  ', vmix_salt
!        if ((vmix_temp.eq.1).or.(vmix_salt.eq.1)) call vmix_fun(un,mix,mode) !old mix_imp.f
        if ((vmix_temp.eq.1).or.(vmix_salt.eq.1)) call vmix_fun(un,mix)
        call cpu_time(time1)
        vmix_time=vmix_time+time1-time0
        write (99,'(a26,f10.3)') 'MIX|           ...rhs done', time1-time0
      endif
! --------------------------------------------------------------------- ATvS-Mix

c      write(99,*) "p0 = ", p0
      B = -Au - mix + Frc - p0*(1- par(RESC))*ures

      if(ires == 0) then
         DO i = 1, n
         DO j = 1, m
         DO k = 1, l
            DO k1 = 1,nun
               row = find_row2(i,j,k,k1)
               B(row) = B(row) * (1 - landm(i,j,k))
            ENDDO
         ENDDO
         ENDDO
         ENDDO
      endif
!     CALL scale_rhs_qz(B)
      end

*****************************************************************************
      SUBROUTINE scale_rhs_qz(B)
      USE m_mat
      implicit none
      include 'usr.com'
      include 'atm.com'
      !
      REAL B(ndim)
      ! LOCAL
      INTEGER row,i,j,k,XX
      REAL sf
      DO row = 1,ndim
         CALL findex(row,i,j,k,XX)
         IF (k.le.l) THEN
!             IF (XX.eq.WW) THEN
!                sf = dfzW(k)
!             ELSE
!                sf = dfzT(k)
!             END IF
!           IF (XX.eq.PP) THEN
!              sf = dfzW(k)
!           ELSE
!              sf = 1
!           END IF
           sf = dfzW(k)
         ELSE
            sf = 1
         END IF
         B(row) = sf*B(row)
      END DO
      
      END SUBROUTINE scale_rhs_qz
*****************************************************************************
      SUBROUTINE scale_mat_qz
      USE m_mat
      implicit none
      include 'usr.com'
      include 'atm.com'
      ! LOCAL
      INTEGER row,i,j,k,XX
      REAL sf
      DO row = 1,ndim
         CALL findex(row,i,j,k,XX)
         IF (k.le.l) THEN
!            IF (XX.eq.WW) THEN
!               sf = dfzW(k)
!            ELSE
!               sf = dfzT(k)
!            END IF
!           IF (XX.eq.PP) THEN
!              sf = dfzW(k)
!           ELSE
!              sf = 1
!           END IF
            sf = dfzW(k)
         ELSE
            sf = 1
         END IF
         coA(begA(row):begA(row+1)-1) = sf*coA(begA(row):begA(row+1)-1)
      END DO
      
      END SUBROUTINE scale_mat_qz
*****************************************************************************
      SUBROUTINE lin
      USE m_mat
*     Thermohaline equations
*     Produce local element matrices for linear operators
      implicit none
      include 'usr.com'
      include 'atm.com'
*     LOCAL
      real    u(n,m,l,np),ub(n,m,l,np),uxx(n,m,l,np),uyy(n,m,l,np),
     &        fv(n,m,l,np),uy(n,m,l,np),uzz(n,m,l,np)
      real    uxc(n,m,l,np),vyc(n,m,l,np),wzc(n,m,l,np) 
      real    px(n,m,l,np),py(n,m,l,np),pz(n,m,l,np) 
      real    ucsi(n,m,l,np),vcsi(n,m,l,np),uxs(n,m,l,np)
     &        ,vxs(n,m,l,np),vzz(n,m,l,np)
      real    v(n,m,l,np),vb(n,m,l,np),vxx(n,m,l,np),vyy(n,m,l,np),
     &        fu(n,m,l,np),vy(n,m,l,np)
      real    txx(n,m,l,np),tyy(n,m,l,np),tzz(n,m,l,np),tyc(n,m,l,np)
      real    tc(n,m,l,np),sc(n,m,l,np),tbc(n,m,l,np),tcb(n,m,l,np)
      real    yc(n,m,la,np),yc2(n,m,la,np),yxx(n,m,la,np),yyy(n,m,la,np)
      real    EH,EV,ph,pv,Ra,lambda, bi,ahcor
      real    hv(n,m,l,np),yadv(n,m,la,np)
      real    uxxc(n,m,l,np),uyyc(n,m,l,np),vxxc(n,m,l,np),vyyc(n,m,l,np)
      real    xes,rintb,rwint !, rintt ! ATvS-Mix
      equivalence (u, v), (uy, vy), (ucsi, vcsi), (uxx, vxx, txx, uxc)
      equivalence (uyy, vyy, tyy, vyc), (uzz, vzz, tzz, wzc)
      equivalence (uxs, vxs, tbc, tyc)
      equivalence (fu, fv, tc), (px, py, pz)

      EV     = par(EK_V) 
      EH     = par(EK_H)
      ph     = (1-par(MIXP))*par(PE_H)  !Isoneutral mixing (MIXP=1) is done in vmix_fun
      pv     = par(PE_V)
      lambda = par(LAMB)
      xes    = par(NLES)
      bi     = par(BIOT) 
      Ra     = par(RAYL)
c      rintt  = par(IFRICT)	! ATvS-Mix

      Al = 0.0
      call uderiv(1,ub)
      call uderiv(2,uxx)
      call uderiv(3,uyy)
      call uderiv(4,uzz)
      call uderiv(5,ucsi)
      call uderiv(6,vxs)
      call uderiv(7,u)
      call coriolis(1,fv)
      call gradp(1,px)
      Al(:,:,1:l,:,UU,UU) = -EH*(uxx+uyy+ucsi) -EV*uzz ! + rintt*u ! ATvS-Mix
      Al(:,:,1:l,:,UU,VV) = -fv - EH*vxs
      Al(:,:,1:l,:,UU,PP) = px

      call vderiv(1,vb )
      call vderiv(2,vxx)
      call vderiv(3,vyy)      
      call vderiv(4,vzz)
      call vderiv(5,vcsi)
      call vderiv(6,uxs)
      call vderiv(7,v)
      call coriolis(2,fu)
      call gradp(2,py)
      Al(:,:,1:l,:,VV,UU) = fu - EH*uxs
      Al(:,:,1:l,:,VV,VV) = -EH*(vxx + vyy+vcsi)  -EV*vzz !+ rintt*v ! ATvS-Mix 
      Al(:,:,1:l,:,VV,PP) = py

      call gradp(3,pz)
      call tderiv(6,tbc)
      Al(:,:,1:l,:,WW,PP) = pz
* 
      Al(:,:,1:l,:,WW,TT) = - Ra *(1. + xes*alpt1) * tbc/2. 
*
      Al(:,:,1:l,:,WW,SS) = lambda * Ra * tbc/2.

      call pderiv(1,uxc)
      call pderiv(2,vyc)
      call pderiv(3,wzc)
      Al(:,:,1:l,:,PP,UU) = uxc
      Al(:,:,1:l,:,PP,VV) = vyc
      Al(:,:,1:l,:,PP,WW) = wzc

      call tderiv(1,tc )
      call tderiv(2,sc )
      call tderiv(3,txx)
      call tderiv(4,tyy)
      call tderiv(5,tzz)
      call tderiv(7,tcb )
      if (la > 0) then
         Al(:,:,1:l,:,TT,TT) = - ph * (txx + tyy) - pv * tzz + Ooa*tc
         call yderiv(1,yc )
         call yderiv(4,yc2 )
         call yderiv(2,yxx)
         call yderiv(3,yyy)
         call yderiv(5,yadv)
         Al(:,:,l+1:l+la,5,UU,UU) = 1.
         Al(:,:,l+1:l+la,5,VV,VV) = 1.
         Al(:,:,l+1:l+la,5,WW,WW) = 1.
         Al(:,:,l+1:l+la,5,PP,PP) = 1.
         Al(:,:,l+1:l+la,5,SS,SS) = 1.
         Al(:,:,l+1:l+la,:,TT,TT) = - Ad * (yxx + yyy) - yc +
     +                              Aa*yadv + bmua*yc2
      else
      Al(:,:,1:l,:,TT,TT) = - ph * (txx + tyy) - pv * tzz + TRES*bi*tc
      endif
      Al(:,:,1:l,:,SS,SS) = - ph * (txx + tyy) - pv * tzz + SRES*bi*sc

      end
*****************************************************************************
      SUBROUTINE nlin_rhs(un)
      USE m_mat
*     Produce local matrices for nonlinear operators for calc of Rhs
      implicit none
      include 'usr.com'
      include 'mix.com'
      include 'atm.com'
*     IMPORT/EXPORT
      real    un(ndim)
*     LOCAL 
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      real    rho(0:n+1,0:m+1,0:l+la+1)
      real    utx(n,m,l,np),vty(n,m,l,np),wtz(n,m,l,np)
      real    usx(n,m,l,np),vsy(n,m,l,np),wsz(n,m,l,np)
      real    t2r(n,m,l,np),t3r(n,m,l,np)
      real    cat1(n,m,l,np),cas1(n,m,l,np)
      real    cat2(n,m,l,np),cas2(n,m,l,np)
      real    uux(n,m,l,np),uvy1(n,m,l,np),uwz(n,m,l,np),uvy2(n,m,l,np)
      real    uvx(n,m,l,np),vvy(n,m,l,np),vwz(n,m,l,np),ut2(n,m,l,np)
      real    bolt(n,m,la,np)
      real    lambda,epsr,Ra,xes,pvc1,pvc2, pv
      equivalence (utx, usx), (vty, vsy), (wtz, wsz)

      lambda = par(LAMB)
      epsr   = par(ROSB)
      Ra = par(RAYL)
      xes = par(NLES) 
      pvc1 = par(P_VC)
      pv = par(PE_V)
      pvc2 = pv*(1.0 - par(ALPC))*par(ENER)
      call usol(un,u,v,w,p,t,s)
      rho    = lambda*s - t *( 1 + xes*alpt1) - 
     +            xes*t*t*alpt2+xes*t*t*t*alpt3

      call unlin(1,uux,u,v,w)
      call unlin(3,uvy1,u,v,w)
      call unlin(5,uwz,u,v,w)
      call unlin(7,uvy2,u,v,w) 
      Al(:,:,1:l,:,UU,UU)  = Al(:,:,1:l,:,UU,UU) + epsr * (uux + uvy1 + uwz + uvy2)    
      call vnlin(1,uvx,u,v,w)
      call vnlin(3,vvy,u,v,w)
      call vnlin(5,vwz,u,v,w)
      call vnlin(7,ut2,u,v,w)
      Al(:,:,1:l,:,VV,UU) = Al(:,:,1:l,:,VV,UU) + epsr *ut2  
      Al(:,:,1:l,:,VV,VV) = Al(:,:,1:l,:,VV,VV) + epsr*(uvx + vvy + vwz)    
      call wnlin(2,t2r,t)
      call wnlin(4,t3r,t)
      Al(:,:,1:l,:,WW,TT) = Al(:,:,1:l,:,WW,TT) - Ra*xes*alpt2*t2r
     +                              + Ra*xes*alpt3*t3r
*
      call tnlin(3,utx,u,v,w,t,rho)
      call tnlin(5,vty,u,v,w,t,rho)
      call tnlin(7,wtz,u,v,w,t,rho)
      Al(:,:,1:l,:,TT,TT) = Al(:,:,1:l,:,TT,TT)+ utx+vty+wtz		! ATvS-Mix
      
      call tnlin(3,usx,u,v,w,s,rho)
      call tnlin(5,vsy,u,v,w,s,rho)
      call tnlin(7,wsz,u,v,w,s,rho)
      Al(:,:,1:l,:,SS,SS) = Al(:,:,1:l,:,SS,SS)+ usx+vsy+wsz		! ATvS-Mix

      write(99,*) 'nlin_rhs done'

      end
*****************************************************************************
      SUBROUTINE nlin_jac(un)
      USE m_mat
*     Produce local matrices for nonlinear operators for calc of Jacobian
      implicit none
      include 'usr.com'
      include 'mix.com'
      include 'atm.com'
*     IMPORT/EXPORT
      real    un(ndim)
*     LOCAL
      real    u(0:n  ,0:m,0:l+la+1), v(0:n, 0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      real    rho(0:n+1,0:m+1,0:l+la+1)
      real    urTx(n,m,l,np),Utrx(n,m,l,np),
     +        vrTy(n,m,l,np),Vtry(n,m,l,np)
      real    wrTz(n,m,l,np),Wtrz(n,m,l,np),
     +        cat1(n,m,l,np),cat2(n,m,l,np),
     +        cat3(n,m,l,np),cat4(n,m,l,np)
      real    t2r(n,m,l,np),t3r(n,m,l,np)
      real    urSx(n,m,l,np),Usrx(n,m,l,np),
     +        vrSy(n,m,l,np),Vsry(n,m,l,np)
      real    wrSz(n,m,l,np),Wsrz(n,m,l,np),
     +        cas1(n,m,l,np),cas2(n,m,l,np),
     +        cas3(n,m,l,np),cas4(n,m,l,np) 
      real    bolt(n,m,la,np)
      real    lambda,epsr,Ra,xes,pvc1,pvc2,pv
      real uux(n,m,l,np),uvy1(n,m,l,np),uwz(n,m,l,np),uvy2(n,m,l,np)
      real uvx(n,m,l,np),vvry(n,m,l,np),vwz(n,m,l,np),wvrz(n,m,l,np)
      real Urux(n,m,l,np),Urvy1(n,m,l,np),Urwz(n,m,l,np),Urvy2(n,m,l,np)
      real uVrx(n,m,l,np),Vrvy(n,m,l,np),Vrwz(n,m,l,np),Urt2(n,m,l,np)
      equivalence (urTx, urSx), (vrTy, vrSy), (wrTz, wrSz)
      equivalence (Utrx, Usrx), (Vtry, Vsry), (Wtrz, Wsrz)

      lambda = par(LAMB)
      epsr   = par(ROSB)
      Ra = par(RAYL) 
      xes = par(NLES)
      pvc1 = par(P_VC)
      pv = par(PE_V)
      pvc2 = pv*(1.0 - par(ALPC))*par(ENER)
      call usol(un,u,v,w,p,t,s)
      rho    = lambda*s - t *( 1 + xes*alpt1) - 
     +            xes*t*t*alpt2+xes*t*t*t*alpt3

      call unlin(2,Urux,u,v,w)
      call unlin(3,uvy1,u,v,w)
      call unlin(4,Urvy1,u,v,w)      
      call unlin(5,uwz,u,v,w)
      call unlin(6,Urwz,u,v,w)
      call unlin(7,uvy2,u,v,w) 
      call unlin(8,Urvy2,u,v,w)
      Al(:,:,1:l,:,UU,UU)  =  Al(:,:,1:l,:,UU,UU)+ epsr *(Urux +uvy1 + uwz + uvy2) 
      Al(:,:,1:l,:,UU,VV)  = Al(:,:,1:l,:,UU,VV)+ epsr * (Urvy1 + Urvy2) 
      Al(:,:,1:l,:,UU,WW)  = Al(:,:,1:l,:,UU,WW) + epsr * Urwz   
            
      call vnlin(1,uvx,u,v,w)
      call vnlin(2,uVrx,u,v,w)
      call vnlin(4,Vrvy,u,v,w)      
      call vnlin(5,vwz,u,v,w)
      call vnlin(6,Vrwz,u,v,w)
      call vnlin(8,Urt2,u,v,w)
      Al(:,:,1:l,:,VV,UU) =   Al(:,:,1:l,:,VV,UU)+ epsr *(Urt2 + uVrx)  
      Al(:,:,1:l,:,VV,VV) =   Al(:,:,1:l,:,VV,VV)+ epsr*(uvx+Vrvy+vwz)  
      Al(:,:,1:l,:,VV,WW) = Al(:,:,1:l,:,VV,WW) + epsr *Vrwz 
*
      call wnlin(1,t2r,t)
      call wnlin(3,t3r,t)
      Al(:,:,1:l,:,WW,TT) = Al(:,:,1:l,:,WW,TT) - Ra*xes*alpt2*t2r 
     +                              + Ra*xes*alpt3*t3r
*
      call tnlin(2,urTx,u,v,w,t,rho)
      call tnlin(3,Utrx,u,v,w,t,rho)
      call tnlin(4,vrTy,u,v,w,t,rho)
      call tnlin(5,Vtry,u,v,w,t,rho)
      call tnlin(6,wrTz,u,v,w,t,rho)
      call tnlin(7,Wtrz,u,v,w,t,rho)
      Al(:,:,1:l,:,TT,UU) = Al(:,:,1:l,:,TT,UU) + urTx
      Al(:,:,1:l,:,TT,VV) = Al(:,:,1:l,:,TT,VV) + vrTy
      Al(:,:,1:l,:,TT,WW) = Al(:,:,1:l,:,TT,WW) + wrTz
      Al(:,:,1:l,:,TT,TT) = Al(:,:,1:l,:,TT,TT) + Utrx + Vtry + Wtrz		! ATvS-Mix


      call tnlin(2,urSx,u,v,w,s,rho)
      call tnlin(3,Usrx,u,v,w,s,rho)
      call tnlin(4,vrSy,u,v,w,s,rho)
      call tnlin(5,VsRy,u,v,w,s,rho)
      call tnlin(6,wrSz,u,v,w,s,rho)
      call tnlin(7,Wsrz,u,v,w,s,rho)

      Al(:,:,1:l,:,SS,UU) = Al(:,:,1:l,:,SS,UU) + urSx
      Al(:,:,1:l,:,SS,VV) = Al(:,:,1:l,:,SS,VV) + vrSy
      Al(:,:,1:l,:,SS,WW) = Al(:,:,1:l,:,SS,WW) + wrSz
      Al(:,:,1:l,:,SS,SS) = Al(:,:,1:l,:,SS,SS) + Usrx + Vsry + Wsrz


      write(99,*) 'nlin_jac done'

      end
*****************************************************************************
      SUBROUTINE usol(un,u,v,w,p,t,s)
*     Go from un to u,v,t,h
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      real    un(ndim)
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
*     LOCAL
      integer i,j,k,row
*     EXTERNAL
      integer find_row2

      u = 0.0
      v = 0.0
      w = 0.0
      p = 0.0
      t = 0.0
      s = 0.0
      do k = 1, l+la 
         do j = 1, m
            do i = 1, n
               u(i,j,k) = un(find_row2(i,j,k,UU))
               v(i,j,k) = un(find_row2(i,j,k,VV))
               w(i,j,k) = un(find_row2(i,j,k,WW))
               p(i,j,k) = un(find_row2(i,j,k,PP))
               t(i,j,k) = un(find_row2(i,j,k,TT))
               s(i,j,k) = un(find_row2(i,j,k,SS))
            enddo
         enddo
      enddo
      do k = 1, l+la 
         do j = 1, m
            if (periodic) then
               u(0,j,k)  = u(N,j,k)
               v(0,j,k)  = v(N,j,k)
               w(N+1,j,k)= w(1,j,k)
               w(0,j,k)  = w(N,j,k)
               p(N+1,j,k)= p(1,j,k)
               p(0,j,k)  = p(N,j,k)
               t(N+1,j,k)= t(1,j,k)
               t(0,j,k)  = t(N,j,k)
               s(N+1,j,k)= s(1,j,k)
               s(0,j,k)  = s(N,j,k)
            else
               u(0,j,k)  = 0.0
               u(N,j,k)  = 0.0
               v(0,j,k)  = 0.0
               v(N,j,k)= 0.0
               p(0,j,k)  = 0.0
               p(N+1,j,k)= 0.0
               t(0,j,k)  = t(1,j,k)
               t(N+1,j,k)= t(N,j,k)
               s(0,j,k)  = s(1,j,k)
               s(N+1,j,k)= s(N,j,k)
            endif
         enddo
      enddo
      do k = 1, l+la 
         do i = 1, n
            u(i,0,k)  = 0.0
            u(i,M,k)= 0.0
            v(i,0,k)  = 0.0
            v(i,M,k)  = 0.0
            p(i,0,k)  = 0.0
            p(i,M+1,k)= 0.0
            t(i,0,k)  = t(i,1,k)
            t(i,M+1,k)= t(i,M,k)
            s(i,0,k)  = s(i,1,k)
            s(i,M+1,k)= s(i,M,k)
         enddo
      enddo
      do j = 1, m 
         do i = 1, n
            u(i,j,0)   = u(i,j,1)
            u(i,j,l+1) = u(i,j,l)
            v(i,j,0)   = v(i,j,1)
            v(i,j,l+1) = v(i,j,l)
            w(i,j,l)   = 0.0
            w(i,j,0)   = 0.0
            p(i,j,l+1) = 0.0
            p(i,j,0)   = 0.0
	    IF (la == 0) t(i,j,l+1) = t(i,j,l) 
            t(i,j,0)   = t(i,j,1) 
            s(i,j,l+1) = s(i,j,l) 
            s(i,j,0)   = s(i,j,1) 
         enddo
      enddo
 18   format(i4,i4,4(g12.4))
      do i=1,n
	do j=1,m
	  do k=1,l
	     if (landm(i,j,k).eq.1) then
		u(i,j,k) = 0.0
		v(i,j,k) = 0.0
		u(i-1,j,k) = 0.0
		v(i-1,j,k) = 0.0
		u(i,j-1,k) = 0.0
		v(i,j-1,k) = 0.0
		u(i-1,j-1,k) = 0.0
		v(i-1,j-1,k) = 0.0
             endif		  
          enddo
        enddo
      enddo

      end
*****************************************************************************
      SUBROUTINE solu(un,u,v,w,p,t,s)
*     Go from u,v,h,t to un
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      real    un(ndim)
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      integer find_row2
*     LOCAL
      integer i,j,k,row

      do k = 1, l+la 
      do j = 1, m
      do i = 1, n
         un(find_row2(i,j,k,UU)) = u(i,j,k)
         un(find_row2(i,j,k,VV)) = v(i,j,k)
         un(find_row2(i,j,k,WW)) = w(i,j,k)
         un(find_row2(i,j,k,PP)) = p(i,j,k)
         un(find_row2(i,j,k,TT)) = t(i,j,k)
         un(find_row2(i,j,k,SS)) = s(i,j,k)
      enddo
      enddo
      enddo

      end
*****************************************************************************
      SUBROUTINE stpnt(un)
      implicit none
      include 'usr.com'
      include 'atm.com'
      integer i,j,k,row,find_row2
      real    un(ndim)

********************************************************
*     PARAMETERS: 
********************************************************
* when data are used, tmax comes from windfit
* otherwise tmax comes from wfun

      par(AL_T)   = 0.1/(2*omegadim*rhodim*hdim*udim*dz*dfzT(l))
      par(RAYL)   = alphaT*gdim*hdim/(2*omegadim*udim*r0dim)       ! Ra
      par(EK_V)   = av/(2*omegadim*hdim*hdim)      ! E_V
      par(EK_H)   = ah/(2*omegadim*r0dim*r0dim)    ! E_H
      par(ROSB)   = udim/(2*omegadim*r0dim)        ! Rossby Number
      par(HMTP)   =  0.0
      par(SUNP)   =  0.0
      par(PE_H)   = kappah/(udim*r0dim)            ! P_H0
      par(PE_V)   = kappav*r0dim/(udim*hdim*hdim)  ! P_V0
      par(P_VC)   =  2.5e+04*par(PE_V)             ! P_VC
      par(LAMB)   = alphaS/alphaT  ! lambda
      par(SALT)   =  1.0e-01       ! gamma
      par(WIND)   =  1.0           ! wind h
      par(TEMP)   =  10.0          ! eta_T
      par(BIOT)   = r0dim/(75.*3600.*24.*udim)     ! nonlinearity in T,S equations
      par(COMB)   =  0.0           ! combined continuation
      par(NLES)   =  0.0
      par(CMPR)   =  0.0 
      par(ALPC)   =  1.0
      par(ENER)   =  1.0e+02
      par(MIXP)   =  0.0           ! neutral physics
      par(MKAP)   =  0.0           ! Gent-McWilliams
      par(SPL1)   =  2.0e+03       ! tanh in convective adjustment
      par(SPL2)   =  0.01          ! maximum admissible slope for neutral physics
      
      un = 0.0
*     do i=1,n
*       do j=1,m
*         do k=1,l
*           row = find_row2(i,j,k,TT)
*           un(row) = -bottem*z(k)
*           un(row) = bottem
*         enddo
*       enddo
*     enddo
*     do i= 5, ndim, 6
*        row = find_row(i,j,k)
*        un(i) = - (arad+brad*t0)/brad
*     enddo

      call vmix_par
      end

*******************************************************************
      SUBROUTINE atmos_coef
      implicit none
      include 'usr.com'
      include 'atm.com'
*     LOCAL
      real muoa,dzne
      integer i,j
      dzne = dz*dfzT(l)
      muoa = rhoa*ch*cpa*uw
      amua = (arad+brad*t0)/muoa
      bmua = brad/muoa
      Aa = uatm*rhoa*cpa*hdima/(r0dim*muoa)
      Ai = rhoa*hdima*cpa*udim/(r0dim*muoa)
      Ad = rhoa*hdima*cpa*d0/(muoa*r0dim*r0dim)
      As = sun0*(1 - c0)/(4*muoa)
      Os = sun0*c0*r0dim/(4*udim*hdim*dzne*rhodim*cp0)
      Ooa = muoa*r0dim/(udim*cp0*rhodim*hdim*dzne)
*      write(79,*) 'Aa   Ad   As   Os  Ooa  amua bmua'
*      write(79,*) Aa,Ad,As,Os,Ooa,amua,bmua
      DO j=1,m
*       albe(j) = 0.15+ 0.05 * cos (y(j))
         albe(j) = 0.3 
         dat(j) =  0.9 + 1.5 * exp(-12*y(j)*y(j)/pi)
         suno(j) = Os*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
         suna(j) = As*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
*         write(79,770) y(j)*180/pi,suna(j),suno(j),dat(j),albe(j)
* 770     format(1x,5(g12.5,1x))
      ENDDO
      
      DO j=0,m   
         davt(j) = 0.9 + 1.5 * exp(-12*yv(j)*yv(j)/pi)  
      ENDDO
      END 


      
