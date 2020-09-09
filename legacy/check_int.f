********************************************************************
      SUBROUTINE outreg(un)
*     write solution
      implicit none
      include 'usr.com'
      include 'atm.com'
*     IMPORT/EXPORT
      real    un(ndim)
*     LOCAL
      integer i,j,k
      real    u(0:n,0:m,0:l+la+1),v(0:n,0:m,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la),p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+1+la),s(0:n+1,0:m+1,0:l+1+la)
      real    rho(0:n+1,0:m+1,0:l+1+la)
      real    cay,tay,c1,c2,dum1,dum2,h1,h2,xnes,lambda
      real    check,swl,swo
      real    hfun, salfun
*
      call usol(un,u,v,w,p,t,s)
      xnes = par(NLES)
      lambda = par(LAMB)
      rho    = lambda*s - t *( 1 + xnes*alpt1) - 
     +            xnes*t*t*alpt2+xnes*t*t*t*alpt3
      call forcing
C----------------------------------------------------------
C    integral check 1: advection of salt
      write(17,*)  'advection of salt'
      check= 0.0
      do i = 1, n
      do j = 1, m
      do k = 1, l
         if( landm(i,j,l) == OCEAN ) then
         check = check + 
     &      (u(i,j,k)+u(i,j-1,k))*(s(i+1,j,k)+s(i,j,k))/(4*dx)-
     &      (u(i-1,j,k)+u(i-1,j-1,k))*(s(i,j,k)+s(i-1,j,k))/(4*dx)+
     &      (v(i,j,k)+v(i-1,j,k))*(s(i,j+1,k)+s(i,j,k))*cos(yv(j))/(4*dy)-
     &      (v(i,j-1,k)+v(i-1,j-1,k))*(s(i,j,k)+s(i,j-1,k))*cos(yv(j-1))/(4*dy)+
     &      w(i,j,k)*(s(i,j,k+1)+s(i,j,k))*cos(y(j))/(2*dz*dfzW(k)) -
     &      w(i,j,k-1)*(s(i,j,k)+s(i,j,k-1))*cos(y(j))/(2*dz*dfzW(k-1))
          endif 
      enddo
      enddo
      enddo
      write(17,*) 'int_V \nabla(u S) d^3x =', check * dx*dy*dz
C    integral check 2: surface integral of freshwater flux
      write(17,*)  'surface salt flux'
      check = 0.0
      do i=1,n
      do j=1,m
         if( landm(i,j,l) == OCEAN ) then 
            if (SRES == 0) then
               check = fs(i,j) * cos(y(j)) + check
            else
               check = par(BIOT)*(par(SALT)*par(COMB)*salfun(x(i),y(j))  - s(i,j,l))
            endif
         endif    
      enddo
      enddo
      write(17,*) 'int_S F_S d^2x =', check * dx * dy
C    integral check 3: conservation of salt      
      write(17,*)  'salt conservation'
      check = 0.0
      do k = 1, l
      do i = 1, n
      do j = 1, m
         if( landm(i,j,k) == OCEAN ) then
            check = s(i,j,k) * cos(y(j)) * dfzT(k) + check
         endif
      enddo
      enddo
      enddo
      write(17,*) 'int_V  S d^3x =', check*dx*dy*dz
C    integral check 4: mixing terms (note: hor and vert both zero)
C    so not PH and PV necessary 
      write(17,*)  'salt diffuson'
      check = 0.0
      do k = 1, l
         h1 = 1./(dfzT(k)*dfzW(k))
         h2 = 1./(dfzT(k)*dfzW(k-1))
         do j = 1, m
            cay = cos(y(j))
            c1 = cos(yv(j))
            c2 = cos(yv(j-1))
            do i = 1, n
               if( landm(i,j,k) == OCEAN ) then
               check =  check + cos(y(j)) * dfzT(k) * (
     &(   s(i+1,j,k)+   s(i-1,j,k)-      2*s(i,j,k))/(dx*dx*cay*cay)+
     &(c1*s(i,j+1,k)+c2*s(i,j-1,k)-(c1+c2)*s(i,j,k))/(dy*dy*cay)    +
     &(h1*s(i,j,k+1)+h2*s(i,j,k-1)-(h1+h2)*s(i,j,k))/(dz*dz)        )
               endif 
            enddo
         enddo
      enddo
      write(17,*) 'int_V  \nabla^2 S d^3x = ', check *dx*dy*dz
C----------------------------------------------------------
C    integral check 5: advection of heat
      write(17,*)  'advection of heat'
      check = 0.0
      do i = 1, n
      do j = 1, m
      do k = 1, l
         if( landm(i,j,k) == OCEAN ) then
         check = check + 
     &      (u(i,j,k)+u(i,j-1,k))*(t(i+1,j,k)+t(i,j,k))/(4*dx)-
     &      (u(i-1,j,k)+u(i-1,j-1,k))*(t(i,j,k)+t(i-1,j,k))/(4*dx)+
     &      (v(i,j,k)+v(i-1,j,k))*(t(i,j+1,k)+t(i,j,k))*cos(yv(j))/(4*dy)-
     &      (v(i,j-1,k)+v(i-1,j-1,k))*(t(i,j,k)+t(i,j-1,k))*cos(yv(j-1))/(4*dy)+
     &      w(i,j,k)*(t(i,j,k+1)+t(i,j,k))*cos(y(j))/(2*dz*dfzW(k)) -
     &      w(i,j,k-1)*(t(i,j,k)+t(i,j,k-1))*cos(y(j))/(2*dz*dfzW(k-1))
          endif 
      enddo
      enddo
      enddo
      write(17,*) 'int_V \nabla(u T) d^3x =', check * dx*dy*dz
C    integral check 6: surface integral of heat flux (uncoupled case) 
      IF (la == 0) THEN
      write(17,*)  'surface heat flux, uncoupled'
      check = 0.0
      do i=1,n
      do j=1,m
         if( landm(i,j,l) == OCEAN ) then 
            if (TRES == 0) then
               check = ft(i,j) * cos(y(j)) + check
            else
               check = par(BIOT)*(par(TEMP)*par(COMB)*hfun(x(i),y(j))  - t(i,j,l))
            endif
         endif    
      enddo
      enddo
      write(17,*) 'int_S Q_T d^2x =', check * dx * dy
C    integral check 7: mixing terms (note: hor and vert both zero)
C    so not PH and PV necessary
      write(17,*)  'heat diffuson'
      check = 0.0
      do k = 1, l
         h1 = 1./(dfzT(k)*dfzW(k))
         h2 = 1./(dfzT(k)*dfzW(k-1))
         do j = 1, m
            cay = cos(y(j))
            c1 = cos(yv(j))
            c2 = cos(yv(j-1))
            do i = 1, n
                if( landm(i,j,l) == OCEAN ) then 
               check =  check + cos(y(j)) * dfzT(k) * (
     &(   t(i+1,j,k)+   t(i-1,j,k)-      2*t(i,j,k))/(dx*dx*cay*cay)+
     &(c1*t(i,j+1,k)+c2*t(i,j-1,k)-(c1+c2)*t(i,j,k))/(dy*dy*cay)    +
     &(h1*t(i,j,k+1)+h2*t(i,j,k-1)-(h1+h2)*t(i,j,k))/(dz*dz)        )
                endif 
            enddo
         enddo
      enddo
      write(17,*) 'int_V \nabla^2 T d^3x = ', check *dx*dy*dz
      ENDIF
C----------------------------------------------------------      
C    integral check 8: mass balance
      write(17,*)  'volume conservation'
      check = 0.0
      do j = 1, m
         cay = cos(y(j))
         c1 = cos(yv(j))
         c2 = cos(yv(j-1))
         do k = 1, l
         do i = 1, n
            if( landm(i,j,l) == OCEAN ) then 
            check = check +dfzT(k)*( (u(i,j,k) +u(i,j-1,k)  - u(i-1,j,k)  - u(i-1,j-1,k)  )/(2*dx)+
     +                     ((v(i,j,k)+v(i-1,j,k))*c1- (v(i,j-1,k)+v(i-1,j-1,k))*c2)/(2*dy)+
     +                     (w(i,j,k)   - w(i,j,k-1))*cay/(dz*dfzT(k)) ) 
             endif 
         enddo
         enddo
      enddo
      write(17,*) 'int_V \nabla u d^3x  =', check * dx*dy*dz
*    integral check 9:  geometric conservation
      write(17,*)  'geometry check'
      check = 0.0
      do k = 1, l
         check = check + dfzT(k)
      enddo
      write(17,*) 'int_z dfz dz - 1 .... =', check * dz - 1.
C---------------------------------------------------------- 
      IF (la == 1) THEN
C    integral check 10: land energy balance
      write(17,*)  'surface land-atm heat flux'
      check = 0.0
      swl = 0.0
      do i=1,n
         do j=1,m
              IF (landm(i,j,l) == LAND) THEN
                   check = check + (t(i,j,l) - t(i,j,l+1))*cos(y(j)) 
                   swl = par(COMB)*par(SUNP) * suno(j)*cos(y(j)) /Ooa + swl
              ENDIF
         enddo
      enddo
      write(17,*) 'int_S_l  (t-tatm)-sw d^2x =', (check-swl)*dx*dy     
C    integral check 11: ocean energy balance
      write(17,*)  'surface ocean-atm heat flux'
      check = 0.0
      swo = 0.0
      do i=1,n
         do j=1,m
              IF (landm(i,j,l) == OCEAN) THEN
                   check = check + (t(i,j,l) - t(i,j,l+1))*cos(y(j)) 
                   swo = par(COMB)*par(SUNP)*suno(j)*cos(y(j)) /Ooa + swo
              ENDIF
         enddo
      enddo
      write(17,*) 'int_S_o (t-tatm)-sw d^2x=', (check-swo)*dx*dy                    
C    integral check 12: global energy balance
      write(17,*)  'surface energy balance'
      check = 0.0
      swo = 0.0
      do i=1,n
         do j=1,m
          check = check + t(i,j,l+1)*cos(y(j)) 
          swo = par(COMB)* par(SUNP)*(suna(j)+suno(j)/Ooa - amua)*cos(y(j))/bmua + swo
         enddo
      enddo
      write(17,*) 'int_S tatm -sw d^2x =', (check-swo)*dx*dy    
C---------------------------------------------------------- 
      ENDIF                  
      END
