      SUBROUTINE mstream(un,psimin,psimax,psisym,isalt,itemp,teq)
      implicit none
      include 'usr.com'
      real    un(ndim),psimin,psimax,psisym
      integer isalt,itemp 
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      integer i, j, k
      real    vs(0:m,l), psim(0:m,0:l), cs, dum, transc,fwf(n,m),teq
      real    salfun,hfun
      real    dum1, dum2

      call usol(un,u,v,w,p,t,s)

*     Calculation atmos equilibrium temperature
      if (la.ge.1) then
         dum1 = 0.0
         dum2 = 0.0
         do i=1,n
            do j=1,m
               cs = cos(y(j))
               dum1 = dum1 + t(i,j,l+la)*cs
               dum2 = dum2 + cs
            enddo
         enddo
          teq = dum1/dum2
       else
          teq = 0.0
       endif

*     Diagnostic atmospheric temperature 
      if (la.ge.1) then
         do j=1,m
            do i=1,n
               write(91,*) i,j,t(i,j,l),t(i,j,l+1)
            enddo  
         enddo
      endif

*     Calculation meridional overturning streamfunction
      vs = 0.0
      DO j=0,m
         DO k=1,l
            dum = 0.0
            DO i=1,n
               dum=0.5*(v(i-1,j,k)+v(i,j,k))*dx+dum
            ENDDO
            vs(j,k) = dum
         ENDDO
      ENDDO
      psim(:,0) = 0.0
      psim(0,:) = 0.0
      DO j=1,m
         cs = cos(yv(j))
         DO k=1,l
            psim(j,k)=-cs*vs(j,k)*dz*dfzT(k)  + psim(j,k-1)
         ENDDO
      ENDDO
      transc = r0dim*hdim*udim/1.0e+06 ! Sverdrups
      psimax = transc*maxval(psim)
      psimin = transc*minval(psim)
      psisym = psimin+ psimax
*
      if (isalt.eq.1)  then ! write freshwater flux -> fort.15 (flux.salt) 
          do j = 1,m
            do i = 1,n
              if (landm(i,j,l).eq.OCEAN) then
                fwf(i,j) = par(BIOT)*(salfun(x(i),y(j)) - s(i,j,l)/
     &                                        (par(COMB)*par(SALT)) )
              endif
            enddo
         enddo
         call write_forcing("flux.salt", fwf,15) 
      endif 
      if (itemp.eq.1)  then ! write heat flux -> fort.15 (flux.heat) 
          do j = 1,m
            do i = 1,n
              if (landm(i,j,l).eq.OCEAN) then
                fwf(i,j) = par(BIOT)*(hfun(x(i),y(j)) - t(i,j,l)/
     &                                      (par(COMB)*par(TEMP)) )
              endif
            enddo
         enddo
         call write_forcing("flux.heat", fwf,15) 
      endif 
      END
*
      SUBROUTINE bstream(un,psimin,psimax,psisym)
      implicit none
      include 'usr.com'
      real    un(ndim),psimin,psimax,psisym
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      integer i, j, k
      real    us(0:n,0:m), psib(0:n,0:m), cs, dum, transc

*     Calculation barotropic  streamfunction

      call usol(un,u,v,w,p,t,s)
      us = 0.0
      DO i=0,n
         DO j=0,m
            dum = 0.0
            DO k=1,l
                dum = u(i,j,k)*dz*dfzT(k)  + dum
            ENDDO
            us(i,j) = dum
         ENDDO
      ENDDO
      psib(:,0) = 0.0
      DO i=0,n
         DO j=1,m
            psib(i,j)=0.5*(us(i,j-1)+us(i,j))*dy + psib(i,j-1)
         ENDDO
      ENDDO
      transc = r0dim*hdim*udim/1.0e+06 ! Sverdrups
      psimax = transc*maxval(psib)
      psimin = transc*minval(psib)
      psisym = psimax + psimin
      END
* 
      SUBROUTINE mstream_basin(un,psimin,psimax,bas)
      implicit none
      include 'usr.com'
      real    un(ndim),psimin,psimax
      character*3 bas
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      integer i, j, k
      real    vs(0:m,l), psim(0:m,0:l), cs, dum, transc,fwf(n,m)
      real    salfun,hfun
      integer sland(0:n+1,0:m+1,0:l+la+1)

*     Calculation meridional overturning streamfunction

      call usol(un,u,v,w,p,t,s)
      call submask(bas,sland) 
      vs = 0.0
      DO j=0,m
         DO k=1,l
            dum = 0.0
            DO i=1,n
               dum=0.5*(v(i-1,j,k)+v(i,j,k))*sland(i,j,k)*dx+dum
            ENDDO
            vs(j,k) = dum
         ENDDO
      ENDDO
      psim(:,0) = 0.0
      psim(0,:) = 0.0
      DO j=1,m
         cs = cos(yv(j))
         DO k=1,l
            psim(j,k)=-cs*vs(j,k)*dz*dfzT(k)  + psim(j,k-1)
         ENDDO
      ENDDO
      transc = r0dim*hdim*udim/1.0e+06 ! Sverdrups
      psimax = transc*maxval(psim)
      psimin = transc*minval(psim)
      END 
******************************************************************
      SUBROUTINE submask(bas,sland)
      implicit none
      include 'usr.com'
      integer sland(0:n+1,0:m+1,0:l+la+1)
      character*3 bas
      integer i,j,k 
      open(unit=50,file='../_mkmask/mask.glo_'//bas) 
      do k = 0, l+la+1
         read(50,*)
         do j = m+1, 0, -1
            read(50,'(180i1)') (sland(i,j,k),i=0,n+1)
         enddo
      enddo
      close(50)
      END
******************************************************************
      SUBROUTINE monitor(itp,idum,ifile)
      implicit none
      include 'bag.com'
      integer  itp,ifile,isalt,itemp,idum
      real       teq
*
* local 
* 
      real psimmax, psimmin, psiamin,psiamax
      real psipmax, psipmin, psimsym
      real psibmax, psibmin, psibsym
      isalt = 0
      itemp = 0 
      call mstream(u,psimmin,psimmax,psimsym,isalt,itemp,teq)
      call bstream(u,psibmin,psibmax,psibsym) 
      write(ifile,'(2i3,9e12.4)') itp,icp,xl,psimmin,psimmax,psibmin,psibmax,teq
      end
******************************************************************
      SUBROUTINE write_ind
      implicit none
      write(7,*) 'itp par    parval    MOC-    MOC+    BSF-    BSF+  TEQ' 
      end 
     

