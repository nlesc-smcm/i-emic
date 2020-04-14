!**********************************************************************

      subroutine levitus_internal(filename,array,type)
      use m_global
      use m_lev
      implicit none
      integer k,klev,kk
      character(len=*),intent(in) :: filename
      character(len=*),intent(in):: type
      real, dimension(n,m,l),intent(out) :: array

      real, dimension(n,m) :: layer
      real dep
      integer :: nlayers

      nlayers=nlev

      write(f99,*) 'levitus temperature interpolation:'
      write(f99,*) 'model: level  depth   levitus: level depth   max'
      do k=l,1,-1
         dep = -z(k)*hdim
         do kk=1,nlayers
            if (depth(kk).le.dep) klev = kk
         enddo
         call levitus_interpol(filename,&
              &                            layer,-5.,50.,k,klev)
         if (type.eq.'TEMP') then
            layer = layer - t0
         else if (type.eq.'SALT') then
            layer=layer-s0
         else
            write(*,*) 'error unknown type please write TEMP or SALT'
         end if

         array(:,:,k) = layer
         !write(f99,999) k,dep,klev,depth(klev),tatmmax
      enddo
      !
      !write(f99,*) 'write levitus temperature field to: ',filename
      !call write_internal_forcing(filename,ftlev,34)

      end subroutine levitus_internal

      !**********************************************************************
      subroutine levitus_sst
        use m_global
        implicit none
        real tatmmax

        !      write(*,*) "ENTER LEVITUS_SST"

        !choose the followings

      write(*,*) '===========SSTforcing============================================'
      write(*,*) 'SST forcing is read in from file '//trim(sstfile)


      call levitus_interpol(locate_file(trim(sstfile)), tatm, &
                                                   -5.,50.,l,1)

      ! call levitus_interpol(locate_file('levitus/new/t00an1'),tatm,&
                                                    !-5.,50.,l,1)
      ! call levitus_interpol(locate_file('cesm/SST_38Ma.txt'),tatm,&
      !                                               -5.,50.,l,1)
      ! call levitus_interpol(locate_file('levitus/new/avtemp'),tatm,&
                                                    !-5.,50.,l,1)

      write(*,*) '===========SSTforcing============================================'

      open(34,FILE=rundir//'fort.34')
      call write_forcing('temp',tatm,34)
      close(34)
! for some reason, this statement doesn't work on Huygens:
!      tatm = tatm - t0
      tatm(1:n,1:m) = tatm(1:n,1:m) - t0
      tatmmax = maxval(tatm)
      write(f99,*) 'fit of levitus temp field done, tatmmax =',tatmmax+t0

      end
      !**********************************************************************


!**********************************************************************
      subroutine levitus_sal
        use m_global
        implicit none
        real emipmax

      write(*,*) '===========SSSforcing============================================'
      write(*,*) 'SSS forcing is read in from file '//trim(sssfile)

      call levitus_interpol(locate_file(trim(sssfile)), emip,&
                                                   20.,40.,l,1)

      ! call levitus_interpol(locate_file('levitus/new/s00an1'),emip,&
      !                                                     30.,40.,l,1)

      ! call levitus_interpol(locate_file('cesm/SSS_38Ma.txt'),emip,&
      !                                                     20.,40.,l,1)
      ! call levitus_interpol(locate_file('levitus/new/avsalt'),emip,&
                                                          !30.,40.,l,1)

      write(*,*) '===========SSSforcing============================================'

      open(33,FILE=rundir//'fort.33')
      call write_forcing('salt',emip,33)
      close(33)
      emip(1:n,1:m) = emip(1:n,1:m) - s0
      emipmax = maxval(emip)
      write(*,*) 'fit of levitus salt field done, emipmax =',emipmax+s0

      end


!**********************************************************************
      subroutine levitus_interpol(filename,forc,lolimit,uplimit,k,klev)
!
! This interpolation routine uses all Levitus points that fall within
! the model gridbox under consideration, rather than using only the
! closests data points.
!
      use m_global
      implicit none
! IMPORT/EXPORT
      character*(*) filename
      real     forc(n,m), lolimit, uplimit
      integer  k, klev
! LOCAL
      real     for
      integer  i, j, kl, nmis
      integer  weight
      real     xilow, xihigh, yjlow, yjhigh
      integer  iilow, iihigh, jjlow, jjhigh
      integer, parameter :: nx  = 360
      integer, parameter :: ny  = 180
      real,    parameter :: missing = -99.9999
      real     dat(0:nx, ny)


      open(1,file=filename,form='formatted',err=123)
      do kl=1,klev
        read(1,'(10f8.4)',err=123) ((dat(i,j),i=1,nx),j=1,ny)
      enddo
      close(1)
!      write(*,*) filename(38:)
!      open(unit=1,file=filename(38:))
!      do j=1,ny
!        do i=1,nx
!          write(1,*) dat(i,j)
!        end do
!      end do
!      close(1)
      dat(0,:) = dat(nx,:)

      do i=0,nx
      do j=1,ny
          if (dat(i,j).gt.(missing+10.0)) then
            if (dat(i,j).gt.uplimit) dat(i,j) = uplimit
            if (dat(i,j).lt.lolimit) dat(i,j) = lolimit
          endif
      enddo
      enddo
!w    where( dat < lolimit ) dat = lolimit
!w    where( dat > uplimit ) dat = uplimit
!
      forc   = missing/20.
      do j = 1, m
        do i = 1, n
          if (landm(i,j,k).eq.OCEAN) then
            nmis   = 0
            xilow  = 180.*(x(i)-0.5*dx)/pi
            xihigh = 180.*(x(i)+0.5*dx)/pi
            iilow  = ceiling(xilow)
            iihigh = floor(xihigh)
            if (iilow .lt.0)  iilow = 0
            if (iihigh.gt.nx) iihigh = nx
!
            yjlow  = 180.*(y(j)-0.5*dy)/pi
            yjhigh = 180.*(y(j)+0.5*dy)/pi
            jjlow  = ceiling(yjlow+90.5)
            jjhigh = floor(yjhigh+90.5)
!
 100        continue
            call interpol(iilow,iihigh,jjlow,jjhigh,dat,for,weight)
!
            if (weight.eq.0) then
               nmis = nmis + 1
               write(*,*)'levitus miss at :',k,i,j, nmis
               iilow  = iilow - 1
               jjlow  = jjlow - 1
               iihigh = iihigh + 1
               jjhigh = jjhigh + 1
               if (iilow .lt.0)  iilow  = 0
               if (jjlow .lt.1)  jjlow  = 1
               if (iihigh.gt.nx) iihigh = nx
               if (jjhigh.gt.ny) jjhigh = ny
               if (nmis.ge.10) stop 'definite levitus miss'
               goto 100
            endif
            forc(i,j) = for/weight

!$$$            write(f99,999) i,j,x(i)*180./pi,y(j)*180./pi,xilow,xihigh,
!$$$     >               yjlow,yjhigh,for, weight
          endif
        enddo
      enddo
!
! An additional smoothing can be done
!
!     call smooth(forc,k)
      return
      
123   call throw_error('failure')

      end
!******************************************************************
      SUBROUTINE interpol(iilow,iihigh,jjlow,jjhigh,dat,for,weight)
      use m_global
      implicit none
!     IMPORT/EXPORT
      integer, parameter :: nx  = 360
      integer, parameter :: ny  = 180
      real,    parameter :: missing = -99.9999
      real     dat(0:nx, ny)
      integer  iilow,iihigh,jjlow,jjhigh,weight
      real     for
!     LOCAL
      integer ii,jj

      weight = 0
      for  = 0.0
      do ii = iilow, iihigh
        do jj = jjlow, jjhigh
          if (dat(ii,jj).gt.(missing+10.0)) then
            for    = for + dat(ii,jj)
            weight = weight + 1
          endif
        enddo
      enddo

      end
!******************************************************************
      SUBROUTINE smooth(forc,k)
      use m_global
      implicit none

      integer i,j,k,iip,iim,jjp,jjm,weight
      real forc(n,m),forc1(n,m)

      forc1 = forc
      do j = 1, m
      do i = 1, n
        if (landm(i,j,k).eq.OCEAN) then
          iip = i+1
          iim = i-1
          jjp = j+1
          jjm = j-1
          if (iip.gt.n) iip = iip-n
          if (iim.lt.1) iim = iim+n
          if (jjp.gt.m) jjp = m
          if (jjm.lt.1) jjm = 1
          weight     = 1
          if (landm(iip,j,k).eq.OCEAN) then
             forc1(i,j) = forc1(i,j) + forc(iip,j)
             weight     = weight + 1
          endif
          if (landm(iim,j,k).eq.OCEAN) then
             forc1(i,j) = forc1(i,j) + forc(iim,j)
             weight     = weight + 1
          endif
          if (landm(i,jjp,k).eq.OCEAN) then
             forc1(i,j) = forc1(i,j) + forc(i,jjp)
             weight     = weight + 1
          endif
          if (landm(i,jjm,k).eq.OCEAN) then
             forc1(i,j) = forc1(i,j) + forc(i,jjm)
             weight     = weight + 1
          endif
          forc1(i,j) = forc1(i,j)/weight
        endif
      enddo
      enddo
      forc = forc1
      end
