!**********************************************************************
!* FILE FORMAT VERSION 0 **********************************************
!**********************************************************************

! NOTE: the subroutines in this file have been decoupled from the rest of THCM
!       (only modules m_par and m_global are used, and m_global is NOT used by anyone
!       else).
!       Some data in m_global is simply set to zero because LOCA
!       takes care of the continuation now. Any data you want to read/write must
!       be passed in as subroutine argument because the I/O may be for the
!       global problem whereas the computations are performed on a subdomain.


!     output function. uloca is expected to be the global solution
!     vector, sorted in 'interleaved' (loca) format. n/m/l
!     specify the grid size and may be different from m,n,l used
!     for the computations.
!     Be careful to call this only by the root process, otherwise
!     files may get severly messed up or data may be lost.
      subroutine write_data(uloca,ofile,lab)

      use m_par
      use m_global

      implicit none

      integer, intent(in)    :: ofile
      integer, intent(inout) :: lab
      integer  i, j, k, XX, row, nskip

      real, dimension(n*m*l*nun)  :: uloca

      if (ndim==0) then
        write(*,*) "WARNING: in write_data: you probably forgot to allocate_global!"
        write(*,*) "(",__FILE__,", line ",__LINE__,")"
      end if

!      write(*,*) "store solution in '"//rundir//"fort.44'"
!      write(*,*) "WARNING: landm array is assumed to be empty!"
!      write(*,*) "(",__FILE__,", line ",__LINE__,")"


!      if ( lab.eq.0 ) call write_geometry(44)
      open(44,file=rundir//'fort.44')
      call write_geometry(44)
      close(44)

      if (ofile.eq.0) then
         return
      endif

      u   = uloca(1:ndim)
      up  = 0
      w   = 0
      xl  = 0
      xlp = 0
      det = 0
      sig = 0

      lab = lab + 1

      open(ofile,file=rundir//'fort.3')

      write(*,*) n,m,l,0

      icp = mod(icp,1000)

      nskip= int((npar-1)/5 + 1) + 1 + nf + ndim*int((nf+1)/10+1)

      write(ofile,"('Version   0',8i4,2i12)") &
                                 lab,icp,npar,nf,n,m,l,nun,ndim,nskip
      write(ofile,'(5(e18.10e3,1X))') (par(i),i=1,npar)
      write(ofile,999) xl,xlp,det
      do j = 1, nf
         write(ofile,999) sig(j,1),sig(j,2)
      enddo
      ! write solution in "old" ordering
      do k = 1, l
         do j = 1, m
            do i = 1, n
               do XX = 1,nun
                  row = find_row2(i,j,k,XX)
                  write(ofile,999) u(row)!,up(row),(w(row,j2),j2=1,nf)
               enddo
            enddo
         enddo
      enddo

      close(ofile)

 999  format(e18.10e3,9(1X,e16.8e3))
      end

!**********************************************************************
      subroutine append_data(ofile,lab)
      use m_par
      use m_global
      implicit none
      integer, intent(in)    :: ofile
      integer, intent(inout) :: lab
      integer  i, j, k, XX, row, j2, nskip

      open(unit=ofile,form='formatted',position='append',&
                                  status='unknown')

      icp = mod(icp,1000)

      nskip= int((npar-1)/5 + 1) + 1 + nf + ndim*int((nf+1)/10+1)

      write(ofile,"('Version   0',8i4,2i12)")&
                                 lab,icp,npar,nf,n,m,l,nun,ndim,nskip
      write(ofile,'(5(e18.10e3,1X))') (par(i),i=1,npar)
      write(ofile,999) xl,xlp,det
      do j = 1, nf
         write(ofile,999) sig(j,1),sig(j,2)
      enddo
      ! write solution in "old" ordering
      do k = 1, l
         do j = 1, m
            do i = 1, n
               do XX = 1,nun
                  row = find_row2(i,j,k,XX)
                  write(ofile,999) u(row),up(row),(w(row,j2),j2=1,nf)
               enddo
            enddo
         enddo
      enddo

      close(ofile)

 999  format(e18.10e3,9(1X,e16.8e3))
      end

!**********************************************************************

      subroutine write_geometry_v0(gfile)
      use m_par
      use m_global
      implicit none
      integer, intent(in) :: gfile
      integer i, j


      write(gfile,"('Version   0',5i4)") n, m, l, nun, SLIP
      write(gfile,999) xmin, xmax, ymin, ymax, hdim
      write(gfile,999) ( x(i), i = 1, n)
      write(gfile,999) ( y(i), i = 1, m)
      write(gfile,999) ( z(i), i = 1, l)
      write(gfile,999) (xu(i), i = 0, n)
      write(gfile,999) (yv(i), i = 0, m)
      write(gfile,999) (zw(i), i = 0, l)
      do j = 1, m
         write(gfile,'(100(i1,x))') (landm(i,j,l), i = 1, n)
      enddo
 999  format(5(e15.7e3,1x))
      end

!**********************************************************************

      subroutine write_geometry_v1(gfile)
      use m_par
      use m_global
      implicit none
      integer, intent(in) :: gfile
      integer i, j, k
      write(gfile,"('Version   1',5i4)") n, m, l, nun, SLIP
      write(gfile,999) xmin, xmax, ymin, ymax, hdim
      write(gfile,999) ( x(i), i = 1, n)
      write(gfile,999) ( y(i), i = 1, m)
      write(gfile,999) ( z(i), i = 1, l)
      write(gfile,999) (xu(i), i = 0, n)
      write(gfile,999) (yv(i), i = 0, m)
      write(gfile,999) (zw(i), i = 0, l)
      do k = 0, l+1
      do j = 0, m+1
         write(gfile,'(100(i1,x))') (landm(i,j,k), i = 0, n+1)
      enddo
      enddo
 999  format(100(e15.7e3,1x))
      end

!**********************************************************************

      subroutine write_geometry(gfile)
      use m_par
      use m_global
      implicit none
      integer, intent(in) :: gfile
      integer i, j, k

      write(gfile,"('Version   2')")
      write(gfile,"(6i4)") n, m, l, 0, nun, SLIP
      write(gfile,999) xmin, xmax, ymin, ymax, hdim
      write(gfile,999) ( x(i), i = 1, n)
      write(gfile,999) ( y(i), i = 1, m)
      write(gfile,999) ( z(i), i = 1, l)
      write(gfile,999) (xu(i), i = 0, n)
      write(gfile,999) (yv(i), i = 0, m)
      write(gfile,999) (zw(i), i = 0, l)
      do k = 0, l+1
      do j = 0, m+1
         write(gfile,'(100(i1,x))') (landm(i,j,k), i = 0, n+1)
      enddo
      enddo
 999  format(5(e15.7e3,1x))
      end

!**********************************************************************
      subroutine read_data(ifile, lab, icpo)
      use m_par
      use m_global
      implicit none
      integer, intent(in)  :: ifile, lab
      integer, intent(out) :: icpo
      character*7   ve
      integer i, j, k, XX, row, j2
      integer vn, irs, npa, nff, nn, mm, ll, nunn, ndi, nskip

      open(ifile)

      read(ifile,*) ve
      if (ve .ne. 'Version') then
         if(ifile.ne.4) stop 'old files should be in fort.4'
         call repnt(lab, icpo)
         return
      endif
      rewind ifile
      do
         read(ifile,'(a7,9i4,2i12)') &
                     ve,vn,irs,icpo,npa,nff,nn,mm,ll,nunn,ndi,nskip
         if (lab.eq.irs) then
            if (npar.ne.npa) stop 'in read_data: wrong no of parametrs'
            read(ifile,'(5(e18.10e3,X))') (par(i),i=1,npar)
            read(ifile,999) xl,xlp,det
            if (irs.ne.0) then
               if(ndim.ne.ndi) stop 'in read_data: wrong dimension'
               if(nf.ne.nff) stop 'in read_data: wrong no of eigenvals'
               do j = 1, nf
                  read(ifile,999) sig(j,1),sig(j,2)
               enddo
               ! read data in "old" ordering
               do k = 1, l
                  do j = 1, m
                     do i = 1, n
                        do XX = 1, nun
                           row = find_row2(i,j,k,XX)
                           read(ifile,999) u(row),up(row),(w(row,j2),j2=1,nf)
                        enddo
                     enddo
                  enddo
               enddo
            else
               u = 0.
               up = 0.
               w  = 0.
            endif
            close(ifile)
            return
         endif
         call skip(ifile,nskip)
      enddo


 999  format(e18.10e3,9(1X,e16.8e3))
      end

!**********************************************************************

      subroutine skip(ifile,nskip)
      implicit none
      integer, intent(in) :: ifile, nskip
      integer  i

      do i = 1, nskip
         read(ifile,'(x)', end = 100)
      enddo
      return

 100  stop 'in skip: unexpected end of file, check the label'
      end

!*******************************************************************************
      SUBROUTINE repnt(irs,icpo)
      use m_par
      use m_global
      implicit none
!     IMPORT/EXPORT
      integer  irs,icpo
!     LOCAL
      integer  i,j,k,XX,row
      integer  LL,nskip,nf1,ndim1,lab
      logical  eof4

      rewind 4
 100  CONTINUE
!         read(4,*,end=200) lab,icpo,nf1,ndim1,nskip,n11,m11,l11,nun11
         read(4,*,end=200) lab,icpo,nf1,ndim1,nskip
         IF (lab.EQ.irs) THEN
!           IF (ndim.NE.ndim1.or.n11.ne.n.or.m11.ne.m.or.l11.ne.l) THEN
           IF (ndim.NE.ndim1) THEN
             write(f99,*) 'dimension of the initial point wrong'
             STOP 'in repnt: dimension of the initial point wrong'
           ENDIF
           IF (nf.NE.nf1) THEN
             write(f99,*) 'Eigenv. No of initial point wrong'
!            STOP 'Eigenv. No of initial point wrong'
           ENDIF
           read(4,999,end=200) xl,xlp,det
           DO LL=1,nf
             read(4,999,end=200) sig(LL,1),sig(LL,2)
           ENDDO
!          DO LL=nf+1,12
!            read(4,999,end=200) dum1,dum2
!          ENDDO
           ! read data in "old" ordering
           do k = 1, l
              do j = 1, m
                 do i = 1, n
                    do XX = 1, nun
                       row = find_row2(i,j,k,XX)
                       read(4,999,end=200) u(row),up(row),(w(row,LL),LL=1,nf)
                    enddo
                 enddo
              enddo
           enddo
!           DO i=1,ndim
!             read(4,999,end=200) u(i),up(i),(w(i,LL),LL=1,nf)
!            read(4,999,end=200) u(i)
!            up(i)=0.
!            DO LL=5,nf
!              w(i,LL)=2*g05caf()-1.0
!            ENDDO
!           ENDDO
           read(4,999,end=200) (par(i),i=1,npar)
         ELSE
            call skip4(nskip,eof4)
            IF (eof4) GOTO 200
            GOTO 100
         END IF
      write(f99,*) 'repnt done'
      RETURN
 200  STOP 'in repnt: error in reading fort.4'
!
 999  format(2X,e19.10e3,7e17.8e3)
!999  format(2X,8e16.8)
      END
!*****************************************************************
      SUBROUTINE skip4(nskip,eof4)
      implicit none
!     IMPORT/EXPORT
      integer  nskip
      logical  eof4
!     LOCAL
      integer  i
!
      eof4=.false.
      DO i=1,nskip
         read(4,999,end=100)
      ENDDO
 999  format(1x)
      return
 100  CONTINUE
      eof4=.true.
      return
      END
!**********************************************************************
      subroutine write_forcing(filename,dat,fno)
      use m_par
      use m_global
      implicit none
      character*(*) filename
      real dat(n, m)
      integer fno, i, j
      write(fno,*) n, m
      open(unit=fno,file=filename)
      do j = 1, m
         write(fno,*) (dat(i,j),i=1,n)
      enddo
      close(fno)
      end
!**********************************************************************
      subroutine write_internal_forcing(filename,dat,fno)
      use m_par
      use m_global
      implicit none
      character*(*) filename
      real dat(n, m, l)
      integer fno, i, j, k
!
      open(unit=fno,file=filename,form='formatted',status='unknown')
      write(fno,'(3i5,a40)') n, m, l, filename
      do k = 1, l
      do j = 1, m
      do i = 1, n
         write(fno,'(e16.8e3)') dat(i,j,k)
      enddo
      enddo
      enddo
      write(6,*) filename,' succesfully written'
!
      end

!**********************************************************************
      subroutine read_internal_forcing(filename,dat,fno)
      use m_par
      use m_global
      implicit none
      character*(*) filename
      real dat(n, m, l)
      integer fno, i, j, k, nn, mm, ll
!
      open(unit=fno,file=filename,form='formatted',status='unknown')
      read(fno,'(3i5)') nn, mm, ll
      if ((nn.ne.n) .or. (mm.ne.m) .or. (ll.ne.l)) then
        write(f99,*) 'Levitus file ',filename,' not compatible'
        stop
      endif
      do k = 1, l
      do j = 1, m
      do i = 1, n
         read(fno,'(e16.8e3)') dat(i,j,k)
      enddo
      enddo
      enddo
      write(6,*) filename,' succesfully read'
!
      end
!**********************************************************************
      subroutine read_forcing(dat,fno)
      use m_par
      use m_global
      implicit none
      real dat(n,m)
      integer fno, i, j, n1 ,m1
      read(fno,*) n1, m1
      do j = 1, m1
         read(fno,*) (dat(i,j),i=1,n1)
      enddo
      end


!! this subroutine can be called from C++ to extract the global grid data from m_global
      subroutine get_grid_data(nn,mm,ll,xx,yy,zz)

      use m_global
      use m_par

      implicit none

      integer nn,mm,ll
      real xx(nn)
      real yy(mm)
      real zz(ll)

      integer i

      if (nn .ne. n .or. mm .ne. m .or. ll .ne. l) then
        write(*,*) 'subroutine get_grid_data was called with mismatched array dimensions!'
        write(*,*) '(inout.F, not returning grid data correctly!)'
        return;
      endif

      do i=1,n
        xx(i) = x(i)
      enddo
      do i=1,m
        yy(i) = y(i)
      enddo
      do i=1,l
        zz(i) = z(i)
      enddo

      end subroutine get_grid_data
