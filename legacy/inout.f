***********************************************************************
** FILE FORMAT VERSION 0 **********************************************
***********************************************************************

      subroutine write_data(ofile,lab)
      implicit none
      include 'bag.com'
      integer, intent(in)    :: ofile
      integer, intent(inout) :: lab
      integer  i, j, k, XX, row, j2, nskip, find_row2

      if ( lab.eq.0 ) call write_geometry(44)
*     call write_geometry(44)

      lab = lab + 1

      open(ofile)

      icp = mod(icp,1000)

      nskip= int((npar-1)/5 + 1) + 1 + nf + ndim*int((nf+1)/10+1)

      write(ofile,"('Version   0',8i4,2i12)") 
     >                            lab,icp,npar,nf,n,m,l,nun,ndim,nskip
      write(ofile,'(5(e18.10e3,X))') (par(i),i=1,npar)
      write(ofile,999) xl,xlp,det
      do j = 1, nf
         write(ofile,999) sig(j,1),sig(j,2)
      enddo
      ! write solution in "old" ordering
      do k = 1, (l+la)
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

 999  format(e18.10e3,9(X,e16.8e3))
      end

***********************************************************************
      subroutine append_data(ofile,lab)
      implicit none
      include 'bag.com'
      integer, intent(in)    :: ofile
      integer, intent(inout) :: lab
      integer  i, j, k, XX, row, j2, find_row2, nskip

      open(unit=ofile,form='formatted',position='append',
     >                             status='unknown')

      icp = mod(icp,1000)

      nskip= int((npar-1)/5 + 1) + 1 + nf + ndim*int((nf+1)/10+1)

      write(ofile,"('Version   0',8i4,2i12)")
     >                            lab,icp,npar,nf,n,m,l,nun,ndim,nskip
      write(ofile,'(5(e18.10e3,X))') (par(i),i=1,npar)
      write(ofile,999) xl,xlp,det
      do j = 1, nf
         write(ofile,999) sig(j,1),sig(j,2)
      enddo
      ! write solution in "old" ordering
      do k = 1, (l+la)
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

 999  format(e18.10e3,9(X,e16.8e3))
      end

***********************************************************************

      subroutine write_geometry_v0(gfile)
      implicit none
      include 'usr.com'
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
 999  format(5(e15.7e3,x))
      end

***********************************************************************

      subroutine write_geometry_v1(gfile)
      implicit none
      include 'usr.com'
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
 999  format(100(e15.7e3,x))
      end

***********************************************************************

      subroutine write_geometry(gfile)
      implicit none
      include 'usr.com'
      integer, intent(in) :: gfile
      integer i, j, k
      write(gfile,"('Version   2')")
      write(gfile,"(6i4)") n, m, l, la, nun, SLIP
      write(gfile,999) xmin, xmax, ymin, ymax, hdim
      write(gfile,999) ( x(i), i = 1, n)
      write(gfile,999) ( y(i), i = 1, m)
      write(gfile,999) ( z(i), i = 1, l)
      write(gfile,999) (xu(i), i = 0, n)
      write(gfile,999) (yv(i), i = 0, m)
      write(gfile,999) (zw(i), i = 0, l)
      do k = 0, l+la+1
      do j = 0, m+1
         write(gfile,'(100(i1,x))') (landm(i,j,k), i = 0, n+1)
      enddo
      enddo
 999  format(5(e15.7e3,x))
      end

***********************************************************************
      subroutine read_data(ifile, lab, icpo)
      implicit none
      include 'bag.com'
      integer, intent(in)  :: ifile, lab
      integer, intent(out) :: icpo
      character*7   ve 
      integer i, j, k, XX, row, j2, find_row2
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
         read(ifile,'(a7,9i4,2i12)') 
     >                ve,vn,irs,icpo,npa,nff,nn,mm,ll,nunn,ndi,nskip
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
               do k = 1, (l+la)
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


 999  format(e18.10e3,9(X,e16.8e3))
      end

***********************************************************************

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

********************************************************************************
      SUBROUTINE repnt(irs,icpo)
      implicit none
      include 'bag.com'
*     IMPORT/EXPORT
      integer  irs,icpo
*     LOCAL
      integer  i,j,k,XX,row,find_row2
      integer  LL,nskip,nf1,ndim1,lab,n11,m11,l11,nun11
      logical  eof4
*     FUNCTION
      real     g05caf,dum1,dum2
*
      rewind 4
 100  CONTINUE   
*         read(4,*,end=200) lab,icpo,nf1,ndim1,nskip,n11,m11,l11,nun11
         read(4,*,end=200) lab,icpo,nf1,ndim1,nskip
         IF (lab.EQ.irs) THEN
*           IF (ndim.NE.ndim1.or.n11.ne.n.or.m11.ne.m.or.l11.ne.l) THEN
           IF (ndim.NE.ndim1) THEN
             write(99,*) 'dimension of the initial point wrong'
             STOP 'in repnt: dimension of the initial point wrong'
           ENDIF
           IF (nf.NE.nf1) THEN
             write(99,*) 'Eigenv. No of initial point wrong'
*            STOP 'Eigenv. No of initial point wrong'
           ENDIF
           read(4,999,end=200) xl,xlp,det
           DO LL=1,nf
             read(4,999,end=200) sig(LL,1),sig(LL,2)
           ENDDO
*          DO LL=nf+1,12
*            read(4,999,end=200) dum1,dum2
*          ENDDO
           ! read data in "old" ordering
           do k = 1, (l+la)
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
*            read(4,999,end=200) u(i)
*            up(i)=0.
*            DO LL=5,nf
*              w(i,LL)=2*g05caf()-1.0
*            ENDDO
!           ENDDO
           read(4,999,end=200) (par(i),i=1,npar)
         ELSE
            call skip4(nskip,eof4)
            IF (eof4) GOTO 200
            GOTO 100
         END IF
      write(99,*) 'repnt done'
      RETURN
 200  STOP 'in repnt: error in reading fort.4'
*
 999  format(2X,e19.10e3,7e17.8e3)
*999  format(2X,8e16.8)
      END
******************************************************************
      SUBROUTINE skip4(nskip,eof4)
      implicit none
*     IMPORT/EXPORT
      integer  nskip
      logical  eof4
*     LOCAL
      integer  i
*
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
***********************************************************************
      subroutine write_forcing(filename,dat,fno)
      implicit none
      character*(*) filename
      include 'par.com'
      real dat(n, m)
      integer fno, i, j
      write(fno,*) n, m
      do j = 1, m
         write(fno,*) (dat(i,j),i=1,n)
      enddo
      end
***********************************************************************
      subroutine write_internal_forcing(filename,dat,fno)
      implicit none
      character*(*) filename
      include 'par.com'
      real dat(n, m, l)
      integer fno, i, j, k
c
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
c
      end

***********************************************************************
      subroutine read_internal_forcing(filename,dat,fno)
      implicit none
      character*(*) filename
      include 'par.com'
      real dat(n, m, l)
      integer fno, i, j, k, nn, mm, ll
c
      open(unit=fno,file=filename,form='formatted',status='unknown')
      read(fno,'(3i5)') nn, mm, ll
      if ((nn.ne.n) .or. (mm.ne.m) .or. (ll.ne.l)) then
        write(99,*) 'Levitus file ',filename,' not compatible'
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
c
      end
***********************************************************************
      subroutine read_forcing(dat,fno)
      implicit none
      include 'par.com'
      real dat(n,m)
      integer fno, i, j, n1 ,m1
      read(fno,*) n1, m1
      do j = 1, m1
         read(fno,*) (dat(i,j),i=1,n1)
      enddo
      end 


