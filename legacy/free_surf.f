      subroutine frs_matrix(un,sig)
      USE m_mat
      implicit none
      include 'usr.com'

      real    un(ndim), sig
      real    u(0:n  ,0:m,0:l+la+1), v(0:n, 0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      real    fs1(n,m,np),fs2(n,m,np),fs3(n,m,np),fs4(n,m,np)
      real    eps_f
      
      eps_f=gdim*hdim/(2*omegadim*r0dim*udim)
      
      call usol(un,u,v,w,p,t,s)

      call frs_nlin(1,fs1,u,v,p)
      call frs_nlin(2,fs2,u,v,p)
      call frs_nlin(3,fs3,u,v,p)
      call frs_nlin(4,fs4,u,v,p)
      Al(:,:,l,:,WW,: ) = 0.0
      Al(:,:,l,5,WW,WW) =-1.0
      Al(:,:,l,:,WW,UU) = Al(:,:,l,:,WW,UU)+fs2*eps_f
      Al(:,:,l,:,WW,VV) = Al(:,:,l,:,WW,VV)+fs4*eps_f
      Al(:,:,l,:,WW,PP) = Al(:,:,l,:,WW,PP)+(fs1+fs3)*eps_f
c      Al(:,:,l,5,WW,PP) = Al(:,:,l,5,WW,PP)+sig*eps_f
      
      call frs_boundaries
      
      end
* ---------------------------------------------------------------------------- *
      subroutine frs_rhs(un)
      ! 1: matrix
      ! 2: rhs
      USE m_mat
      implicit none
      include 'usr.com'

      real    un(ndim)
      real    u(0:n  ,0:m,0:l+la+1), v(0:n, 0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)
      real    fs2(n,m,np),fs4(n,m,np)
      real    eps_f
      
      eps_f=gdim*hdim/(2.0*omegadim*r0dim*udim)

      call usol(un,u,v,w,p,t,s)

      call frs_nlin(2,fs2,u,v,p)
      call frs_nlin(4,fs4,u,v,p)
      Al(:,:,l,:,WW,: ) = 0.0
      Al(:,:,l,5,WW,WW) =-1.0
      Al(:,:,l,:,WW,UU) = Al(:,:,l,:,WW,UU)+fs2*eps_f
      Al(:,:,l,:,WW,VV) = Al(:,:,l,:,WW,VV)+fs4*eps_f

      call frs_boundaries
      call frs_forcing

      end
* ---------------------------------------------------------------------------- *
      subroutine frs_nlin(type,atom,u,v,p)
      implicit none
      include 'usr.com'
      ! 1: Uhx
      ! 2: uHx
      ! 3: Vhx
      ! 4: vHx
      
      integer type, i,j
      real atom(n,m,np)
      real u(0:n,0:m,0:l+la+1),v(0:n,0:m,0:l+la+1),p(0:n+1,0:m+1,0:l+la)
      real cos2i(m), udum, vdum, pdum
      
      atom=0.0
      select case(type)
        case(1) !Uhx
	  cos2i=8.0*cos(y)*dx
	  do j=1,m
	    do i=1,n
	      udum=u(i-1,j-1,l)+u(i-1,j,l)+u(i,j-1,l)+u(i,j,l)
	      atom(i,j, 2)=-1.5*udum/cos2i(j)
	      atom(i,j, 8)= 1.5*udum/cos2i(j)
	      atom(i,j,11)= 0.5*udum/cos2i(j)
	      atom(i,j,17)=-0.5*udum/cos2i(j)
	    enddo
	  enddo
	case(2) !uHx
	  cos2i=8.0*cos(y)*dx
	  do j=1,m
	    do i=1,n
	      pdum=1.5*p(i+1,j,l)-0.5*p(i+1,j,l-1)
     +	          -1.5*p(i-1,j,l)+0.5*p(i-1,j,l-1)
	      atom(i,j,1)= pdum/cos2i(j)
	      atom(i,j,2)= pdum/cos2i(j)
	      atom(i,j,4)= pdum/cos2i(j)
	      atom(i,j,5)= pdum/cos2i(j)
	    enddo
	  enddo
	case(3) ! Vhy
	  do j=1,m
	    do i=1,n
	      vdum=u(i-1,j-1,l)+u(i-1,j,l)+u(i,j-1,l)+u(i,j,l)
	      atom(i,j, 4)=-1.5*vdum/(8.0*dy)
	      atom(i,j, 6)= 1.5*vdum/(8.0*dy)
	      atom(i,j,13)= 0.5*vdum/(8.0*dy)
	      atom(i,j,15)=-0.5*vdum/(8.0*dy)
	    enddo
	  enddo
	case(4) ! vHy
	  do j=1,m
	    do i=1,n
	      pdum=1.5*p(i,j+1,l)-0.5*p(i,j+1,l-1)
     +	          -1.5*p(i,j-1,l)+0.5*p(i,j-1,l-1)
	      atom(i,j,1)= pdum/(8.0*dy)
	      atom(i,j,2)= pdum/(8.0*dy)
	      atom(i,j,4)= pdum/(8.0*dy)
	      atom(i,j,5)= pdum/(8.0*dy)
	    enddo
	  enddo
      end select      
      end
* ---------------------------------------------------------------------------- *
      subroutine frs_boundaries
      USE m_mat
      implicit none
      include 'usr.com'
      
      integer i,j,find_row2,row
      
      
      do i = 1, n
        do j = 1, m
          if (landm(i,j,l).eq.OCEAN) then
            if (landm(i-1,j-1,l  ).eq.LAND) then ! 1
              Al(i,j,l, 1,WW,UU) = 0.0
              Al(i,j,l, 1,WW,VV) = 0.0
            endif
            if (landm(i-1,j  ,l  ).eq.LAND) then ! 2
              Al(i,j,l, 2,WW,: ) = 0.0
              Al(i,j,l, 1,WW,UU) = 0.0
              Al(i,j,l, 1,WW,VV) = 0.0
            endif
            if (landm(i-1,j+1,l  ).eq.LAND) then ! 3
              Al(i,j,l, 2,WW,UU) = 0.0
              Al(i,j,l, 2,WW,VV) = 0.0
            endif
            if (landm(i  ,j-1,l  ).eq.LAND) then ! 4 
              Al(i,j,l, 4,WW,: ) = 0.0
              Al(i,j,l, 1,WW,UU) = 0.0
              Al(i,j,l, 1,WW,VV) = 0.0        
            endif
            if (landm(i  ,j+1,l  ).eq.LAND) then ! 6 
              Al(i,j,l, 6,WW,: ) = 0.0
            endif
            if (landm(i+1,j-1,l  ).eq.LAND) then ! 7
              Al(i,j,l, 4,WW,: ) = 0.0
            endif
            if (landm(i+1,j  ,l  ).eq.LAND) then ! 8 
              Al(i,j,l, 8,WW,: ) = 0.0
            endif
            if (landm(i+1,j+1,l  ).eq.LAND) then	! 9
              Al(i,j,l, 5,WW,UU) = 0.0
              Al(i,j,l, 5,WW,VV) = 0.0
            endif
          endif
          if ((landm(i,j,l).ne.OCEAN).and.(landm(i,j,l).ne.ATMOS)) then
            Al(i,j,l,:,WW,:) = 0.0
            row=find_row2(i,j,l,WW)
            Frc(row) = 0.0
            Al(i,j,l,5,WW,WW) = 1.0
          endif
        enddo
      enddo
      
      end
* ---------------------------------------------------------------------------- *
      subroutine frs_forcing
      USE m_mat
      implicit none
      include 'usr.com'
      
      integer i,j,find_row2,row
      
      do i=1,n
        do j=1,m
	  row=find_row2(i,j,l,WW)
	  Frc(row)=0.0
        enddo
      enddo	   
      
      end
* ---------------------------------------------------------------------------- *
      subroutine frs_test
      USE m_mat
      implicit none
      include 'usr.com'
      
      real un(ndim), rh(ndim), g05caf
      integer i,j,k,row, find_row2
      integer i2,j2,k2,e2
      integer i3,j3,k3,e3
      
      do i=1,ndim
!        un(i) = 1.0+0.1*g05caf()
      enddo

      call rhs(un,rh)
      call matrix(un,0.1)

      i=find_row2(8,8,16,WW)
      do j=bega(i),bega(i+1)-1
        call findex(jcoa(j),i2,j2,k2,e2)
        write (20,'(4i5,es16.8)') i2, j2, k2, e2, coa(j)
      enddo
      
      stop

      do i=1,ndim
	call findex(i,i3,j3,k3,e3)
          do k=bega(i),bega(i+1)-1
	    call findex(jcoa(k),i2,j2,k2,e2)
	    if (k2.eq.l) then
	    if (e2.eq.3)
     +	    write (20+e2,'(10i8,es16.8)') i,i3,j3,k3,e3,jcoa(k),i2,j2,k2,e2,coa(k)
	    if (e2.eq.4)
     +	    write (20+e2,'(10i8,es16.8)') i,i3,j3,k3,e3,jcoa(k),i2,j2,k2,e2,coa(k)
            endif
	  enddo  
      enddo
      stop
      
      end
* ---------------------------------------------------------------------------- *
      subroutine frs_test2
      USE m_mat
      implicit none
      include 'usr.com'
      
      integer i,j,k,k2,e,row, find_row2
      
      do i=1,ndim
          do k=bega(i),bega(i+1)-1
	    write (20,'(2i8,2es16.8)') i,jcoa(k),coa(k)
	  enddo  
      enddo
      stop      
      end
* ---------------------------------------------------------------------------- *
      subroutine write_frs(un)
      USE m_mat
      implicit none
      include 'usr.com'
      
      real un(ndim)
      integer i,j, find_row2, row
      
      do j=1,m
        do i=1,n
	  row=find_row2(i,j,l,WW)
	  write (20,*) un(row)
	enddo
      enddo

      end
      
      
      
      
      
      
      
      
      
      
