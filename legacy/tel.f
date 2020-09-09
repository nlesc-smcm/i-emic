      program tel
      integer n,m,l
      parameter(n=96,m=38,l=13)
      integer i,j,k,irow,icol,iu,itel
      write(6,*) 'give irow'
      read(5,*) irow
      write(6,*) 'give icol'
      read(5,*) icol
      do i=1,n
	 do j=1,m
	   do k=1,l
	      do iu=1,6
	      itel = 6*((k-1)*n*m+n*(j-1)+i-1) +iu 
	      if (itel.eq.irow) then
		 write(6,*) i,j,k,iu,itel
              endif 
	      if (itel.eq.icol) then
		 write(6,*) i,j,k,iu,itel
              endif 
              enddo
            enddo
           enddo
      enddo
      end
