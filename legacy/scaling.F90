      MODULE m_scaling
!
      ! scaling vectors
      REAL, DIMENSION(:), ALLOCATABLE :: rs, cs, rsa, csa
!
      CONTAINS
!**************************************************************
      subroutine rcscale_a(rl, rs, cs, rsa, csa)
      USE m_mat
      implicit none
      include 'usr.com'
!      include 'mat.com'
      real RL(ndim)
      real rs(nun), cs(nun), rsa(nun), csa(nun)
      integer i,j,k,ii,jj,kk,ee,ndimo,ix,iy,iz,jx,jy,jz

!     row and column scaling simultaneously
      do i = 1, ndim
         call findex(i,ix,iy,iz,ii)
         if (landm(ix,iy,iz) == OCEAN) then
            do j = begA(i), begA(i+1)-1
               call findex(jcoA(j),jx,jy,jz,jj)
               if (landm(jx,jy,jz) == OCEAN) then
                  coA(j) = rs(ii)*coA(j)*cs(jj)
               elseif (landm(jx,jy,jz) == ATMOS) then
                  coA(j) = rs(ii)*coA(j)*csa(jj)
	       else
                  coA(j) = rs(ii)*coA(j)
	       end if
            end do
            rl(i) = rs(ii)*rl(i)
         elseif (landm(ix,iy,iz) == ATMOS) then ! l < iz <= l+la
            do j = begA(i), begA(i+1)-1
               call findex(jcoA(j),jx,jy,jz,jj)
               if (landm(jx,jy,jz) == OCEAN) then
                  coA(j) = rsa(ii)*coA(j)*cs(jj)
	       elseif (landm(jx,jy,jz) == ATMOS) then
                  coA(j) = rsa(ii)*coA(j)*csa(jj)
               else
                  coA(j) = rsa(ii)*coA(j)
               end if
            end do
            rl(i) = rsa(ii)*rl(i)
         else 
            do j = begA(i), begA(i+1)-1
               call findex(jcoA(j),jx,jy,jz,jj)
               if (landm(jx,jy,jz) == OCEAN) then
                  coA(j) = coA(j)*cs(jj)
               elseif (landm(jx,jy,jz) ==ATMOS) then
                  coA(j) = coA(j)*csa(jj)
               end if
            end do
         end if
      end do

!      ndimo = nun*n*m*l
!     row scaling
!      do i = 0, ndimo-1, nun
!         call findex(i+1,ix,jy,kz,en)
!         if(landm(ix,jy,kz) == OCEAN) then
!            do ii = 1, nun 
!               do k = begA(i+ii), begA(i+ii+1)-1
!                  coA(k) = rs(ii) * coA(k)
!               enddo
!               rl(i+ii) = rs(ii) * rl(i+ii)
!            enddo
!         endif
!      enddo
!      do i = ndimo+1, ndim
!         ee = mod(i-1,nun) + 1
!         do k = begA(i), begA(i+1)-1
!            coA(k) = rsa(ee) * coA(k)
!         enddo
!         rl(i) = rsa(ee) * rl(i)
!      enddo
!
!     column scaling
!      do i = 1, ndim
!         do k = begA(i), begA(i+1)-1
!            call findex(jcoA(k),ix,jy,kz,en)
!            if(landm(ix,jy,kz) == OCEAN) then
!               jj = mod(jcoA(k)-1,nun) + 1
!               coA(k) = cs(jj) * coA(k)
!            else
!               if(jcoA(k) .gt. ndimo) then
!                  jj = mod(jcoA(k)-1,nun) + 1
!                  coA(k) = csa(jj) * coA(k)
!               endif
!            endif
!        enddo
!      enddo

      end subroutine rcscale_a

!****************************************
      subroutine scalesol_a(rl, cs, csa)
      implicit none
      include 'usr.com'
      real RL(ndim)
      real cs(nun), csa(nun)
      integer i,ii,jj,kk,ee,ndimo,ix,iy,iz,en

      do i = 1, ndim
         call findex(i,ix,iy,iz,ii)
         if (landm(ix,iy,iz) == OCEAN) then
            rl(i) = cs(ii) * rl(i)
         elseif (landm(ix,iy,iz) == ATMOS) then ! l < iz < l+la
            rl(i) = csa(ii) * rl(i)
         end if
      end do
!      ndimo = nun*n*m*l
!      do i = 1, ndimo
!         call findex(i,ix,jy,kz,ii)
!         if (landm(ix,jy,kz) == OCEAN) then
!            rl(i) = cs(ii) * rl(i)
!         endif
!      enddo
!      do i = ndimo+1, ndim
!         ee = mod(i-1,nun) + 1
!         rl(i) = csa(ee) * rl(i)
!      enddo
      end subroutine scalesol_a

!****************************************
      subroutine scale_eivec_a(eivec, cs, csa, nf)
      USE m_mat
      implicit none
      include 'par.com'
!      include 'mat.com'
      integer nf
      complex eivec(ndim,nf)
      real    cs(nun), csa(nun)
      integer i,ii,j,ndimo
      ndimo = nun*n*m*l
      do i = 1, ndimo
         ii = mod(i-1,nun) + 1
         do j = 1, nf 
            eivec(i,j) = cs(ii) * eivec(i,j)
         enddo
      enddo
      do i = ndimo+1, ndim
         ii = mod(i-1,nun) + 1
         do j = 1, nf 
            eivec(i,j) = csa(ii) * eivec(i,j)
         enddo
      enddo
      end subroutine scale_eivec_a

!****************************************
      subroutine rcscale(rl, rs, cs)
      USE m_mat
      implicit none
      include 'par.com'
!      include 'mat.com'
      real RL(ndim)
      real rs(nun), cs(nun)
      integer i,j,k,ii,jj,kk,ee
!     row scaling
      do i = 1, ndim
         ee = mod(i-1,nun) + 1
         do k = begA(i), begA(i+1)-1
            coA(k) = rs(ee) * coA(k)
         enddo
         rl(i) = rs(ee) * rl(i)
      enddo
!     column scaling
      do i = 1, ndim
         do k = begA(i), begA(i+1)-1
            jj = mod(jcoA(k)-1,nun) + 1
            coA(k) = cs(jj) * coA(k)
         enddo
      enddo
      end subroutine rcscale

!****************************************
      subroutine scalesol(rl, cs)
      implicit none
      include 'par.com'
      real RL(ndim)
      real cs(nun)
      integer i,ii,jj,kk,ee
      do i = 1, ndim
         ee = mod(i-1,nun) + 1
         rl(i) = cs(ee) * rl(i)
      enddo
      end subroutine scalesol

!****************************************
      subroutine scale_eivec(eivec, cs, nf)
      USE m_mat
      implicit none
      include 'par.com'
!      include 'mat.com'
      integer nf
      complex eivec(ndim,nf)
      real    cs(nun)
      integer i,ii,j
      do i = 1, ndim
         ii = mod(i-1,nun) + 1
         do j = 1, nf 
            eivec(i,j) = cs(ii) * eivec(i,j)
         enddo
      enddo
      end subroutine scale_eivec

!****************************************
      subroutine scaling(ndim, nun, begA, jcoA, coA, rs, cs)
      implicit none
      integer ndim, nun
      integer begA(*), jcoA(*)
      real    coA(*)
      real    db(nun,nun)
      integer i,j,k,ii,jj
      real    rs(nun), cs(nun)
      db = 0.
!     find the average diagonal block*
      do i = 1, ndim
         ii = mod(i-1,nun) + 1
         do j = i-ii+1,i-ii+nun
            jj = mod(j-1,nun) + 1
            do k = begA(i), begA(i+1)-1
               if( jcoA(k) == j ) then
                  db(ii,jj) = db(ii,jj) + coA(k)
                  exit
               endif
            enddo
         enddo
      enddo
      db = nun * db / ndim
      print*, ' Average Block '
      call printmat(db, nun)
      call scal(db,nun,rs,cs)
      write(*,'(a14,6(f10.4),a1)') 'Column Scale [', cs, ']'
      write(*,'(a14,6(f10.4),a1)') 'Row    Scale [', rs, ']'
      end subroutine scaling

!****************************************
      subroutine scaling_land(begA, jcoA, coA, rs, cs)
      implicit none
      include 'usr.com'
      integer begA(*), jcoA(*)
      real    coA(*)
      real    db(nun,nun)
      integer i,j,k,ii,jj,ix,iy,iz,jx,jy,jz,en,nl
      real    rs(nun), cs(nun)
      db = 0.
      nl = 0
!     find the average diagonal block excluding land and atmos.
      do i = 1, ndim
         call findex(i,ix,iy,iz,ii)
         if (landm(ix,iy,iz) == OCEAN) then
            if (ii.EQ.PP) nl = nl+1
            do j = begA(i), begA(i+1)-1
               call findex(jcoA(j),jx,jy,jz,jj)
               if ((ix.eq.jx).and.(iy.eq.jy).and.(iz.eq.jz)) then
                  db(ii,jj) = db(ii,jj) + coA(j)
               end if
            end do
         end if
      end do
            
!         call findex(i+1,ix,jy,kz,en)
!         if(landm(ix,jy,kz) == OCEAN) then
!            nl = nl+1
!            do ii = 1, nun
!               do j = i+1,i+nun
!                  jj = mod(j-1,nun) + 1
!                  do k = begA(i+ii), begA(i+ii+1)-1
!                     if( jcoA(k) == j ) then
!                        db(ii,jj) = db(ii,jj) + coA(k)
!                        exit
!                     endif
!                  enddo
!               enddo
!            enddo
!         endif
!      enddo
      db = db / nl
      print*, ' Average Block (nl = ', nl, ')'
      call printmat(db, nun)
      call scal(db,nun,rs,cs)
      write(*,'(a14,6(f10.4),a1)') 'Column Scale [', cs, ']'
      write(*,'(a14,6(f10.4),a1)') 'Row    Scale [', rs, ']'
      end subroutine scaling_land

!****************************************
      SUBROUTINE scal(mat,size,dr,dc)
!     On output dr and dc contain the row and columnscaling, 
!     respectively
      USE m_dgeco
      USE m_dgedi
      INTEGER          size,psize
      DOUBLE PRECISION mat(size,size),dr(size),dc(size)
      PARAMETER (psize=15)
      INTEGER          iwork(psize),i, pvt(size)
      DOUBLE PRECISION Det(2), rcond, rwork(psize),idc,idr
      DOUBLE PRECISION a,b,c

      a = (mat(5,5)+mat(6,6))/2D0; b = abs(mat(5,6)); c = abs(mat(6,5));

      DO i = 1,size
         pvt(i) = i
      END DO
      IF (size.gt. psize) STOP "increase psize in scale.F"
!     
!     Factor the matrix and estimate condition number:
      CALL dgeco (mat, size, pvt, Rcond)
!
!     Error, if matrix exactly singular or Rcond underflows
      IF (Rcond .LE. 0.0D0) STOP "matrix singular??"
!
!     Compute the inverse of a matrix only "job = 01"
!      CALL dgedi (mat,size, pvt, Det, 01)
      CALL dgedi (mat,size, pvt)
!     mat contains the inverse now!

!     The scaling below is special for oceanography problem
      dr(1)=1.0
      dc(1)=1.0
      idc=sqrt(mat(1,1)/mat(2,2))
      dr(2)=1/idc;
      dc(2)=dr(2);
      idr=sqrt(abs(mat(1,1)/mat(4,4)));
      dr(4)=1/idr;
      dc(4)=dr(4)
!     two possibilities
!     idc(3)*idr(3)*mat(3,3)=mat(1,1)
!     or
!     idr(3)*idc(4)*mat(4,3)=mat(1,1)
      if (abs(mat(4,3)) .GT. abs(mat(3,3))) THEN
!     idr contains currently idc(4) 
         idr=mat(1,1)/(idr*mat(4,3));
      else
         idr=sqrt(abs(mat(1,1)/mat(3,3)))
      endif
      dr(3)=2/idr
      dc(3)=dr(3)
      if (abs(mat(4,5)*mat(5,4)).LT..01*abs(mat(4,4)*mat(5,5))) then
         mat(4,5)=1;mat(5,4)=1;
      endif
      idc=sqrt(abs(mat(1,1)*mat(4,5)/(mat(5,4)*mat(5,5))))
      idr=mat(1,1)/(idc*mat(5,5))
      dr(5)=1/idr
      dc(5)=1/idc
      if (abs(mat(4,6)*mat(6,4)).LT..01*abs(mat(4,4)*mat(6,6))) then
         mat(4,6)=1;mat(6,4)=1;
      endif
      idc=sqrt(abs(mat(1,1)*mat(4,6)/(mat(6,4)*mat(6,6))))
      idr=mat(1,1)/(idc*mat(6,6))
      dr(6)=1/idr
      dc(6)=1/idc
!ACdN overrule scaling of T and S 
!      IF ((a.gt.0).and.(b.gt.0)) THEN
!	WRITE(*,*) "overruling scaling"
!      	dr(5) = 1D0/b
!      	dr(6) = 1D0/a
!      	dc(5) = b/a
!      	dc(6) = 1D0
!      ENDIF      
      END SUBROUTINE scal

!*************************************

      END MODULE m_scaling
