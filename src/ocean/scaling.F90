#include "fdebug.h"

      !! module to compute row and column scaling vectors for the        
      !! Jacobian (to be called from from C++). The procedure is:        
      !! - call 'compute_cell' to obtain a 'sample scaling' for  
      !!   one grid cell                                                 
      !! - broadcast it to everyone so it is the same on all pids        
      !! - call compute_grid to get the actual vectors                
      !!                                                                 
      !! note: the vectors returned are the inverse of those to be       
      !! passed to Trilinos!                                             
      MODULE m_scaling
!

      use m_par, only: nun
      
      private
      public :: average_block, compute
      
!
      CONTAINS
      

!! compute average 6x6 block on the subdomain,           
!! afterwards you should average it over all processes   
!! and call the compute subroutine, passing the globally 
!! averaged block in. The averaging is performed on the  
!! jacobian in m_mat, not the one in Trilinos.           
subroutine average_block(db)
  
use m_usr
use m_mat
implicit none

real, dimension(nun,nun), intent(inout) :: db
integer :: nl, i, ii, j, jj
integer :: ix, iy, iz, jx, jy, jz

      db = 0.
      nl = 0

!     find the average diagonal block excluding land.

! we ignore the top layer of grid cells. If you use TRES=1,SRES=0,
! the top layer leads to different scaling for T and S, which in  
! turn leads to problems when solving ATS.
      do i = 1, ndim !-n*m*nun
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
    if (nl.gt.0) then
      db = db / nl
    end if


end subroutine average_block

!**************************************************************
      
      !! compute scaling, 
      !! construct arrays for full grid
      subroutine compute(db, row_scaling, col_scaling)
      USE m_usr
      implicit none

      real, dimension(nun,nun) :: db
      real, dimension(ndim) :: row_scaling, col_scaling
      
      real, dimension(nun) :: rs,cs,rsa,csa

      integer i,j,k,ii,jj,kk,ee,ndimo,ix,iy,iz

      call scal(db,nun,rs,cs)
                                    
      _DEBUG2_("rs: ",rs)
      _DEBUG2_("cs: ",cs)
      do i = 1, ndim
        call findex(i,ix,iy,iz,ii)

        if (landm(ix,iy,iz) == OCEAN) then
           
          row_scaling(i) = rs(ii)
          col_scaling(i) = cs(ii)

       else ! land cell or something weird
          row_scaling(i) = 1.0
          col_scaling(i) = 1.0
        end if
      end do

      end subroutine compute


!#ifdef HAVE_MRILU
#if 0
!****************************************
! original subroutine scal from scaling.F90 in THCM 7.0
! uses linpack routines contained in MRILU7
      SUBROUTINE scal(mat,size,dr,dc)
!     On output dr and dc contain the row and columnscaling,
!     respectively
      USE m_dgeco
      USE m_dgedi
      INTEGER            :: size
      INTEGER, PARAMETER :: psize=15
      REAL, DIMENSION(size,size) :: mat
      REAL, DIMENSION(size) :: dr,dc
      INTEGER, DIMENSION(psize) :: iwork
      INTEGER :: i
      INTEGER, DIMENSION(size) :: pvt
      REAL, DIMENSION(2) :: Det
      REAL, DIMENSION(psize) :: rwork
      REAL :: rcond, idc,idr
      REAL :: a,b,c

      a = (mat(5,5)+mat(6,6))/2D0; b = abs(mat(5,6)); c = abs(mat(6,5));

      DO i = 1,size
         pvt(i) = i
      END DO
      IF (size.gt. psize) STOP "increase psize in scaling.F90"
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
!       WRITE(*,*) "overruling scaling"
!       dr(5) = 1D0/b
!       dr(6) = 1D0/a
!       dc(5) = b/a
!       dc(6) = 1D0
!      ENDIF
      END SUBROUTINE scal

#else
!****************************************
! modified subroutine scal from scaling.F90 in THCM 6.0
! this routine is a bit messy since it still cotained a lot of linpack
! which is now replaced by lapack
      SUBROUTINE scal(mat,size,dr,dc)
!     On output dr and dc contain the row and columnscaling, 
!     respectively
      implicit none
      
      integer, intent(in) :: size
      real, dimension(size,size), intent(inout) :: mat      
      real, dimension(size), intent(out) :: dr,dc

      integer, parameter :: psize=15
      integer, dimension(psize) :: iwork
      integer, dimension(size) :: pvt
      real, dimension(2) :: Det
      real               :: rcond, Anorm
      real, dimension(psize) :: rwork
      real :: idc,idr
      integer :: lwork !=4*size
      real, dimension(4*size) :: work
      integer :: i
      integer :: lda,info
      integer, dimension(size) :: ipvt

      lwork =4*size
      
! compute infinity-norm of matrix
      Anorm = 0.0
      do i=1,size
        Anorm = max(Anorm,sum(abs(mat(i,:))))
      end do
        

      DO i = 1,size
         pvt(i) = i
      END DO
      IF (size.gt. psize) STOP "increase psize in scale.F"
!     
!     Factor the matrix and estimate condition number:
!      CALL dgeco (mat, size, pvt, Rcond)
       call dgetrf(size,size,mat,size,ipvt,info)
       call dgecon('I', size, mat, size, Anorm, rcond, work, iwork, info )
!
!     Error, if matrix exactly singular or Rcond underflows
      if (1.0+rcond==1.0) then
         write(*,*) ' thcm scaling: diagonal block is singular up to working precision'
         return
      end if
!     Compute the inverse of a matrix only "job = 01"
!      CALL dgedi (mat,size, pvt, Det, 01)
!      call dgedi (mat,size, pvt)
!     mat contains the inverse now!
      call dgetri(size,mat,size,ipvt,work,lwork,info)

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
      END SUBROUTINE scal
#endif

!*************************************

      END MODULE m_scaling
