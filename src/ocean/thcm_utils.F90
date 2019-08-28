#if 0
/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#endif
      !! some remains from THCM, i.e. bgs_utils, postp etc that implement tasks
      !! which require, for instance, access to the staggered grid arrays yv etc.
      !! These remains are called from C++ to avoid having to re-implement them
      !! (they are pretty much independent on anything except standard THCM)
      MODULE m_thcm_utils
      
      use m_usr
      implicit none
      
      real, dimension(:,:,:), pointer :: aliasU, aliasV
      real, dimension(:,:,:), pointer :: aliasW
      real, dimension(:,:,:), pointer :: aliasP, aliasT, aliasS
!
      CONTAINS


!**************************************************************

!**************************************************************
!! build the depth-averaging matrix M_{uv} for the velocity
!! in/out: CSR arrays for Mzuv
      SUBROUTINE build_Mzuv(nrows, ncols, nnz, begM,jcoM,coM)
!     constructs the depth-averaging matrix Mzuv
!     such that Mzuv*Mzuv'=Id

      use m_usr
      IMPLICIT none

!     IMPORT/EXPORT
      integer :: nrows, ncols, nnz
      integer, dimension(nrows+1) :: begM
      integer, dimension(nnz) :: jcoM
      real, dimension(nnz) :: coM

!     LOCAL
      INTEGER                        :: nr,nl,ier,i,j,k, dumuv,svs
      INTEGER                        :: rowu, rowv, colu, colv, idx
!     

      nr = ncols   ! real dimension of the problem
      nl = nr/l    ! dimension of one layer     
      begM(nl+1) = nr+1
      DO i = 1,nl
         begM(i) = (i-1)*l + 1
         DO k = 1,(l)
            idx = (i-1)*l+k
            jcoM(idx) = i + (k-1)*nl 
         END DO
      END DO
      
      coM = 1
      ! set coefficients of landpoints to 0
      DO i = 1,n-1
      DO j = 1,m-1
      loop:DO k = 1,l-1
         IF ((landm(i,j,k)==OCEAN).AND.(landm(i+1,j,k)==OCEAN)) THEN
         IF ((landm(i+1,j+1,k)==OCEAN).AND.(landm(i,j+1,k)==OCEAN)) THEN
            CYCLE loop  ! point is ok, surrounded by OCEAN cells, do not change
         END IF
         END IF
         ! if not interupted, then one of the points is not OCEAN, so set co zero
         rowu = ((j-1)*n +(i-1))*2 + 1
         rowv = ((j-1)*n +(i-1))*2 + 2
         colu = rowu*(l) + k
         colv = rowv*(l) + k
         coM(colu) = 0
         coM(colv) = 0
      END DO loop
      END DO
      END DO
      ! make sure Mzuv Mzuv' = I?
      DO i = 1,nl
         idx = begM(i)
         k = SUM(coM(idx:idx+l-1))
         IF (k.EQ.0) WRITE(6,*) "HELP deling door 0", i
         coM(idx:idx+l-1) = coM(idx:idx+l-1)/SQRT(REAL(k))
      ENDDO

      END SUBROUTINE build_Mzuv
       

!**************************************************************


     !! compute meridional streamfunction. vs is v integrated over x-direction,
     !! the arguments are 1D-pointers representing vs(0:m,1:l) and psiM(0:m,0:l)
     subroutine compute_psim(vs,psim)

#define VS_(j,k) vs((m+1)*((k)-1)+(j)+1)
#define PSIM_(j,k) PsiM((m+1)*(k)+(j)+1)

     use m_usr
     implicit none
     
     real, dimension((m+1)*l) :: vs
     real, dimension((m+1)*(l+1)) :: PsiM
     real :: cs
     integer :: j,k
     
     PsiM = 0.0
     do j=1,m
         cs = cos(yv(j))
         do k=1,l
            if ((z(k)*hdim).lt.(-500)) then
               PSIM_(j,k)=-cs*VS_(j,k)*dz*dfzT(k)  + PSIM_(j,k-1)
            endif
         end do
      end do
     
    end subroutine compute_psim

   !! depth-integrate u-velocity (used by OceanGrid::recomputePsiB())
   subroutine depth_int_u(u,us)

#define US_(i,j) us((n+1)*(j)+(i)+1)
#define U_(i,j,k) u( ((k)*(m+1)+(j))*(n+1) + (i)+1 )

   use m_usr
   implicit none
     
   real, dimension((m+1)*(n+1)) :: us
   real, dimension((m+1)*(n+1)*(l+2)) :: u
   integer :: i,j,k
   real :: dum
     
   ! taken from thcm file bstream.f, subroutine bstream
   us(1:(m+1)*(n+1)) = 0.0
   DO i=0,n
      DO j=0,m
         !write(*,*) '============================================'
         dum = 0.0
         DO k=1,l
            !write(*,*) i,j,k
            dum = U_(i,j,k)*dz*dfzT(k)  + dum
         END DO
         !write(*,*) '==================',i,j,'==================='
         US_(i,j) = dum
         !write(*,*) '============================================'
      END DO
   END DO                                                                          
     
   end subroutine depth_int_u

  !! this is a trick to make 1D C-style arrays look 3D-Fortran-like.
  !! this sub is only used internally.
  subroutine  aliasGrid(u,v,w,p,T,S)

  use m_usr
  implicit none

  real, dimension(1:n+1,1:m+1,1:l+2), target :: u,v
  real, dimension(1:n+2,1:m+2,1:l+1), target :: w
  real, dimension(1:n+2,1:m+2,1:l+2), target :: p,T,S

  ! write(*,*) "============"
  ! write(*,*) u(8,8,8)
  ! write(*,*) "============"
  
  aliasU => u(1:n+1,1:m+1,1:l+2)
  aliasV => v(1:n+1,1:m+1,1:l+2)
  aliasW => w(1:n+2,1:m+2,1:l+1)
  aliasP => p(1:n+2,1:m+2,1:l+2)
  aliasT => T(1:n+2,1:m+2,1:l+2)
  aliasS => S(1:n+2,1:m+2,1:l+2)

!  real, dimension(0:n,0:m,0:l+1), target :: u,v
!  real, dimension(0:n+1,0:m+1,0:l), target :: w
!  real, dimension(0:n+1,0:m+1,0:l+1), target :: p,T,S

!  aliasU => u(0:n,0:m,0:l+1)
!  aliasV => v(0:n,0:m,0:l+1)
!  aliasW => w(0:n+1,0:m+1,0:l)
!  aliasP => p(0:n+1,0:m+1,0:l+1)
!  aliasT => T(0:n+1,0:m+1,0:l+1)
!  aliasS => S(0:n+1,0:m+1,0:l+1)

END SUBROUTINE aliasGrid


      !! calls usol for 1D arrays, i.e. the input (uvwpTS)-vector
      !! is reshaped into 1D arrays u,v,w,p,T,S
      subroutine usol1D(un,u,v,w,p,T,S)
      
      use m_usr
      implicit none

                  
      real, dimension(ndim) :: un
      real, dimension((n+1)*(m+1)*(l+2)), target :: u,v
      real, dimension((n+2)*(m+2)*(l+1)), target :: w
      real, dimension((n+2)*(m+2)*(l+2)), target :: p,T,S

!      real, dimension(:) :: un,u,v,w,p,T,S
      

#ifdef DEBUGGING
      integer :: pos,i,j,k 
#endif      

!            write(*,*) 'un: ',un
!            write(*,*) 'u:  ',U
!            write(*,*) 'v:  ',V

      call aliasGrid(u,v,w,p,T,S)
#ifdef DEBUGGING      
      pos = 1
      write(*,*) 'U/V-ARRAY:'
      do k=1,l+2
        do j=1,m+1
          do i=1,n+1
            write(*,*) aliasU(i,j,k),U(pos)
            write(*,*) aliasV(i,j,k),V(pos)
            write(*,*) ''
            pos=pos+1
          end do
        end do
      end do
      write(*,*) 'W-ARRAY:'
      pos=1
      do k=0,l
        do j=0,m+1
          do i=0,n+1
            write(*,*) aliasW(i,j,k),W(pos)
            pos=pos+1
          end do
        end do
      end do
      write(*,*) 'P/T/S-ARRAY:'
      pos=1
      do k=0,l+1
        do j=0,m+1
          do i=0,n+1
            write(*,*) aliasP(i,j,k),P(pos)
            write(*,*) aliasT(i,j,k),T(pos)
            write(*,*) aliasS(i,j,k),S(pos)
            write(*,*) ''
            pos=pos+1
          end do
        end do
      end do
      WRITE(*,*) 'call solu...'
#endif      
      call usol(un,aliasU,aliasV,aliasW,aliasP,aliasT,aliasS)      
      !WRITE(*,*) 'done!'           
      
      end subroutine usol1D

  ! extract a copy of the landm array.
  ! returns landm with full index range,
  ! (0:n+1,0:m+1,0:l+1)
  subroutine get_landm(cland)

  implicit none

  integer, dimension((n+2)*(m+2)*(l+2)) :: cland

  integer :: i,j,k,pos

  pos = 1
  do k=0,l+1
    do j=0,m+1
      do i=0,n+1
        cland(pos) = landm(i,j,k)
        pos=pos+1
      end do
    end do
  end do

  end subroutine get_landm

!! subroutine to get scaling for integral condition

!! If SRES==0 we need to satisfy \int S = 0 which amounts
!! to subtracting a constant from S after solving the ATS
!! system inside the preconditioner. For the integral we 
!! need some scaling factors because of metric terms.    
subroutine intcond_scaling(val,ind,len)

implicit none

  integer :: find_row2

  real, dimension(n*m*l) :: val
  integer, dimension(n*m*l) :: ind
  integer :: len
  
  integer :: i,j,k,v

      ! this code is copied from assemble.f, subroutine intcond
      v = 0
      do k = 1, l
         do j = 1, m
            do i = 1, n
               if ( landm(i,j,k) == OCEAN ) then
                  v = v + 1
                  val(v) = cos(y(j)) * dfzT(k)
                  ind(v) = FIND_ROW2(i,j,k,SS)
               endif
            enddo
         enddo
      enddo
      len=v
      
end subroutine intcond_scaling

!! subroutine that computes weights for load-balancing our application.
!! it does this as follows:                                            
!!   INPUT/OUTPUT: array, double n times m array which will contain    
!! suggested weights for each 'water-column', i.e cells (i,j,1:l).     
!! the weights are computed by assigning the following weights to indi-
!! vidual cells:                                                       
!! 0 - land cells                                                      
!! 1 - standard ocean cells                                            
!! X - extra weights for mixing and convective adjustment, computed by 
!!     mix_imp.f and weighted additionally by the user input to this   
!!     subroutine.                                                     
subroutine loadbal_weights(array,fac_ntrphys,fac_consmix,fac_convadj)

use m_usr
use m_mix

implicit none

real, dimension(n,m) :: array
real :: fac_ntrphys,fac_consmix,fac_convadj
integer :: i,j

! count land as 0, ocean as 1 
do j=1,m
  do i=1,n
    array(i,j) = count(mask=(landm(i,j,:)==OCEAN))
  end do
end do

! currently we use weights of 1 for 
! - neutral physics (index 1)
! - cons. mixing    (index 2)
! - conv. adjustm.  (index 3)
array = array + fac_ntrphys*vmix_counts(1,:,:)
array = array + fac_consmix*vmix_counts(2,:,:)
array = array + fac_convadj*vmix_counts(3,:,:)

array = array/l

end subroutine loadbal_weights

END MODULE m_thcm_utils
