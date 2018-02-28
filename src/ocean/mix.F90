! F90 replacement for 'mix.com'

module m_mix

      real vmix_time
      integer,allocatable, dimension(:) :: vmix_row, vmix_col
      integer :: vmix_dim
      integer, allocatable,dimension(:) ::  vmix_ngrp
      integer :: vmix_mingrp, vmix_maxgrp
      integer,allocatable,dimension(:) :: vmix_ipntr, vmix_jpntr
      integer vmix_flag, vmix_temp, vmix_salt
      integer vmix_fix, vmix_out, vmix_diff
      
      ! global number of grid-cells, required to get 
      ! scaling of conv. adj and mixing right in vmix_fun
      integer :: nmlglob
      
      !! stores a value indicating amount of mixing and conv. adj. in each water column
      !! Trilinos-THCM uses this information for load-balancing.
      real, dimension(:,:,:), allocatable :: vmix_counts


contains

subroutine allocate_mix()

   use m_usr
   implicit none

   allocate(vmix_row(ndim*(nun*np+1))); ! could be taken smaller: ndim*(2*15+1)
   allocate(vmix_col(ndim*(nun*np+1))); ! could be taken smaller: ndim*(2*15+1)
   allocate(vmix_ngrp(ndim));
   allocate(vmix_ipntr(ndim+1), vmix_jpntr(ndim+1));
   allocate(vmix_counts(3,n,m));
   vmix_counts=0.0
end subroutine allocate_mix

subroutine deallocate_mix()

   use m_usr
   implicit none

   deallocate(vmix_row); 
   deallocate(vmix_col); 
   deallocate(vmix_ngrp);
   deallocate(vmix_counts);
   deallocate(vmix_ipntr, vmix_jpntr);

end subroutine deallocate_mix

! set the vmix_fix flag (required for controlling mixing and convective adjustment
! during continuation/time-stepping)
subroutine set_vmix_fix(vm_fix)

implicit none

integer :: vm_fix
vmix_fix = vm_fix

end subroutine set_vmix_fix

end module m_mix
