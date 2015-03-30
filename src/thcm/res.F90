!! former common block file 'res.com'
!! This seems to be for something called
!! 'residue continuation', I don't know
!! if we need this in Trilinos-THCM, currently
!! it's not used.
module m_res

      integer ires
      real    p0
      real,allocatable, dimension(:) :: ures

    contains      
      
      subroutine allocate_res(ndim)
      
      implicit none
      integer :: ndim
      
      ires=0
      p0=0.0
      allocate(ures(ndim));
      ures=0.d0
      
      end subroutine allocate_res


      subroutine deallocate_res()
      
      implicit none
      
      deallocate(ures);
      
      end subroutine deallocate_res
      
end module m_res
