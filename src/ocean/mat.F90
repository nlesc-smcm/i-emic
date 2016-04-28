#include "fdebug.h"

MODULE m_mat

  use, intrinsic :: iso_c_binding

  ! defines the location of the matrix
  ! replaces old common block file "mat.com" 

  ! originally in usr.com:
  real,    dimension(:,:,:,:,:,:), ALLOCATABLE :: Al

  ! originally in mat.com: now allocated in C++ via the
  ! subroutines get_array_sizes and set_pointers
  real(c_double), dimension(:), POINTER :: coA
  integer(c_int), dimension(:), POINTER :: jcoA 
  integer(c_int), dimension(:), POINTER :: begA

  integer :: maxnnz !! allocated memory for jacobian matrix entries

  real(c_double), dimension(:), POINTER :: coB        

  logical, dimension(:,:,:), ALLOCATABLE :: active
  ! used to keep track on active couplings

contains

  !! allocates the Al array and 'active'. The
  !! CRS matrix A and B are allocated by C++ via 'allocate_crs' (below)
  subroutine allocate_mat

    use m_usr
    implicit none

    allocate(Al(n,m,l+la,np,nun,nun))
    allocate(active(np,nun,nun))

  end subroutine allocate_mat

  subroutine deallocate_mat

    use m_usr
    implicit none

    deallocate(Al)
    deallocate(active)

  end subroutine deallocate_mat

  !! ask for the dimensions of the CSR arrays
  subroutine get_array_sizes(nrows,nnz)

    use m_usr

    implicit none

    integer :: nrows,nnz

    nrows = ndim
    nnz = ndim*(nun*np+1)
    maxnnz = nnz

  end subroutine get_array_sizes


  !! set pointers in fortran to arrays in C++
  !! so the caller will have access to the data
  subroutine set_pointers(nrows,nnz,begC,jcoC,coC,coBC)

    implicit none

    integer(c_int) :: nrows,nnz
    integer(c_int), dimension(nrows+1),target :: begC
    integer(c_int), dimension(nnz),target :: jcoC
    real(c_double), dimension(nnz),target :: coC
    real(c_double), dimension(nrows),target :: coBC

    _DEBUG_("module m_mat: set CRS pointers to external memory...")
    _DEBUG2_("nrows: ",nrows)
    _DEBUG2_("nnz: ",nnz)

    begA=>begC(1:nrows+1)
    jcoA=>jcoC(1:nnz)
    coA=>coC(1:nnz)
    coB=>coBC(1:nrows)

  end subroutine set_pointers


END MODULE m_mat
