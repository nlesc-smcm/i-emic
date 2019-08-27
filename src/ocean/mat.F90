#include "fdebug.h"

MODULE m_mat

  use, intrinsic :: iso_c_binding

  ! defines the location of the matrix
  ! replaces old common block file "mat.com"

  ! originally in usr.com:
  real,    dimension(:,:,:,:,:,:), ALLOCATABLE :: Al, An
  real,    dimension(:,:,:), ALLOCATABLE :: Alocal

  ! originally in mat.com: now allocated in C++ via the
  ! subroutines get_array_sizes and set_pointers
  real(c_double), dimension(:), POINTER :: coA
  integer(c_int), dimension(:), POINTER :: jcoA
  integer(c_int), dimension(:), POINTER :: begA

  integer :: maxnnz !! allocated memory for jacobian matrix entries

  real(c_double), dimension(:), POINTER :: coB

  ! Used for the stochastic forcing
  real(c_double), dimension(:), POINTER :: coF
  integer(c_int), dimension(:), POINTER :: jcoF
  integer(c_int), dimension(:), POINTER :: begF

contains

  !! allocates the Al and An arrays. The
  !! CRS matrix A and B are allocated by C++ via 'allocate_crs' (below)
  subroutine allocate_mat

    use m_usr
    implicit none

    allocate(Al(np,nun,nun,n,m,l))
    allocate(An(np,nun,nun,n,m,l))
    allocate(Alocal(np,nun,nun))

  end subroutine allocate_mat

  subroutine deallocate_mat

    use m_usr
    implicit none

    deallocate(Al)
    deallocate(An)
    deallocate(Alocal)

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
  subroutine set_pointers(nrows,nnz,begC,jcoC,coC,coBC,begFC,jcoFC,coFC)

    implicit none

    integer(c_int) :: nrows,nnz
    integer(c_int), dimension(nrows+1),target :: begC
    integer(c_int), dimension(nnz),target :: jcoC
    real(c_double), dimension(nnz),target :: coC

    real(c_double), dimension(nrows),target :: coBC

    integer(c_int), dimension(nrows+1),target :: begFC
    integer(c_int), dimension(nnz),target :: jcoFC
    real(c_double), dimension(nnz),target :: coFC

    _DEBUG_("module m_mat: set CRS pointers to external memory...")
    _DEBUG2_("nrows: ",nrows)
    _DEBUG2_("nnz: ",nnz)

    begA=>begC(1:nrows+1)
    jcoA=>jcoC(1:nnz)
    coA=>coC(1:nnz)

    coB=>coBC(1:nrows)

    begF=>begFC(1:nrows+1)
    jcoF=>jcoFC(1:nrows)
    coF=>coFC(1:nrows)

  end subroutine set_pointers


END MODULE m_mat
