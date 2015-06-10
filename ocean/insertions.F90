#include "fdefs.h"

module m_insertions

  use m_usr

  implicit none
  
contains

  subroutine insert_atmosphere(inserted_atmos)

    use, intrinsic :: iso_c_binding
    use m_par  ! for pi?
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_atmos
    real    :: temfun ! from forcing.F90
    integer :: i,j,pos

    if (ite.eq.2) then
       pos = 1
       do j = 1,m
          do i = 1,n
             tatm(i,j) = inserted_atmos(pos)
             ! tatm(i,j) = temfun(x(i),y(j))
             pos = pos + 1
          end do
       end do
    else
       _INFO2_("Unable to insert atmosphere, ite: ", ite)
    end if
  end subroutine insert_atmosphere

end module m_insertions
