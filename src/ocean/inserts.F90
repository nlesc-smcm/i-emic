#include "fdefs.h"

module m_inserts

  use m_usr

  implicit none

contains

  subroutine insert_atmosphere_t(inserted_atmos_temp)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_atmos_temp
    integer :: i,j,pos

    if (coupled_atm.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             tatm(i,j) = inserted_atmos_temp(pos)
             pos = pos + 1
          end do
       end do
       !_INFO2_("++ Inserting T: tatm(5,5) = ", tatm(5,5))
    else
       _INFO2_("Not inserting atmosphere T : coupled_atm=", coupled_atm)
    end if
  end subroutine insert_atmosphere_t

  subroutine insert_atmosphere_ep(inserted_atmos_ep)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_atmos_ep
    integer :: i,j,pos

    if (coupled_atm.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             epfield(i,j) = inserted_atmos_ep(pos)
             pos = pos + 1
          end do
       end do
       !_INFO2_("++ Inserting EP: epfield(5,5) = ", epfield(5,5))
    else
       _INFO2_("Not inserting atmosphere EP : coupled_atm=", coupled_atm)
    end if
  end subroutine insert_atmosphere_ep

end module m_inserts
