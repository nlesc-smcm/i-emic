#include "fdefs.h"

module m_probe

  use m_usr

  implicit none

contains

  !!------------------------------------------------------------------
  subroutine get_atmosphere_t(atmos_temp)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: atmos_temp
    integer :: i,j,pos

    if (coupled_atm.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             atmos_temp(pos) = tatm(i,j)
             pos = pos + 1
          end do
       end do
    else
       _INFO_("No coupling, not obtaining any data")
    end if
  end subroutine get_atmosphere_t

  !!------------------------------------------------------------------
  subroutine get_atmosphere_q(atmos_q)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: atmos_q
    integer :: i,j,pos

    if (coupled_atm.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             atmos_q(pos) = qatm(i,j)
             pos = pos + 1
          end do
       end do
    else
       _INFO_("No coupling, not obtaining any data")
    end if
  end subroutine get_atmosphere_q
  
  !!------------------------------------------------------------------
  subroutine get_atmosphere_p(atmos_p)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: atmos_p
    integer :: i,j,pos

    if (coupled_atm.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             atmos_p(pos) = pfield(i,j)
             pos = pos + 1
          end do
       end do
    else
       _INFO_("No coupling, not obtaining any data")
    end if
  end subroutine get_atmosphere_p

end module m_probe
