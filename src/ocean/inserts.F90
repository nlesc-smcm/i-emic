#include "fdefs.h"

module m_inserts

  use m_usr

  implicit none

contains

  !!------------------------------------------------------------------
  subroutine insert_atmosphere_t(inserted_atmos_t)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_atmos_t
    integer :: i,j,pos

    if (coupled_atm.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             tatm(i,j) = inserted_atmos_t(pos)
             pos = pos + 1
          end do
       end do
       !_INFO2_("++ Inserting T: tatm(5,5) = ", tatm(5,5))
    else
       _INFO2_("Not inserting atmosphere T : coupled_atm=", coupled_atm)
    end if
  end subroutine insert_atmosphere_t

  !!------------------------------------------------------------------
  subroutine insert_atmosphere_q(inserted_atmos_q)
    
    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_atmos_q
    integer :: i,j,pos

    if (coupled_atm.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             qatm(i,j) = inserted_atmos_q(pos)
             pos = pos + 1
          end do
       end do
    else
       _INFO2_("Not inserting atmosphere q : coupled_atm =", coupled_atm)
    end if
  end subroutine insert_atmosphere_q

  !!------------------------------------------------------------------
  subroutine insert_atmosphere_p(inserted_atmos_p)
    
    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_atmos_p
    integer :: i,j,pos

    if (coupled_atm.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             pfield(i,j) = inserted_atmos_p(pos)
             pos = pos + 1
          end do
       end do
    else
       _INFO2_("Not inserting atmosphere P : coupled_atm=", coupled_atm)
    end if
  end subroutine insert_atmosphere_p

  !!------------------------------------------------------------------
  subroutine insert_emip(inserted_emip)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_emip
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          emip(i,j) = inserted_emip(pos)
          pos = pos + 1
       end do
    end do

  end subroutine insert_emip

  !!------------------------------------------------------------------
  subroutine insert_tatm(inserted_tatm)
    
    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr
    
    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_tatm
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          tatm(i,j) = inserted_tatm(pos)
          pos = pos + 1
       end do
    end do

  end subroutine insert_tatm


end module m_inserts
