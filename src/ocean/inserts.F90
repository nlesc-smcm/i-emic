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

    if (coupled_T.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             tatm(i,j) = inserted_atmos_t(pos)
             pos = pos + 1
          end do
       end do
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

    if ((coupled_T.eq.1).or.(coupled_S.eq.1)) then
       pos = 1
       do j = 1,m
          do i = 1,n
             qatm(i,j) = inserted_atmos_q(pos)
             pos = pos + 1
          end do
       end do
    end if
  end subroutine insert_atmosphere_q

  !!------------------------------------------------------------------
  subroutine insert_atmosphere_a(inserted_atmos_a)
    
    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_atmos_a
    integer :: i,j,pos

    if (coupled_T.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             albe(i,j) = inserted_atmos_a(pos)
             pos = pos + 1
          end do
       end do
    end if
  end subroutine insert_atmosphere_a
  
  !!------------------------------------------------------------------
  subroutine insert_atmosphere_p(inserted_atmos_p)
    
    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_atmos_p
    integer :: i,j,pos

    if (coupled_S.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             patm(i,j) = inserted_atmos_p(pos)
             pos = pos + 1
          end do
       end do
    end if
  end subroutine insert_atmosphere_p

  !!------------------------------------------------------------------
  subroutine insert_seaice_q(inserted_seaice_q)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_seaice_q
    integer :: i,j,pos
    pos = 1
    do j = 1,m
       do i = 1,n
          qsa(i,j) = inserted_seaice_q(pos)
          pos = pos + 1
       end do
    end do

  end subroutine insert_seaice_q

  !!------------------------------------------------------------------
  subroutine insert_seaice_m(inserted_seaice_m)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_seaice_m
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          msi(i,j) = inserted_seaice_m(pos)
          pos = pos + 1
       end do
    end do

  end subroutine insert_seaice_m

  !!------------------------------------------------------------------
  subroutine insert_seaice_g(inserted_seaice_g)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_seaice_g
    integer :: i,j,pos

    if (coupled_S.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             gsi(i,j) = inserted_seaice_g(pos)
             pos = pos + 1
          end do
       end do
    end if
  end subroutine insert_seaice_g

  !!------------------------------------------------------------------
  !!FIXME superfluous?
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
          emip(i,j) = inserted_emip(pos) * (1 - landm(i,j,l))
          pos = pos + 1
       end do
    end do

  end subroutine insert_emip

  !!------------------------------------------------------------------
  subroutine insert_adapted_emip(inserted_aEmip)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_aEmip
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          adapted_emip(i,j) = inserted_aEmip(pos) * (1 - landm(i,j,l))
          pos = pos + 1
       end do
    end do

  end subroutine insert_adapted_emip

  !!------------------------------------------------------------------
  subroutine insert_emip_pert(inserted_pEmip)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n), intent(in) :: inserted_pEmip
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          spert(i,j) = inserted_pEmip(pos) * (1 - landm(i,j,l))
          pos = pos + 1
       end do
    end do

  end subroutine insert_emip_pert

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
