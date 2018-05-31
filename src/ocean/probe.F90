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

    ! Coupling with T (sensible heat flux) would need the atmospheric
    ! temperature tatm.
    if (coupled_T.eq.1) then
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

    ! Coupling with T and S would both need atmospheric humidity qatm
    ! (latent heat flux in T eq., E-P in S eq.).
    if ((coupled_T.eq.1).or.(coupled_S.eq.1)) then
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

    ! Coupling with S would need precipitation P (E-P).
    if (coupled_S.eq.1) then
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

  !!------------------------------------------------------------------
  subroutine compute_evap(ocean_evap, un)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr
    use m_atm

    implicit none

    ! import/export
    real(c_double), dimension(m*n)  :: ocean_evap
    real(c_double), dimension(ndim) :: un

    ! local
    integer :: i,j,pos, ctr

    ! external
    integer find_row2

    ! reset vector
    ocean_evap = 0

    ctr = 0;
    if ((coupled_T.eq.1).or.(coupled_S.eq.1)) then
       pos = 1
       do j = 1,m
          do i = 1,n
             if (landm(i,j,l).eq.0) then
                ocean_evap(pos) = ((deltat / qdim) * &
                     dqso * un(find_row2(i,j,l,TT)) - &
                     qatm(i,j))
                ctr = ctr + 1
             endif
             pos = pos + 1
             
          end do
       end do
       
       _INFO2_("  thcm probe:  ctr    = ", ctr)
       _INFO2_("  thcm probe:  eta    = ", eta)
       _INFO2_("  thcm probe:  deltat = ", deltat)
       _INFO2_("  thcm probe:  qdim   = ", qdim)
       _INFO2_("  thcm probe:  dqso   = ", dqso)

    else
       _INFO_("No coupling, not obtaining any data")
    end if

  end subroutine compute_evap

  !!------------------------------------------------------------------
  subroutine get_emip(cemip)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: cemip
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          cemip(pos) = emip(i,j) * (1 - landm(i,j,l))
          pos = pos + 1
       end do
    end do
  end subroutine get_emip

  !!------------------------------------------------------------------
  subroutine get_adapted_emip(cemip)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: cemip
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          cemip(pos) = adapted_emip(i,j)
          pos = pos + 1
       end do
    end do
  end subroutine get_adapted_emip

  !!------------------------------------------------------------------
  subroutine get_emip_pert(cspert)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: cspert
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          cspert(pos) = spert(i,j)
          pos = pos + 1
       end do
    end do
  end subroutine get_emip_pert
  
  !!------------------------------------------------------------------
  subroutine get_salflux(un, salflux)
    use, intrinsic :: iso_c_binding
    use m_par
    use m_usr
    use m_atm
    implicit none
    real(c_double), dimension(m*n) :: salflux
    integer :: i,j,pos
    real gamma, dedt

    !     IMPORT/EXPORT
    real(c_double),dimension(ndim) ::    un
    
    !     LOCAL
    real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
    real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
    real    T(0:n+1,0:m+1,0:l+la+1), S(0:n+1,0:m+1,0:l+la+1)
    
    gamma = par(COMB) * par(SALT) 
    dedt  = (deltat / qdim) * dqso

    call usol(un,u,v,w,p,T,S)
    
    pos = 1
    do j = 1,m
       do i = 1,n
          if (landm(i,j,l).eq.OCEAN) then
             if (coupled_S.eq.1) then
                salflux(pos) = nus * dedt * T(i,j,l) / gamma -  &
                     nus * ( qatm(i,j) + pfield(i,j) ) / gamma
             else
                salflux(pos) = (1 - SRES + SRES*par(BIOT)) * emip(i,j) - &
                     SRES * par(BIOT) * S(i,j,l) / gamma
             endif
          endif
          pos = pos + 1
       end do
    end do

  end subroutine get_salflux

  !!------------------------------------------------------------------
  subroutine get_temflux(un, temflux)
    use, intrinsic :: iso_c_binding
    use m_par
    use m_usr
    implicit none
    real(c_double), dimension(m*n) :: temflux
    integer :: i,j,pos
    real etabi

    !     IMPORT/EXPORT
    real(c_double),dimension(ndim) ::    un

    !     LOCAL
    real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
    real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
    real    T(0:n+1,0:m+1,0:l+la+1), S(0:n+1,0:m+1,0:l+la+1)

    etabi = par(COMB)*par(TEMP)

    call usol(un,u,v,w,p,T,S)

    pos = 1
    do j = 1,m
       do i = 1,n
          if (landm(i,j,l).eq.OCEAN) then
             if (coupled_T.eq.0) then
                temflux(pos) = (1 - TRES + TRES*par(BIOT)) * tatm(i,j) - &
                     TRES * par(BIOT) * T(i,j,l) / etabi
             else
                ! not sure what this would be
             endif
          endif
          pos = pos + 1
       end do
    end do

  end subroutine get_temflux

end module m_probe
