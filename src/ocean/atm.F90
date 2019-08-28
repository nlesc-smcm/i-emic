! Atmosphere parameters and arrays

module m_atm

      real, parameter ::  hdima = 8400.      ! atmospheric scale height \[H_a\]
      real, parameter ::   rhoa = 1.25       ! atmospheric density \[\rho_a\]
      real, parameter ::   uatm = 0.0        ! advection?? not in the paper
      real, parameter ::     ce = 1.3e-03    ! exchange coefficient (Dalton)
      real, parameter ::     ch = 0.94 * ce  ! exchange coefficient \[C_H\]
      real, parameter ::    cpa = 1000.      ! heat capacity \[C_{pa}\]
      real, parameter ::     uw = 8.5        ! mean atmospheric surface wind speed
                                             ! \[|V_a|\] --> mu \approx 13
      real, parameter ::     d0 = 3.1e+06    ! constant eddy diffusivity \[D_0\]
      real, parameter ::     c0 = 0.43       ! atmospheric absorption coefficient
      real, parameter :: sigmab = 5.67e-08   ! 
      real, parameter ::   arad = 216.0      ! radiative flux param A
      real, parameter ::   brad = 1.5        ! radiative flux param B
      real, parameter ::   sun0 = 1360.      ! solar constant \[\Sigma_0\]
      real, parameter ::     lv = 2.5e+06    ! latent heat of vaporization \[L_v\]
      
      real    qdim, nuq, nus, eta, dqso, eo0, albe0, albed, lvsc
      
      real    Ai, Ad, As, Aa, Aoa, amua, bmua, scorr

      real ::  Ooa = 1.0 ! default value to avoid initialization
                         ! issues, set in atmos_coeff
      real ::   Os = 1.0 ! default value to avoid initialization
                         ! issues, set in atmos_coeff
      real, allocatable, dimension(:) ::  dat, davt, suna, suno,upa

contains

subroutine allocate_atm(m)

      allocate(dat(m),davt(0:m),suna(m),suno(m),upa(m));

end subroutine allocate_atm

subroutine deallocate_atm()

      deallocate(dat,davt,suna,suno,upa);

end subroutine deallocate_atm

end module m_atm
