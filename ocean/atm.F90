
module m_atm


      real, parameter ::  hdima = 8400.
      real, parameter ::   rhoa = 1.25
      real, parameter ::   uatm = 0.0
      real, parameter ::     ce = 1.3e-03
      real, parameter ::     ch = 0.94 * ce  ! exchange coefficient 
      real, parameter ::    cpa = 1000.
      real, parameter ::     uw = 8.5        ! -> mu=2.5
      real, parameter ::     d0 = 3.1e+06    ! constant eddy diffusivity \[D_0\]
      real, parameter ::     c0 = 0.43       ! atmospheric absorption coefficient
      real, parameter :: sigmab = 5.67e-08   ! 
      real, parameter ::   arad = 216.0      ! radiative flux param A
      real, parameter ::   brad = 1.5        ! radiative flux param B
      real, parameter ::   sun0 = 1360.      ! solar constant \[\Sigma_0\]

      real    Ai, Ad, As, Aa, Os, Aoa, Ooa, amua, bmua
      real, allocatable, dimension(:) ::    dat, davt, albe, suna, suno,upa

contains

subroutine allocate_atm(n,m,l)

      allocate(dat(m),davt(0:m),albe(m),suna(m),suno(m),upa(m));

end subroutine allocate_atm

subroutine deallocate_atm()

      deallocate(dat,davt,albe,suna,suno,upa);

end subroutine deallocate_atm

end module m_atm
