! Here we collect a few ice params. These will eventually be overwritten
! during synchronization in usrc::set_seaice_parameters.
 
module m_ice
  
  real*8 :: zeta = 0.0        ! combination of sea ice parameters
  real*8 :: a0   = -0.0575    ! freezing temperature S sensitivity
  real*8 :: Lf   = 3.347e+05  ! latent heat of fusion of ice
  real*8 :: Qvar = 0.0        ! 
  real*8 :: Q0   = 0.0        ! 
 
contains
end module m_ice

