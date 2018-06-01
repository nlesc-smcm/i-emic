! Here we collect a few ice params. These will eventually be overwritten
! during synchronization in usrc::set_seaice_parameters.
 
module m_ice
  
  real :: zeta = 0.0      ! combination of sea ice parameters
  real :: a0   = -0.0575  ! freezing temperature S sensitivity
  real :: Lf   = 3.347e5  ! latent heat of fusion of ice
  real :: Qvar = 0.0      ! latent heat of fusion of ice
  real :: Q0   = 0.0      ! latent heat of fusion of ice

 
contains

end module m_ice
