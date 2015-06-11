#include "Atmosphere.H"

Atmosphere::Atmosphere()
{
	pars_.rhoa   = 1.25      ; //! atmospheric density 
	pars_.hdima  = 8400.     ; //! atmospheric scale height 
	pars_.cpa    = 1000.     ; //! heat capacity 
	pars_.d0     = 3.1e+06   ; //! constant eddy diffusivity 
	pars_.arad   = 216.0     ; //! radiative flux param A
	pars_.brad   = 1.5       ; //! radiative flux param B
	pars_.sun0   = 1360.     ; //! solar constant 
	pars_.c0     = 0.43      ; //! atmospheric absorption coefficient
	pars_.ce     = 1.3e-03   ; //! exchange coefficient 
	pars_.ch     = 0.94 * ce ; //! exchange coefficient 
	pars_.uw     = 8.5       ; //! mean atmospheric surface wind speed
	pars_.t0     = 15.0      ; //! reference temperature
	pars_.udim   = 0.1e+00   ; //! typical horizontal velocity of the ocean
	pars_.r0dim  = 6.37e+06  ; //! radius of the earth
}
