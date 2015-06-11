#include "Atmosphere.H"
#include "math.h"
#define PI_ 3.14159265358979323846

Atmosphere::Atmosphere()
{
	
	// Filling the parameters
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

	// Filling the coefficients
	{
		using Parameters;
		coeff_.muoa =  rhoa * ch * cpa * uw;
		coeff_.amua = (arad + brad * t0) / muoa;
		coeff_.bmua =  brad / muoa;
		coeff_.Ai   =  rhoa * hdima * cpa * udim / (r0dim * muoa);
		coeff_.Ad   =  rhoa * hdima * cpa * d0 / (muoa * r0dim * r0dim);
		coeff_.As   =  sun0 * (1 - c0) / (4 * muoa);		
	}

	// Set problem size and domain limits (temporary!)
	n = 16;
	m = 16;
	l = 1;

	xmin = 286 * PI_ / 180;
	xmax = 350 * PI_ / 180;
	ymin = 10  * PI_ / 180;
	ymax = 74  * PI_ / 180;

	dx = (xmax - xmin) / n;
	dy = (ymax - ymin) / m;
	
	// Fill x
	for (int i = 0; i != n+1; ++i)
	{
		xu.push_back(xmin + i * dx);
		xc.push_back(xmin + (i - 0.5) * dx);
	}
	
	// Fill y and latitude-based arrays
	for (int j = 0; j != m+1; ++j)
	{
		yv.push_back( ymin + j * dy );
		yc.push_back( ymin + (j - 0.5) * dy );
		albe.push_back(0.3);
		datc.push_back(0.9 + 1.5 * exp(-12 * yc(j) * yc(j) / PI_ ));
		datv.push_back(0.9 + 1.5 * exp(-12 * yv(j) * yv(j) / PI_ ));
		suna.push_back(As*(1 - .482 * (3 * pow(sin(yc(j)), 2) - 1.) / 2.) *
					   (1 - albe(j)));
	}
}
