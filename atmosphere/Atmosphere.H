//!------------------------------------------------------------------------
/*
  Atmosphere is a model within the I-EMIC. It is a reimplementation of the
  energy balance model already available in the old THCM code.

  Author: Erik -> t.e.mulder@uu.nl
*/
//!------------------------------------------------------------------------

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <vector>

class Atmosphere
{
	int n; //! east-west   (x) direction
	int m; //! north-south (y) direction
	int l; //! z direction

	struct Parameters
	{
		double   rhoa; //! atmospheric density \[\rho_a\]
		double  hdima; //! atmospheric scale height \[H_a\]
		double    cpa; //! heat capacity \[C_{pa}\]
		double     d0; //! constant eddy diffusivity \[D_0\]
		double   arad; //! radiative flux param A
		double   brad; //! radiative flux param B
		double   sun0; //! solar constant \[\Sigma_0\]        
		double     c0; //! atmospheric absorption coefficient
		double     ce; //! exchange coefficient 
		double     ch; //! exchange coefficient \[C_H\]
		double     uw; //! mean atmospheric surface wind speed \[|V_a|\] 
		double     t0; //! reference temperature

		// necessary for the nondimensionalization of time:
		double   udim; //! typical horizontal velocity of the ocean
		double  r0dim; //! radius of the earth
	};
	
	Parameters pars_;

	struct Coefficients
	{
		double muoa;  //! constant exchange coefficient
		double amua;  //! (arad+brad*t0)/muoa
		double bmua;  //!  brad/muoa
		double   Ai;  //!  rhoa*hdima*cpa*udim/(r0dim*muoa)
		double   Ad;  //!  rhoa*hdima*cpa*d0/(muoa*r0dim*r0dim)
		double   As;  //!  sun0*(1 - c0)/(4*muoa)
	};

	Coefficients coeff_;

    //! latitude-dependent albedo
	std::vector<double> albe;

	//! latitudinal dependence of eddy diffusivity:
	//!  0.9 + 1.5 * exp(-12*y(j)*y(j)/pi)
	std::vector<double> dat;

	//! latitudinal dependence of short-wave radiative flux:
	//!  As*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
	std::vector<double> suna; 
	
public:
	// constructor
	Atmosphere();
	// destructor
	~Atmosphere();
	
};



#endif
