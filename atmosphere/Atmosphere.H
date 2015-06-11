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
	int n; //! #grid-points in east-west (x) direction
	int m; //! #grid-points in north-south (y) direction

	double xmin, xmax; //! limits in x-direction 
	double ymin, ymax; //! limits in y-direction 

	double dx, dy; //! grid increments
	
	std::vector<double> xc;  //! x-grid (cell centers)
	std::vector<double> xu;  //! x-grid ([u,v]-nodes)
	std::vector<double> yc;  //! y-grid (cell centers)
	std::vector<double> yv;  //! y-grid ([u,v]-nodes)

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
	//!  0.9 + 1.5 * exp(-12*yc(j)*yc(j)/pi)
	std::vector<double> datc; //! (defined on cell centers)
	std::vector<double> datv; //! (defined on [u,v]-nodes)
	
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
