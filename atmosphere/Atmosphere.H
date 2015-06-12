//!------------------------------------------------------------------------
/*
  Atmosphere is a model within the I-EMIC. It is a reimplementation of the
  energy balance model already available in the old THCM code.

  Although this is a 2D model we keep some notation of the 3D ocean model.
  
  Author: Erik -> t.e.mulder@uu.nl
*/
//!------------------------------------------------------------------------

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <vector>

class Atmosphere
{
	int n_; //! #grid-points in east-west (x) direction
	int m_; //! #grid-points in north-south (y) direction
	int l_; //! #grid-points in vertical (z) direction

	//! #neighbours, for example in a 3D grid:
	//!   stencil np = 27:
	//!   +----------++-------++----------+
    //!   | 12 15 18 || 3 6 9 || 21 24 27 |
    //!   | 11 14 17 || 2 5 8 || 20 23 26 |
    //!   | 10 13 16 || 1 4 7 || 19 22 25 |
    //!   |  below   || center||  above   |
    //!   +----------++-------++----------+
	int np_;  
	int nun_; //! #unknowns

	//! A multidimensional array describing the dependencies among the unknowns:
	//! Al(i,j,k,21,U,V) = c   <=>   d/dt U|(i,j,k) = ... + c * V|(i-1,j+1,k+1)
	//! Al(i,j,k,13,U,V) = d   <=>   d/dt U|(i,j,k) = ... + d * V|(i,j-1,k-1)
	double ******Al_;
	
	double xmin_, xmax_; //! limits in x-direction 
	double ymin_, ymax_; //! limits in y-direction 

	double dx_, dy_; //! grid increments
	
	std::vector<double> xc_;  //! x-grid ( cell centers)
	std::vector<double> xu_;  //! x-grid ( [u,v]-nodes)
	std::vector<double> yc_;  //! y-grid ( cell centers)
	std::vector<double> yv_;  //! y-grid ( [u,v]-nodes)

	// Parameters	
	double   rhoa_; //! atmospheric density \[\rho_a\]
	double  hdima_; //! atmospheric scale height \[H_a\]
	double    cpa_; //! heat capacity \[C_{pa}\]
	double     d0_; //! constant eddy diffusivity \[D_0\]
	double   arad_; //! radiative flux param A
	double   brad_; //! radiative flux param B
	double   sun0_; //! solar constant \[\Sigma_0\]        
	double     c0_; //! atmospheric absorption coefficient
	double     ce_; //! exchange coefficient 
	double     ch_; //! exchange coefficient \[C_H\]
	double     uw_; //! mean atmospheric surface wind speed \[|V_a|\] 
	double     t0_; //! reference temperature
	
	// Necessary for the nondimensionalization of time:
	double   udim_; //! typical horizontal velocity of the ocean
	double  r0dim_; //! radius of the earth

	// Coefficients
	double muoa_;  //! constant exchange coefficient
	double amua_;  //! (arad+brad*t0)/muoa
	double bmua_;  //!  brad/muoa
	double   Ai_;  //!  rhoa*hdima*cpa*udim/(r0dim*muoa)
	double   Ad_;  //!  rhoa*hdima*cpa*d0/(muoa*r0dim*r0dim)
	double   As_;  //!  sun0*(1 - c0)/(4*muoa)

    //! Latitude-dependent albedo
	std::vector<double> albe_;

	//! Latitudinal dependence of eddy diffusivity:
	//!  0.9 + 1.5 * exp(-12*yc(j)*yc(j)/pi)
	std::vector<double> datc_; //! (defined on cell centers)
	std::vector<double> datv_; //! (defined on [u,v]-nodes)
	
	//! Latitudinal dependence of short-wave radiative flux:
	//!  As*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
	std::vector<double> suna_;
	
public:
	// constructor
	Atmosphere();
	// destructor
	~Atmosphere(){}

	void test();
};



#endif
