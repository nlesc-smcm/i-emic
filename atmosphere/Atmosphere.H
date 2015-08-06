//!------------------------------------------------------------------------
/*
  Atmosphere is a model within the I-EMIC. It is a reimplementation of the
  energy balance model already available in the old THCM code.

  Although this is a 2D model we maintain notation of the 3D ocean model.

  --> THIS IMPLEMENTATION SHOULD MOVE TO 0-BASED! better for everyone
  
  Author: Erik -> t.e.mulder@uu.nl
*/
//!------------------------------------------------------------------------

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H


#include <vector>
#include <memory>
#include <ctime> // Timer

#include "SuperVector.H"
#include "MultiArray.H"
#include "AtmosphereDefinitions.H"

class DependencyGrid;
class Atom;
class Timer;

using std::array;
using std::vector;

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
	int np_;  //! #neighbouring points
	int nun_; //! #unknowns

	std::shared_ptr<DependencyGrid> Al_;
	
	double xmin_, xmax_; //! limits in x-direction 
	double ymin_, ymax_; //! limits in y-direction 

	double dx_, dy_; //! grid increments
	
	std::vector<double> xc_;  //! x-grid ( cell centers)
	std::vector<double> xu_;  //! x-grid ( [u,v]-nodes)
	std::vector<double> yc_;  //! y-grid ( cell centers)
	std::vector<double> yv_;  //! y-grid ( [u,v]-nodes)

	// Ocean surface temperature
	std::vector<double> oceanTemp_;

	// Forcing vector
	std::vector<double> frc_;

	// State of the atmosphere,
	// wrapped in a shared_ptr so we can safely hand out
	// a view of this member.
	std::shared_ptr<std::vector<double> > state_;

	// Right hand side of the discretized system of ODE's
	std::shared_ptr<std::vector<double> > rhs_;

	// Solution vector
	std::shared_ptr<std::vector<double> > sol_;

	// CRS matrix arrays storing the Jacobian
	std::vector<double> ico_;
	std::vector<double> jco_;
	std::vector<double> beg_;

	std::shared_ptr<Timer> timer_;

	// Dense matrix storage
	std::vector<double> denseA_;
	// Pivot array for lapack
	int *ipiv_;

	// Continuation parameter
	double    ampl_; //! amplitude of forcing 
	double amplEnd_; //! max amplitude forcing

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
	
	// Necessary for nondimensionalizations:
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
	~Atmosphere();

	void solve(std::shared_ptr<SuperVector> rhs = nullptr);
	void computeRHS();
	void computeJacobian();

	std::shared_ptr<SuperVector> getSolution(char mode = 'C');
	std::shared_ptr<SuperVector> getState(char mode = 'C');
	std::shared_ptr<SuperVector> getRHS(char mode = 'C');
	
	void setState(std::shared_ptr<SuperVector> state)
		{ *state_ = *(state->getAtmosVector()); }

	void setRHS(std::shared_ptr<SuperVector> rhs)
		{ *rhs_   = *(rhs->getAtmosVector()); }

	void   setPar(double value) { ampl_ = value;   }
	double getPar()             { return ampl_;    }
	double getParDestination()  { return amplEnd_; }

	void test();
	void setOceanTemperature(std::vector<double> &sst,
							 double relaxation);

	void dumpState() { writeAll(); }

	// ! Put default fields in atmosphere and ocean
	void idealizedOcean();
	void idealizedState();
	void zeroOcean();
	void zeroState();

	// ! -------------------------------------------------------
	std::shared_ptr<std::vector<double> > getOceanBlock();
	std::shared_ptr<std::vector<int> >    getAtmosRows();

	
private:
	std::shared_ptr<SuperVector> getVector
	(char mode, std::shared_ptr<std::vector<double> > vec);
	
	// ! Apply local discretization
	void discretize(int type, Atom &atom);

	// ! Apply boundary conditions in dependency grid
	void boundaries();

    // ! Create forcing vector
	void forcing();

	
	// ! Perform the dot product of a single row with the state vector
	double matvec(int row);

    //! Defines location of neighbouring grid points
	//! +----------++-------++----------+
	//! | 12 15 18 || 3 6 9 || 21 24 27 |
	//! | 11 14 17 || 2 5 8 || 20 23 26 |
	//! | 10 13 16 || 1 4 7 || 19 22 25 |
	//! |  below   || center||  above   |
	//! +----------++-------++----------+
	//! shift(i,j,k,i2,j2,k2,loc) sets the neighbour at location loc
	//!  w.r.t. the center of the stencil (5) defined above
	void shift(int i, int j, int k,
			   int &i2, int &j2, int &k2, int loc);

	// ! Find row in vector or matrix corresponding to unknown
	// ! XX at grid point (i,j,k).
	// ! Input and output are expected to be 1-based.
	// ! +-----------------------------------------------------+
	// ! | Example for nun=6:                                  |
	// ! | grid point (i,j,k) | row: (  u,  v,  w,  p,  T,  S) |
	// ! |--------------------+--------------------------------|
	// ! |            (1,1,1) |      (  1,  2,  3,  4,  5,  6) |
	// ! |            (2,1,1) |      (  7,  8,  9, 10, 11, 12) |
	// ! |                    |                                |
	// ! |            (n,1,1) |      ( 6*(n-1)+1,  ...   ,6*n) |
	// ! |                    |                                |
	// ! |            (1,2,1) |      ( 6*n+1, ...      ,6*n+6) |
	// ! |            (2,2,1) |      ( 6*n+7,  ...    ,6*n+12) |
	// ! |                    |                                |
	// ! |            (1,1,2) |      ( 6*n*m+1, ...  ,6*n*m+6) |
	// ! +-----------------------------------------------------+
	int find_row(int i, int j, int k, int XX);

	// ! Assemble the dependency grid into a matrix in crs form
	void assemble();

	void buildDenseA();

	// ! Write everything to output files.
	void writeAll();
};

//==================================================================
//! Dependency grid:
//! A multidimensional array describing the dependencies among the unknowns:
//! grid(i,j,k,21,U,V) = c   <=>   d/dt U|(i,j,k) = ... + c * V|(i-1,j+1,k+1)
//! grid(i,j,k,13,U,V) = d   <=>   d/dt U|(i,j,k) = ... + d * V|(i,j-1,k-1)
//==================================================================
class DependencyGrid
{
    MultiArray<double, 6> grid_;
	int n_, m_, l_, np_, nun_;
	
public:
	DependencyGrid(int n, int m, int l, int np, int nun);
	~DependencyGrid();

	double get(int i, int j, int k, int loc, int A, int B);
	void   set(int i, int j, int k, int loc, int A, int B, double value);
	void   set(int const (&range)[8], int A, int B, Atom &atom);
	void   add(double scalar, Atom &atom);

};

//==================================================================
//! Atom:
//! A multidimensional array describing anonymous dependencies among neighbours:
//! atom(i,j,k,21) = c   <=>   d/dt {}|(i,j,k) = ... + c * {}|(i-1,j+1,k+1)
//! atom(i,j,k,13) = d   <=>   d/dt {}|(i,j,k) = ... + d * {}|(i,j-1,k-1)
//==================================================================
class Atom
{
	MultiArray<double, 4> atom_;
	int n_, m_, l_, np_;
	
public:
	Atom(int n, int m, int l, int np);
	~Atom();

	double get(int i, int j, int k, int loc);
	void   set(int i, int j, int k, int loc, double value);
	void   set(int const (&range)[6], int loc, double value);

	void update(double scalarThis,
				double scalarA, Atom &A,
				double scalarB, Atom &B,
				double scalarC, Atom &C);
};

//-----------------------------------------------------------------------------
// Timer class, compatible with macros in GlobalDefinitions.H
class Timer
{
	std::clock_t startTime_;
public:
	Timer()
		:
		startTime_(std::clock())
		{}

	void ResetStartTime()
		{ startTime_ = std::clock();	}

	double ElapsedTime()
		{ return (double) (std::clock() - startTime_) / CLOCKS_PER_SEC; }
};

#endif
