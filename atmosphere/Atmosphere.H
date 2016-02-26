//!------------------------------------------------------------------------
/*
  Atmosphere is a model within the I-EMIC. It is a reimplementation of the
  energy balance model already available in the old THCM code. However,
  this implementation is able to communicate with Jonas Thies' parallel
  Trilinos-THCM code. Communication is done by the CoupledModel class.

  Although this is a 2D model we maintain notation of the 3D ocean model.

  Several options for solving the Jacobian will be implemented, right now
  we have dense 'D' and banded 'B' 

  --> THIS IMPLEMENTATION SHOULD MOVE TO 0-BASED! better for everyone
  
  Author: Erik -> t.e.mulder@uu.nl
*/
//!------------------------------------------------------------------------

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H


#include <vector>
#include <memory>
#include <array>
#include <map>
#include "SuperVector.H"
#include "MultiArray.H"
#include "AtmosphereDefinitions.H"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "GMRESSolver.H"

class DependencyGrid;
class Atom;

using std::array;
using std::vector;

class Atmosphere
{
public:
	using VectorPtr     = std::shared_ptr<SuperVector>;
	using CRSMat        = std::map<std::string, vector<double> > ; 
	using ParameterList = Teuchos::RCP<Teuchos::ParameterList>;

private:
	ParameterList params_;
	
	int n_;   //! #grid-points in east-west (x) direction
	int m_;   //! #grid-points in north-south (y) direction
	int l_;   //! #grid-points in vertical (z) direction
	int dim_; //! full dimension of the problem

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

	// Periodic boundary conditions in x-direction
	bool periodic_;

	// Specify whether we should use a land mask
	bool use_landmask_;
	
	// Surface mask
	std::shared_ptr<std::vector<int> > surfmask_;
	
	// Ocean surface temperature
	std::vector<double> surfaceTemp_;

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

	// Complete jacobian in CRS form
	std::shared_ptr<CRSMat> jac_;

	// Dense matrix storage --> superfluous
	std::vector<double> denseA_;

 	// Pivot array for lapack
	int *ipiv_;

	std::vector<double> bandedA_;      // Banded matrix storage
	int ksub_;                         // Number of subdiagonals
	int ksup_;                         // Number of superdiagonals
	int ldimA_;                        // Leading dimension of banded storage
	bool buildLU_;

	// Solving scheme: 'D': dense
	//                 'B': banded
	//                 'G': GMRES
	char solvingScheme_;

	// Preconditioner: 'J': Jacobi
	//                 'D': Full solve using native solving scheme
	//                 'N': None
	char preconditioner_;
	
	// GMRES solver object
	GMRESSolver<Atmosphere, std::shared_ptr<SuperVector> > gmresSolver_;
	bool gmresInitialized_;


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
	double  Ooa_;  //!  ...? something in the ocean

    //! Latitude-dependent albedo
	std::vector<double> albe_;

	//! Latitudinal dependence of eddy diffusivity:
	//!  0.9 + 1.5 * exp(-12*yc(j)*yc(j)/pi)
	std::vector<double> datc_; //! (defined on cell centers)
	std::vector<double> datv_; //! (defined on [u,v]-nodes)
	
	//! Latitudinal dependence of short-wave radiative flux:
	//!  As*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
	std::vector<double> suna_;

	//!  Os*(1-.482*(3*sin(y(j))**2-1.)/2.)*(1-albe(j))
	std::vector<double> suno_;

	// -------------------------------------------------------
	// Continuation 
	// -------------------------------------------------------
	std::vector<std::string> allParameters_;
	std::string parName_;
	double      comb_; //! combined forcing
	double      sunp_; //! solar forcing    

	// HDF5 input/output files
	std::string const inputFile_;
	std::string const outputFile_;
	
	// Flags for loading/saving states from/to hfd5 file
	bool loadState_;
	bool saveState_;

	// Enable saving every state in multiple formats
	bool storeEverything_;	
	
public:
	// Constructor, specify horizontal grid dimensions
	//  > the number of vertical layers is assumed 1 at this point
	Atmosphere(ParameterList params);

	// Destructor
	virtual ~Atmosphere();

	// Solve the matrix
	// Needs some factorization
	void solve(std::shared_ptr<SuperVector> rhs = std::shared_ptr<SuperVector>());
	void solve(SuperVector const &rhs, SuperVector &out);
		
	// Compute the right hand side
	void computeRHS();

	// Compute the Jacobian matrix
	void computeJacobian();

	// Return the Jacobian matrix in CRS format
	std::shared_ptr<CRSMat> getJacobian();

	std::shared_ptr<SuperVector> getSolution(char mode = 'C');
	std::shared_ptr<SuperVector> getState(char mode = 'C');
	std::shared_ptr<SuperVector> getRHS(char mode = 'C');

	// Apply the Jacobian to a vector and return the result u = J*v
	std::shared_ptr<SuperVector> applyMatrix(SuperVector const &v);
	void applyMatrix(SuperVector const &v, SuperVector &out);

	// Apply the preconditioner inverse to a vector and return
	// the result u = P^{-1}*v
	std::shared_ptr<SuperVector> applyPrecon(SuperVector const &v);
	void applyPrecon(SuperVector const &v, SuperVector &out);
	
	// Getters and setter for the continuation parameters
	void    setPar(double value);
	void    setPar(std::string const &parName, double value);
	void    setParName(std::string const &parName) {parName_ = parName; }
	double  getPar();
	double  getPar(std::string const &parName);

	void test();
	void setOceanTemperature(std::vector<double> const &sst);

	// ! Pre and post-convergence processing
	void preProcess();
	void postProcess();

	// ! Put default fields in atmosphere and ocean
	void idealizedOcean();
	void idealizedState();
	void zeroOcean();
	void zeroState();

	// ! -------------------------------------------------------
	void getOceanBlock(std::vector<double> &values,
					   std::vector<int> &rows);
	
	//! --------------------------------------------------------
	//! Set an external 2D surface mask in the atmosphere.
	//! The size of the vector should be n*m.
	void setSurfaceMask(std::shared_ptr<std::vector<int> > surfm);
	
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

	// GMRES routines
	void initializeGMRES();
	void GMRESSolve(std::shared_ptr<SuperVector> rhs);
	void GMRESSolve(std::shared_ptr<SuperVector> rhs,
					std::shared_ptr<SuperVector> out);
	
	// ! Compute residual
	double computeResidual(std::shared_ptr<SuperVector> rhs);
	double computeResidual(SuperVector const &rhs,
						   SuperVector const &x);

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

	// Create dense matrix storage to use with lapack functions
	void buildDenseA();

	// Write vector to output file
	void write(std::vector<double> &vector, const std::string &filename);

public:
	// ! Write everything to output files.
	void writeAll();
	void writeJacobian(std::string const filename)
		{
			std::stringstream ico_fname, jco_fname, beg_fname;
			ico_fname << filename << ".ico";
			jco_fname << filename << ".jco";
		    beg_fname << filename << ".beg";
			write(ico_, ico_fname.str());
			write(jco_, jco_fname.str());
			write(beg_, beg_fname.str());
		}

	//! public HDF5 load/save routines
	int loadState() { return loadStateFromFile(inputFile_); }	
	int saveState() { return saveStateToFile(outputFile_); }
	
private:
	
	//! HDF5 load/save routines
	int saveStateToFile(std::string const &filename);
	int loadStateFromFile(std::string const &filename);
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

#endif
