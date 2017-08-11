//!------------------------------------------------------------------------
/*
  Atmosphere is a model within the I-EMIC. It is a reimplementation of the
  energy balance model already available in the old THCM code. However,
  this implementation is able to be run within a domain decomposition.

  Although this is a 2D model we maintain notation of the 3D ocean model.

  Author: Erik -> t.e.mulder@uu.nl
*/
//!------------------------------------------------------------------------

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <vector>
#include <memory>
#include <array>
#include <map>

#include "MultiArray.H"
#include "Utils.H"

#include "AtmosphereDefinitions.H"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

class DependencyGrid;
class Atom;

class Atmosphere
{
    //! Parameterlist
    Teuchos::RCP<Teuchos::ParameterList> params_;

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

    // Find out whether we are part of a domain decomposition
    bool parallel_;

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

    // Integral condition coefficients (only in serial case)
    std::shared_ptr<std::vector<double> > intcondCoeff_;

    // Integral condition row (only in serial case)
    int rowIntCon_;

    // CRS matrix arrays storing the Jacobian
    std::vector<double> co_;
    std::vector<int> jco_;
    std::vector<int> beg_;

    // Pivot array for lapack
    std::vector<int> ipiv_;

    std::vector<double> bandedA_; // Banded matrix storage
    int  ksub_;                   // Number of subdiagonals
    int  ksup_;                   // Number of superdiagonals
    int  ldimA_;                  // Leading dimension of banded storage
    bool buildLU_;

    // Preconditioner: 'J': Jacobi
    //                 'D': Full solve using native solving scheme
    //                 'N': None
    char preconditioner_;

    // Parameters
    double   rhoa_; //! atmospheric density \[\rho_a\]
    double   rhoo_; //! oceanic density \[\rho_o\]
    double  hdima_; //! atmospheric scale height \[H_a\]
    double  hdimq_; //! humidity scale height \[H_q\]
    double   hdim_; //! vertical length scale \[H\]
    double    cpa_; //! heat capacity \[C_{pa}\]
    double     d0_; //! constant eddy diffusivity \[D_0\]
    double  kappa_; //! humidity eddy diffusivity \[\kappa\]    
    double   arad_; //! radiative flux param A
    double   brad_; //! radiative flux param B
    double   sun0_; //! solar constant \[\Sigma_0\]
    double     c0_; //! atmospheric absorption coefficient
    double     ce_; //! exchange coefficient (Dalton number)
    double     ch_; //! exchange coefficient \[C_H\] (Stanton number)
    double     uw_; //! mean atmospheric surface wind speed \[|V_a|\]
    double    t0a_; //! reference temperature atmosphere
    double    t0o_; //! reference temperature ocean
    double    t0i_; //! reference temperature ice
    double     q0_; //! reference humidity (kg / kg)
    double   qdim_; //! humidity scale (kg / kg)

    // Necessary for nondimensionalizations:
    double   udim_; //! typical horizontal velocity of the ocean
    double  r0dim_; //! radius of the earth

    double Phv_;     //! horizontal Peclet number for vaport transport
    double nuqeta_;  //! (hdim_ / hdimq_) * ce_ * uw_

    double dqso_;  //! derivative sat. humidity w.r.t. T, at t0o (ocean).
    double dqsi_;  //! derivative sat. humidity w.r.t. T, at t0i (ice).

    double P_;     //! Idealized precipitation

    // Coefficients
    double muoa_;  //! constant exchange coefficient
    double amua_;  //! (arad+brad*t0)/muoa
    double bmua_;  //!  brad/muoa
    double   Ai_;  //!  rhoa*hdima*cpa*udim/(r0dim*muoa)
    double   Ad_;  //!  rhoa*hdima*cpa*d0/(muoa*r0dim*r0dim)
    double   As_;  //!  sun0*(1 - c0)/(4*muoa)
    double  Ooa_;  //!  ...? something in the ocean
    double   Os_;  //!  ...? something in the ocean

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

    //! -------------------------------------------------------
    //! Continuation
    //! -------------------------------------------------------
    std::vector<std::string> allParameters_;
    std::string parName_;

    double      comb_; //! combined forcing
    double      sunp_; //! solar forcing
    double      humf_; //! humidity forcing

public:
    //! Constructor for use with parallel atmosphere
    Atmosphere(int n, int m, int l, bool periodic,
               double xmin, double xmax, double ymin, double ymax,
               Teuchos::RCP<Teuchos::ParameterList> params);

    //! Constructor
    Atmosphere(Teuchos::RCP<Teuchos::ParameterList> params);

    //! General setup things for both constructors
    void setup();

    //! Destructor
    virtual ~Atmosphere();

    //! Solve the matrix
    void solve(std::shared_ptr<std::vector<double> > const &rhs);

    //! Compute and get current residual vector
    std::shared_ptr<std::vector<double> >
    getCurrResVec(std::shared_ptr<std::vector<double> > const &x,
                  std::shared_ptr<std::vector<double> > const &rhs);

    //! Compute the right hand side
    void computeRHS();

    //! Compute the Jacobian matrix
    void computeJacobian();

    //! Return the Jacobian matrix in CRS format
    std::shared_ptr<Utils::CRSMat> getJacobian();

    std::shared_ptr<std::vector<double> > getSolution (char mode = 'C');
    std::shared_ptr<std::vector<double> > getState    (char Mode = 'C');
    std::shared_ptr<std::vector<double> > getRHS      (char mode = 'C');

    void setState(std::shared_ptr<std::vector<double> > in) { state_ = in; }

    //! Apply the Jacobian to a vector and return the result u = J*v
    void applyMatrix(std::vector<double>  const &v,
                     std::vector<double>  &out);

    //! Getters and setter for the continuation parameters
    void    setPar(double value);
    void    setPar(std::string const &parName, double value);
    void    setParName(std::string const &parName) {parName_ = parName; }
    double  getPar();
    double  getPar(std::string const &parName);
    std::string const getParName() { return parName_; }

    //! Accept ocean temperature vector
    void setOceanTemperature(std::vector<double> const &sst);

    //! Return problem dimension
    int  dim() { return dim_; }

    //! Return number of continuation parameters
    int  npar() { return (int) allParameters_.size(); }

    //! Return parameter name, supply index
    std::string const int2par(int ind);

    //! Return parameter index, supply name
    int const par2int(std::string const &name);

    //! Put default fields in atmosphere and ocean
    void idealizedOcean();
    void idealizedState();
    void idealized();
    void zeroOcean();
    void zeroState();

    //! Get the dependency of q on ocean temperature
    double getDqDTo() { return comb_ * humf_ *
            nuqeta_ * dqso_ / qdim_; }

    //! Get the dA_{ij} coefficients for surface integrals with q
    void intcondCoeff(std::vector<double> &val, std::vector<double> &ind);

    //! --------------------------------------------------------
    //! Set an external 2D surface mask in the atmosphere.
    //! The size of the vector should be n*m.
    void setSurfaceMask(std::shared_ptr<std::vector<int> > surfm);

private:

    std::shared_ptr<std::vector<double> > getVector
    (char mode, std::shared_ptr<std::vector<double> > vec);

    //! Apply local discretization
    void discretize(int type, Atom &atom);

    //! Apply boundary conditions in dependency grid
    void boundaries();

    //! Create forcing vector
    void forcing();

    //! Perform the dot product of a single row with the state vector
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

    //! Assemble the dependency grid into a matrix in CRS form
    void assemble();

    //! Write vector to output file
    void write(std::vector<double> &vector, const std::string &filename);

public:
    //! Find row in vector or matrix corresponding to unknown
    //! XX at grid point (i,j,k).
    //! Input and output are expected to be 1-based.
    //! +-----------------------------------------------------+
    //! | Example for nun=6:                                  |
    //! | grid point (i,j,k) | row: (  u,  v,  w,  p,  T,  S) |
    //! |--------------------+--------------------------------|
    //! |            (1,1,1) |      (  1,  2,  3,  4,  5,  6) |
    //! |            (2,1,1) |      (  7,  8,  9, 10, 11, 12) |
    //! |                    |                                |
    //! |            (n,1,1) |      ( 6*(n-1)+1,  ...   ,6*n) |
    //! |                    |                                |
    //! |            (1,2,1) |      ( 6*n+1, ...      ,6*n+6) |
    //! |            (2,2,1) |      ( 6*n+7,  ...    ,6*n+12) |
    //! |                    |                                |
    //! |            (1,1,2) |      ( 6*n*m+1, ...  ,6*n*m+6) |
    //! +-----------------------------------------------------+
    int find_row(int i, int j, int k, int XX);
    int find_surface_row(int i, int j, int k, int XX);

    //! Write everything to output files.
    void writeAll();

private:

};

//==================================================================
//! Dependency grid:
//! A multidimensional array describing the dependencies among the unknowns:
//! grid(i,j,k,21,U,V) = c   <=>   (...) * d/dt U|(i,j,k) = ... + c * V|(i-1,j+1,k+1)
//! grid(i,j,k,13,U,V) = d   <=>   (...) * d/dt U|(i,j,k) = ... + d * V|(i,j-1,k-1)
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
//! atom(i,j,k,21) = c   <=>   (...) * d/dt {}|(i,j,k) = ... + c * {}|(i-1,j+1,k+1)
//! atom(i,j,k,13) = d   <=>   (...) * d/dt {}|(i,j,k) = ... + d * {}|(i,j-1,k-1)
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
                double scalarB, Atom &B);

    void update(double scalarThis,
                double scalarA, Atom &A,
                double scalarB, Atom &B,
                double scalarC, Atom &C);

    void scale(double scalarThis);

};

#endif
