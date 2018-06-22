//!------------------------------------------------------------------------
/*
  AtmosLocal is the serial component of the atmosphere model in the
  I-EMIC. It is based on the energy balance model in the old THCM
  code. The atmosphere is extended with a vapor anomaly model and
  latent heat fluxes. This implementation is meant to run within a
  parallel domain decomposition.

  Although this is a 2D model we maintain some notation of the 3D
  ocean model.

  Author: Erik -> t.e.mulder@uu.nl
*/
//!------------------------------------------------------------------------

#ifndef ATMOSLOCAL_H
#define ATMOSLOCAL_H

#include <vector>
#include <memory>
#include <array>
#include <map>

#include "Utils.H"
#include "DependencyGrid.H"

#include "AtmosphereDefinitions.H"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

class AtmosLocal
{
public:

    //! Parameters we need to communicate (see datamembers for explanation).
    struct CommPars
    {
        double tdim;
        double qdim;
        double nuq;
        double eta;
        double dqso;
        double dqsi;
        double dqdt;
        double Eo0;
        double t0o;
        double t0i;
        double da;
        double tauf;
    };

private:
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
    int aux_; //! number of auxiliary unknowns    

    std::shared_ptr<DependencyGrid> Al_;

    double xmin_, xmax_;           //! local  limits in x-direction
    double ymin_, ymax_;           //! local  limits in y-direction
    double xmin_glob_, xmax_glob_; //! global limits in x-direction
    double ymin_glob_, ymax_glob_; //! global limits in y-direction

    double dx_, dy_; //! grid increments

    std::vector<double> xc_;  //! x-grid ( cell centers)
    std::vector<double> xu_;  //! x-grid ( [u,v]-nodes)
    std::vector<double> yc_;  //! y-grid ( cell centers)
    std::vector<double> yv_;  //! y-grid ( [u,v]-nodes)

    //! Periodic boundary conditions in x-direction
    bool periodic_;

    //! Find out whether we are part of a domain decomposition
    bool parallel_;

    //! Surface mask
    std::shared_ptr<std::vector<int> > surfmask_;
    
    //! Forcing vector
    std::vector<double> frc_;

    //! State of the atmosphere,
    //! wrapped in a shared_ptr so we can safely hand out
    //! a view of this member.
    std::shared_ptr<std::vector<double> > state_;

    //! Right hand side of the discretized system of ODE's
    std::shared_ptr<std::vector<double> > rhs_;

    //! Solution vector
    std::shared_ptr<std::vector<double> > sol_;

    //! Mass matrix diagonal
    std::shared_ptr<std::vector<double> > diagB_;    

    //! Land surface temperature
    std::shared_ptr<std::vector<double> > lst_;

    //! Ocean surface temperature
    std::shared_ptr<std::vector<double> > sst_;

    //! Sea ice surface temperature
    std::shared_ptr<std::vector<double> > sit_;

    //! Sea ice mask. In the local atmosphere model the sea ice mask
    //! is used to select surface temperature sst or sit. See
    //! forcing() and computeEvaporation().
    std::shared_ptr<std::vector<double> > Msi_;

    //! Integral condition row (only in serial case)
    int rowIntCon_;

    //! Integral condition coefficients (only in serial case)
    std::shared_ptr<std::vector<double> > intcondCoeff_;

    //! Precipitation integral coefficients (only in serial case)
    std::shared_ptr<std::vector<double> > pIntCoeff_;

    //! Total sea surface area
    double totalArea_;

    //! CRS matrix arrays storing the Jacobian
    std::vector<double> co_;
    std::vector<int> jco_;
    std::vector<int> beg_;

    //! Pivot array for lapack
    std::vector<int> ipiv_;

    std::vector<double> bandedA_; //! Banded matrix storage
    int  ksub_;                   //! Number of subdiagonals
    int  ksup_;                   //! Number of superdiagonals
    int  ldimA_;                  //! Leading dimension of banded storage
    bool buildLU_;

    //! Parameters
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
    double   tdim_; //! temperature scale (not used)
    double     q0_; //! reference humidity (kg / kg)
    double   qdim_; //! humidity scale (kg / kg)
    double     lv_; //! latent heat of vaporization (J / kg)

    //! Necessary for nondimensionalizations:
    double   udim_;  //! typical horizontal velocity of the ocean
    double  r0dim_;  //! radius of the earth

    //! Albedo parameters
    double     a0_; //! FIXME todo: reference albedo
    double     da_; //! FIXME todo: albedo excursion
    
    double   tauf_; //! FIXME todo
    double   tauc_; //! FIXME todo

    double     Tm_; //! FIXME todo
    double     Tr_; //! FIXME todo
    double     Pa_; //! FIXME todo
    double    epm_; //! FIXME todo
    double    epr_; //! FIXME todo
    double    epa_; //! FIXME todo

    //! Finite difference increment in albedo function derivatives
    double df_ = 1e-6;

    double Phv_;    //! horizontal Peclet number for vaport transport
    double nuq_;    //! qdim_ * (rhoo_ / rhoa_) * (hdim_ / hdimq_)
    double eta_;    //! (rhoa_ / rhoo_) * ce_ * uw_;

    double qso_;    //! ocean sat. humidity [Bolton, 1980]
    double qsi_;    //! ice   sat. humidity [Bolton, 1980]
    double Eo0_;    //! background evaporation over ocean
    double Ei0_;    //! background evaporation over sea ice
    double Po0_;    //! background precipitation
    double dqso_;   //! derivative sat. humidity w.r.t. T, at t0o (ocean).
    double dqsi_;   //! derivative sat. humidity w.r.t. T, at t0i (ice).

    double lvscale_; //! latent heat due to P coeff

    //! Evaporation field
    std::shared_ptr<std::vector<double> > E_;

    //! Precipitation
    std::shared_ptr<std::vector<double> > P_;

    //! Coefficients
    double muoa_;   //! constant exchange coefficient
    double amua_;   //! (arad+brad*t0)/muoa
    double bmua_;   //!  brad/muoa
    double   Ai_;   //!  rhoa*hdima*cpa*udim/(r0dim*muoa)
    double   Ad_;   //!  rhoa*hdima*cpa*d0/(muoa*r0dim*r0dim)
    double   As_;   //!  sun0*(1 - c0)/(4*muoa)
    double  Ooa_;   //!  ...? some combination of ocean coeffs, see THCM
    double   Os_;   //!  ...? some combination of ocean coeffs, see THCM

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
    double      latf_; //! latent heat forcing

public:
    //! Constructor for use with parallel atmosphere
    AtmosLocal(int n, int m, int l, bool periodic,
               double xmin, double xmax, double ymin, double ymax,
               Teuchos::RCP<Teuchos::ParameterList> params);

    //! Constructor
    AtmosLocal(Teuchos::RCP<Teuchos::ParameterList> params);

    //! Parameter input routine
    void setParameters(Teuchos::RCP<Teuchos::ParameterList> params);
        
    //! General setup things for both constructors
    void setup();

    //! Setup integral coefficients
    void setupIntCoeff();

    //! Destructor
    virtual ~AtmosLocal();

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

    //!  Compute mass matrix
    void computeMassMat();    

    void computeEvaporation();
    void computePrecipitation();

    //! Return the Jacobian matrix in CRS format
    std::shared_ptr<Utils::CRSMat> getJacobian();

    //! Get routines
    std::shared_ptr<std::vector<double> > getSolution    (char mode = 'C');
    std::shared_ptr<std::vector<double> > getState       (char Mode = 'C');
    std::shared_ptr<std::vector<double> > getRHS         (char mode = 'C');
    std::shared_ptr<std::vector<double> > getDiagB       (char mode = 'C');
    std::shared_ptr<std::vector<double> > getSST         (char Mode = 'C');
    std::shared_ptr<std::vector<double> > getE           (char mode = 'C');
    std::shared_ptr<std::vector<double> > getP           (char mode = 'C');
    std::shared_ptr<std::vector<double> > getPIntCoeff   (char mode = 'C');

    void setState(std::shared_ptr<std::vector<double> > in) { state_ = in; }
    void setPrecipitation(std::shared_ptr<std::vector<double> > in) { P_ = in; }

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

    //! Accept sea ice temperature vector
    void setSeaIceTemperature(std::vector<double> const &sit);

    //! Get land temperature set in AtmosLocal::forcing()
    std::shared_ptr<std::vector<double> > getLandTemperature() { return lst_; }

    //! Accept sea ice mask
    void setSeaIceMask(std::vector<double> const &Msi);

    //! Return problem dimension (mainly for testing)
    int  dim() { return dim_; } 

    //! Return number of continuation parameters
    int  npar() { return (int) allParameters_.size(); }

    //! Return parameter name, supply index
    std::string const int2par(int ind);

    //! Return parameter index, supply name
    int const par2int(std::string const &name);

    //! Put default fields in atmosphere and ocean, accepts a value
    //! for the homogeneous precipitation.
    void idealized(double defP = 1.0);

    //! Zero atmosphere state and sst
    void zeroOcean();
    void zeroState();
    
    //! Obtain parameters
    void getCommPars(CommPars &parStruct);
    
    //! Get the dA_{ij} coefficients for surface integrals with q. We
    //! ignore land points by default. Also for the q integral.
    void integralCoeff(std::vector<double> &val,
                       std::vector<double> &ind,
                       int nun = ATMOS_NUN_, bool ignoreLand = true);

    //! --------------------------------------------------------
    //! Set an external 2D surface mask in the atmosphere.
    //! The size of the vector should be n*m.
    void setSurfaceMask(std::shared_ptr<std::vector<int> > surfm);

    //! Compute land temperature based on atmospheric temperature and
    //! shortwave radiative flux. FIXME todo
    // void landTemperature(double Ta, double )

private:

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

    //! For serial use, apply integral condition after assemble
    void intcond();

    //! Write vector to output file
    void write(std::vector<double> &vector, const std::string &filename);

    //! Land temperature is a function of albedo, atmospheric
    //! temperature and latitude grid-box, expecting albedo A, Ta
    //! anomaly and grid point j.
    double Tl(double A, double Ta, int j)
        { return Ta + comb_ * sunp_ * suno_[j] * (1 - A) / Ooa_; }
    
    //! derivative of land temperature w.r.t. albedo
    double dTldA(int j)
        { return -comb_ * sunp_ * suno_[j] / Ooa_; }

    //! derivative of land temperature w.r.t. atmospheric temperature
    double dTldTa()
        { return 1.0; }
    
    //! Snow/ice albedo parametrization expecting albedo A,
    //! atmospheric temperature anomaly Ta, global precipitation
    //! anomaly P and latitude grid-box j
    double aF(double A, double Ta, double P, int j)
        {
            // return H(t0o_ + Tm_ - Tl(A,Ta,j), epm_) *
            //     H(t0o_ + Tr_ - Tl(A,Ta,j), epr_) *
            //     H(Po0_ + eta_ * qdim_ * P - Pa_, epa_);

            // return 5*Ta + 4*A + 3*P;
            
            // return H(Tm_ - Tl(A,Ta,j), epm_) *
            //     H(Tr_ - Tl(A,Ta,j), epr_) *
            //     H(P - Pa_, epa_);

            return H(Tm_ - Tl(A,Ta,j), epm_) *
                H(Tr_ - Tl(A,Ta,j), epr_) * H(P - Pa_, epa_);
        }

    //! aF derivative w.r.t. T, lazy: finite difference
    double daFdA(double A, double Ta, double P, int j)
        { return ( aF(A+df_,Ta,P,j) - aF(A,Ta,P,j) )/df_; }            

    //! aF derivative w.r.t. P, lazy: finite difference
    double daFdP(double A, double Ta, double P, int j)
        { return ( aF(A,Ta,P+df_,j) - aF(A,Ta,P,j) )/df_; }            

    //! aF derivative w.r.t. Ta, lazy: finite difference
    double daFdT(double A, double Ta, double P, int j)
        { return ( aF(A,Ta+df_,P,j) - aF(A,Ta,P,j) )/df_; }            

    //! derivative of atmos temperature equation w.r.t. albedo
    double dTadA(int j)
        { return -comb_ * sunp_ * suna_[j]; }

    //! Heavyside approximation
    double H(double x, double eps)
        { return (1./2.) *  (1.0 + tanh(x / eps)); }

    //! Heavyside approximation derivative
    double dHdx(double x, double eps)
        { return (1./2.) *  (1.0 + tanh(x / eps)); }

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
    

    //! Write everything to output files.
    void writeAll();

private:

};

#endif
