//!---------------------------------------------------------------
/*

  SeaIce is a simple 2D thermodynamic sea ice model relying on the
  ocean and atmosphere states. For testing we allow standalone
  capabilities. Setup is similar to the Atmosphere class, framework is
  partly defined by Model.H.

*/
//!---------------------------------------------------------------

#ifndef SEAICE_H
#define SEAICE_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_CrsGraph.h>
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>

#include "Model.H"
#include "TRIOS_Domain.H"
#include "GlobalDefinitions.H"
#include "SeaIceDefinitions.H"
#include "Utils.H"

// exp, pow, sin
#include <math.h>

#include <functional>

// forward decl
class Ocean;
class Atmosphere;
class DependencyGrid;

class SeaIce : public Model
{

public:

    //! Parameters we need to communicate (see datamembers for
    //! explanation). The same struct is defined in the ocean fortran
    //! code (usrc.F90) to accomodate communication.
    struct CommPars
    {
        double zeta;
        double a0;
        double Lf;
        double s0;
        double rhoo;
        double Qvar;
        double Q0;
    };

    enum FluxID
    {
        _QSW = 0, // Shortwave radiative heat flux
        _QSH,     // Sensible heat flux
        _QLH      // Latent heat flux
    };

protected:
    //! parameterlist
    ParameterList params_;

    //! global grid size
    int nGlob_,mGlob_, dimGlob_;

    //! auxiliary unknowns (integral correction gamma)
    int aux_;

    //! local grid size
    int nLoc_, mLoc_, dimLoc_;

    //! Flux containers
    std::vector<double> QSos_;
    std::vector<double> EmiP_;

    //! info on periodicity is needed to get the correct domain
    //! decomposition
    bool periodic_;

    //! degrees of freedom
    int dof_;

    //! global grid dimensions
    double xmin_,xmax_,ymin_,ymax_;

    //! local grid dimensions
    double xminLoc_,xmaxLoc_,yminLoc_,ymaxLoc_;

    //! local grid increments
    double dx_, dy_;

    //! total integration area
    double totalArea_;

    //! Finite difference increment for lazy mask function derivative
    double df_ = 1e-6;

    //! local grid
    std::vector<double> x_, y_;

    //! TRIOS domain object for parallelization strategy
    Teuchos::RCP<TRIOS::Domain> domain_;

    //! standard map, without ghost nodes (non-overlapping).
    Teuchos::RCP<Epetra_Map> standardMap_;

    //! standard, single unknown surface map, without ghost nodes
    //! (non-overlapping) .
    Teuchos::RCP<Epetra_Map> standardSurfaceMap_;

    //! assembly map, with ghost nodes (overlapping).
    Teuchos::RCP<Epetra_Map> assemblyMap_;

    //! assembly, single unknown surface map, with ghost nodes
    //! (overlapping).
    Teuchos::RCP<Epetra_Map> assemblySurfaceMap_;

    //! non-overlapping rhs
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! non-overlapping mass matrix
    Teuchos::RCP<Epetra_Vector> diagB_;

    //! non-overlapping solution vector
    Teuchos::RCP<Epetra_Vector> sol_;

    //! non-overlapping sst
    Teuchos::RCP<Epetra_Vector> sst_;

    //! non-overlapping sss
    Teuchos::RCP<Epetra_Vector> sss_;

    //! non-overlapping tatm
    Teuchos::RCP<Epetra_Vector> tatm_;

    //! non-overlapping qatm
    Teuchos::RCP<Epetra_Vector> qatm_;

    //! non-overlapping patm
    Teuchos::RCP<Epetra_Vector> patm_;

    //! non-overlapping albedo
    Teuchos::RCP<Epetra_Vector> albe_;


    //! overlapping localState
    Teuchos::RCP<Epetra_Vector> localState_;

    //! overlapping localRHS
    Teuchos::RCP<Epetra_Vector> localRHS_;

    //! overlapping localDiagB
    Teuchos::RCP<Epetra_Vector> localDiagB_;

    //! overlapping localSol
    Teuchos::RCP<Epetra_Vector> localSol_;

    //! overlapping localSST
    Teuchos::RCP<Epetra_Vector> localSST_;

    //! overlapping localSSS
    Teuchos::RCP<Epetra_Vector> localSSS_;

    //! overlapping localAtmosT
    Teuchos::RCP<Epetra_Vector> localAtmosT_;

    //! overlapping localAtmosQ
    Teuchos::RCP<Epetra_Vector> localAtmosQ_;

    //! overlapping localAtmosP
    Teuchos::RCP<Epetra_Vector> localAtmosP_;

    //! overlapping localAtmosA
    Teuchos::RCP<Epetra_Vector> localAtmosA_;

    //! non-overlapping integral coefficients
    Teuchos::RCP<Epetra_Vector> intCoeff_;

    //! non-overlapping integral coefficients
    Teuchos::RCP<Epetra_Vector> localIntCoeff_;

    //! gathered global integral coefficients;
    Teuchos::RCP<Epetra_Vector> globalIntCoeff_;

    //! global surface land  mask
    std::shared_ptr<std::vector<int> > surfmask_;

    //! local surface land mask
    Teuchos::RCP<Epetra_IntVector> localSurfmask_;

    //! State component maps
    std::map<int, Teuchos::RCP<Epetra_Map> > Maps_;

    //! State component importers
    std::map<int, Teuchos::RCP<Epetra_Import> > Imps_;

    //! matrix dependency grid
    std::shared_ptr<DependencyGrid> Al_;

    //! matrix graph, necessary to build Jacobian matrix
    Teuchos::RCP<Epetra_CrsGraph> matrixGraph_;

    //! Jacobian matrix
    Teuchos::RCP<Epetra_CrsMatrix> jac_;

    //! ifpack preconditioner object
    Teuchos::RCP<Ifpack_Preconditioner> precPtr_;

    //! CRS matrix arrays storing the Jacobian
    std::vector<double> co_;
    std::vector<int> jco_;
    std::vector<int> beg_;

    //! preconditioning initialization flag
    bool precInitialized_;

    //! preconditioning computation flag
    bool recomputePrec_;

    //! mass matrix computation flag
    bool recompMassMat_;

    double taus_;     //! threshold ice thickness
    double epsilon_;  //! approximation steepness

    //! background mean values
    double t0o_;  //! background ocean temp
    double t0i_;  //! background seaice temp t0i
    double t0a_;  //! background atmos temp t0a
    double s0_;   //! ocean background salinity s0
    double q0_;   //! atmos background humidity q0
    double qdim_; //! atmos humidity scale
    double tdim_; //! temperature scale
    double H0_;   //! seaice background thickness H0
    double M0_;   //! seaice background mask M0

    //! ice formation constants
    double ch_;   //! empirical constant
    double utau_; //! skin friction velocity, ms^{-1}
    double rhoo_; //! sea water density, kg m^{-3}
    double rhoi_; //! ice density, kg m^{-3}
    double rhoa_; //! atmospheric density, kg m^{-3}
    double cpo_;  //! sea water heat capacity, W s kg^{-1} K^{-1]
    double Lf_;   //! latent heat of fusion of ice, J kg^{-1}
    double Ls_;   //! latent heat of sublimation of ice, J kg^{-1}
    double Ic_;   //! constant ice conductivity, W m^{-1} K^{-1}
    double a0_;   //! freezing temperature sensitivity

    //! combined parameter
    double zeta_;

    //! ocean model
    double r0dim_;
    double udim_;

    //! background heat flux variation
    double Qvar_;

    //! background heat flux
    double Q0_;

    //! sublimation constants, parameters for saturation humidity over
    //! ice
    double c1_;
    double c2_;
    double c3_;
    double c4_;
    double c5_;

    //! function to compute background sublimation over ice
    std::function<double(double)> qsi_;

    //! function to compute background sublimation over water
    std::function<double(double)> qso_;

    //! function to compute background sublimation derivatve over ice
    std::function<double(double)> dqsi_;

    //! function to compute background sublimation derivative over water
    std::function<double(double)> dqso_;

    double ce_; //! Dalton number
    double uw_; //! mean atmospheric surface wind speed, ms^{-1}

    //! typical vertical velocity
    double eta_;

    //! Shortwave radiation constants and functions
    double albe0_; //! reference albedo
    double albed_; //! albedo excursion
    double sun0_;  //! solar constant
    double c0_;    //! atmospheric absorption coefficient
    double Ch_;    //! Ch
    double cpa_;   //! heat capacity

    //! exchange coefficient
    double muoa_;

    //! Background sublimation and derivatives
    double E0i_;  // sublimation
    double E0o_;  // evaporation
    double dEdT_; // sublimation derivative to T
    double dEdq_; // sublimation derivative to q

    //! Ocean nondimensionalization prefactor
    double pQSnd_;

    //! Continuation parameters
    std::vector<std::string> allParameters_;

    double      comb_; //! combined forcing
    double      sunp_; //! solar forcing
    double      latf_; //! latent heat forcing
    double     maskf_; //! mask forcing
    double       shf_; //! sensible heat forcing

public:
    //! Constructor based on communicator and parameterlist
    SeaIce(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);

    //! destructor
    ~SeaIce() {};

    std::string name() const { return "seaice"; }
    virtual int modelIdent() const { return 2; }

    //! compute QSos and EmiP fluxes locally
    void computeLocalFluxes(double *state, double *sss, double *sst,
                            double *qatm, double *patm);

    std::vector<Teuchos::RCP<Epetra_Vector> > getFluxes();

    //! compute right hand side
    void computeRHS();

    //! compute jacobian
    void computeJacobian();

    //! compute diagonal mass matrix
    void computeMassMat();

    void applyMatrix(Epetra_MultiVector const &in,
                     Epetra_MultiVector &out);

    void applyMassMat(Epetra_MultiVector const &in,
                      Epetra_MultiVector &out);

    void applyPrecon(Epetra_MultiVector const &in,
                     Epetra_MultiVector &out);

    void initializePrec();

    void initializeState();

    void solve(Teuchos::RCP<Epetra_MultiVector> const &b);

    //! set idealized forcing (idealized external model states)
    void idealizedForcing();

    //! get pointer to jacobian matrix
    Teuchos::RCP<Epetra_CrsMatrix> getJacobian() { return jac_; }

    //! get pointer to CRS struct local Jacobian
    std::shared_ptr<Utils::CRSMat> getLocalJacobian();

    VectorPtr getState(char mode = 'C')
        { return Utils::getVector(mode, state_); }
    VectorPtr getSolution(char mode = 'C')
        { return Utils::getVector(mode, sol_); }
    VectorPtr getRHS(char mode = 'C')
        { return Utils::getVector(mode, rhs_); }
    VectorPtr getMassMat(char mode = 'C')
        { return Utils::getVector(mode, diagB_); }

    void getCommPars(CommPars &parStruct);

    double  getPar(std::string const &parName);
    void    setPar(std::string const &parName, double value);

    //! get total number of continuation parameters
    int npar();

    //! convert index to parameter name
    std::string int2par(int ind) const;

    Utils::MaskStruct getLandMask()
        {
            WARNING("Not implemented", __FILE__, __LINE__);
            return Utils::MaskStruct();
        }

    void setLandMask(Utils::MaskStruct const &mask);

    Teuchos::RCP<TRIOS::Domain> getDomain() { return domain_; }

    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<SeaIce> seaice)
        { return std::make_shared<Utils::CRSMat>(); }

    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Atmosphere> atmos);
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Ocean> ocean);

    void synchronize(std::shared_ptr<SeaIce> seaice) {}
    void synchronize(std::shared_ptr<Atmosphere> atmos);
    void synchronize(std::shared_ptr<Ocean> ocean);

    void pressureProjection(Teuchos::RCP<Epetra_Vector> vec) {}

    Teuchos::RCP<Epetra_Vector> interface(Teuchos::RCP<Epetra_Vector> vec, int XX) const;

    Teuchos::RCP<Epetra_Vector> interfaceH() const;
    Teuchos::RCP<Epetra_Vector> interfaceQ() const;
    Teuchos::RCP<Epetra_Vector> interfaceM() const;
    Teuchos::RCP<Epetra_Vector> interfaceT() const;
    Teuchos::RCP<Epetra_Vector> interfaceG() const;

    // global interface row
    int interface_row(int i, int j, int XX)
        {
            return find_row0(nGlob_, mGlob_, i, j, XX);
        }

    int dof() { return dof_; }

    void buildPreconditioner() {}

    void preProcess();
    void postProcess();

    std::string writeData(bool describe) const;

private:
    //! create x_ and y_
    void createGrid();

    //! create integral coefficients
    void createIntCoeff();

    //! create matrix graph
    void createMatrixGraph();

    //! create local jacobian
    void computeLocalJacobian();

    //! Assemble dependency grid into CRS matrix
    void assemble();

    //! latitudinal dependence shortwave radiation
    //! --> similar to atmos impl: can be factorized
    double shortwaveS(double y)
        { return (1. - .482 * (3. * pow(sin(y), 2) - 1.) / 2.); }

    //! dimensional freezing temperature (dominant term)
    double freezingT(double S)
        { return a0_ * (S + s0_); }

    //! ice surface temperature deviation (linear)
    double iceSurfT(double Q, double H, double S)
        { return freezingT(S) - t0i_ +
                (Q0_*H0_ + H0_*Qvar_*Q + Q0_*H) / Ic_; }

    // //! ice surface temperature deviation (nonlinear)
    // double iceSurfT(double Q, double H, double S)
    //     { return freezingT(S) - t0i_ +
    //             (Q0_*H0_ + H0_*Qvar_*Q + Q0_*H + Q*H) / Ic_; }

    //! sea ice mask as a function thickness anomaly H
    double maskFun(double H)
        { return (1. / 2.) * (1. + tanh( H / epsilon_ )); }

    //! sea ice mask function lazy derivative
    double dMdH(double H)
        { return ( maskFun(H+df_) - maskFun(H) )/df_; }

    //! 0-based local find_row function
    int find_row0(int n, int m, int i, int j, int XX);

    //! 1-based local find_row function
    int find_row1(int n, int m, int i, int j, int XX);

    // HDF5-based save and load functions components other than the
    // state and parameters
    void additionalImports(EpetraExt::HDF5 &HDF5, std::string const &filename){}

    void additionalExports(EpetraExt::HDF5 &HDF5, std::string const &filename);

};

#endif
