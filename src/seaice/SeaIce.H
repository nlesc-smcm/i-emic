//!---------------------------------------------------------------
/*
  
  SeaIce is a simple 2D thermodynamic sea ice model relying on the
  ocean and atmosphere states. For testing we allow standalone
  capabilities. Setup is similar to the AtmospherePar class.
  
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
#include "Utils.H"
#include "DependencyGrid.H"

// exp, pow, sin
#include <math.h>

class SeaIce : public Model
{

public:
    using ParameterList = Teuchos::RCP<Teuchos::ParameterList>;
    using VectorPtr     = Teuchos::RCP<Epetra_Vector>;

private:
    //! parameterlist
    ParameterList params_;

    //! communicator object
    Teuchos::RCP<Epetra_Comm> comm_;

    //! global grid size
    int nGlob_,mGlob_, dimGlob_;

    //! local grid size
    int nLoc_, mLoc_, dimLoc_;

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

    //! non-overlapping state
    Teuchos::RCP<Epetra_Vector> state_;

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
    
    double taus_;     //! threshold ice thickness
    double epsilon_;  //! approximation steepness
    
    //! background mean values
    double t0o_; //! background ocean temp 
    double t0a_; //! background atmos temp t0a
    double t0i_; //! background seaice temp t0i
    double tvar_; //! ocean temp variation tvar
    double s0_; //! ocean background salinity s0
    double svar_; //! ocean salinity variation svar
    double q0_; //! atmos background humidity q0
    double qvar_; //! atmos humidity variation qvar
    double H0_; //! seaice background thickness H0
    double M0_; //! seaice background mask M0

    //! ice formation constants
    double ch_; //! empirical constant
    double utau_; //! skin friction velocity, ms^{-1}
    double rhoo_; //! sea water density, kg m^{-3}
    double rhoi_; //! ice density, kg m^{-3}
    double rhoa_; //! atmospheric density, kg m^{-3}
    double cpo_; //! sea water heat capacity, W s kg^{-1} K^{-1]
    double Lf_; //! latent heat of fusion of ice, J kg^{-1}
    double Ls_; //! latent heat of sublimation of ice, J kg^{-1}
    double Ic_; //! constant ice conductivity, W m^{-1} K^{-1}

    //! combined parameter
    double zeta_;
    
    //! background heat flux variation
    double Qvar_;

    //! background heat flux
    double Q0_;

    //! sublimation constants, parameters for saturation humidity over
    //! ice
    double c1_;
    double c2_;
    double c3_;
    
    double ce_; //! Dalton number
    double uw_; //! mean atmospheric surface wind speed, ms^{-1}

    //! typical vertical velocity
    double eta_;

    //! Shortwave radiation constants and functions
    double alpha_; //! albedo
    double sun0_;  //! solar constant
    double c0_;    //! atmospheric absorption coefficient
    double Ch_;    //! Ch
    double cpa_;   //! heat capacity

    //! exchange coefficient 
    double muoa_;

    //! Background sublimation and derivatives
    double E0_;
    double dEdT_;
    double dEdq_;

public:
    //! Constructor based on communicator and parameterlist
    SeaIce(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);

    //! destructor
    ~SeaIce() {};

    //! compute right hand side
    void computeRHS();

    //! compute jacobian
    void computeJacobian();

    void applyMatrix(Epetra_MultiVector const &in,
                     Epetra_MultiVector &out);

    void applyMassMat(Epetra_MultiVector const &in,
                      Epetra_MultiVector &out);

    void applyPrecon(Epetra_MultiVector &in,
                     Epetra_MultiVector &out);

    void initializePrec();
    
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


private:
    //! create x_ and y_
    void createGrid();

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

    //! freezing temperature (dominant term)
    double freezingT(double S)
        { return -0.0575 * (S + s0_); }

    //! ice surface temperature (linearized)
    double iceSurfT(double Q, double H, double S)
        { return freezingT(S) - t0i_ + (Q0_*H0_ + H0_*Qvar_*Q + Q0_*H) / Ic_; }

    //! 0-based local find_row function
    int find_row0(int n, int m, int i, int j, int XX)
        { return dof_ * ( j*n + i) + XX - 1; }
    
    //! 1-based local find_row function
    int find_row1(int n, int m, int i, int j, int XX)
        { return dof_ * ( (j-1)*n + (i-1)) + XX; }
};
    
#endif

