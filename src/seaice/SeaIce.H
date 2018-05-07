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

#include "TRIOS_Domain.H"
#include "GlobalDefinitions.H"

class SeaIce
{

public:
    using ParameterList = Teuchos::RCP<Teuchos::ParameterList>;
    using VectorPtr     = Teuchos::RCP<Epetra_Vector>;

private:
    //! parameterlist
    ParameterList params_;

    //! communicator object
    Teuchos::RCP<Epetra_Comm> comm_;

    //! grid size
    int n_,m_,l_, dim_;

    //! info on periodicity is needed to get the correct domain
    //! decomposition
    bool periodic_;
    
    //! degrees of freedom
    int dof_;

    //! grid dimensions
    double xmin_,xmax_,ymin_,ymax_;

    //! TRIOS domain object for parallelization strategy
    Teuchos::RCP<TRIOS::Domain> domain_;

    //! standard map, without ghost nodes (non-overlapping).
    Teuchos::RCP<Epetra_Map> standardMap_;

    //! assembly map, with ghost nodes (overlapping).
    Teuchos::RCP<Epetra_Map> assemblyMap_;

    //! parallel atmosphere state (non-overlapping)
    Teuchos::RCP<Epetra_Vector> state_;

    //! parallel atmosphere state (overlapping)
    Teuchos::RCP<Epetra_Vector> localState_;

    double taus;     // threshold ice thickness
    double epsilon;  // approximation steepness
    
    // background mean values
    double t0o; // background ocean temp 
    double t0a; // background atmos temp t0a
    double t0i; // background seaice temp t0i
    double tvar; // ocean temp variation tvar
    double s0; // ocean background salinity s0
    double svar; // ocean salinity variation svar
    double q0; // atmos background humidity q0
    double qvar; // atmos humidity variation qvar
    double H0; // seaice background thickness H0
    double M0; // seaice background mask M0

    // ice formation constants
    double ch; // empirical constant
    double utau; // skin friction velocity, ms^{-1}
    double rhoo; // sea water density, kg m^{-3}
    double rhoi; // ice density, kg m^{-3}
    double rhoa; // atmospheric density, kg m^{-3}
    double cpo; // sea water heat capacity, W s kg^{-1} K^{-1]
    double Lf; // latent heat of fusion of ice, J kg^{-1}
    double Ls; // latent heat of sublimation of ice, J kg^{-1}
    double Ic; // constant ice conductivity, W m^{-1} K^{-1}

    // combined parameter
    double zeta;
    
    // sublimation constants, parameters for saturation humidity over
    // ice
    double c1;
    double c2;
    double c3;


public:
    //! Constructor based on communicator and parameterlist
    SeaIce(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);

    //! destructor
    ~SeaIce() {};

    void computeRHS();

private:
    
    void computeLocalRHS();
    
    
};

#endif

