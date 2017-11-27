//!
//! Ocean is a model within the I-EMIC. It is a wrapper for Jonas Thies'
//!   Trilinos-THCM, using only a fraction of the functionality that is
//!   available in there.
//!
//!   A model is required to have the following member functions:
//!
//!     VectorPtr getState()
//!     VectorPtr getRHS()
//!     VectorPtr getSolution()
//!
//!     double getPar()
//!     void   setPar()
//!
//!     void preProcess()
//!     void postProcess()
//!
//!     void computeRHS()
//!     void computeJacobian()
//!     void solve()
//!
//!   The VectorPtr type that is returned from the get{.} functions will have
//!   to be known by the interface (e.g. Continuation) that calls this class.


#ifndef OCEAN_H
#define OCEAN_H

#include <Teuchos_RCP.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack_Preconditioner.h>

#include "GlobalDefinitions.H"
#include "AtmospherePar.H"
#include "THCMdefs.H"
#include "OceanGrid.H"
#include "Combined_MultiVec.H"
#include "Utils.H"

#include <string>

// forward declarations
class THCM;

namespace TRIOS
{ class Domain; }

class Epetra_Comm;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_MultiVector;
class Epetra_Operator;
class Epetra_CrsMatrix;


class Ocean
{
public:

    using VectorPtr = Teuchos::RCP<Epetra_Vector>;
    using MatrixPtr = Teuchos::RCP<Epetra_CrsMatrix>;
    using PreconPtr = Teuchos::RCP<Ifpack_Preconditioner>;
    using LandMask  = Utils::MaskStruct;

protected:

    //
    Teuchos::RCP<THCM> thcm_;

    VectorPtr state_;
    VectorPtr rhs_;
    VectorPtr sst_;
    VectorPtr massMatrix_;

    MatrixPtr jac_;
    Teuchos::RCP<Epetra_Comm> comm_;

    VectorPtr rowScaling_;
    VectorPtr colScaling_;
    VectorPtr rowScalingRecipr_;
    VectorPtr colScalingRecipr_;

    // 'F' Flexible GMRES (Belos/Trilinos)
    // 'I' IDR
    char solverType_;

    // Tolerance for recomputing preconditioner
    double recompTol_;

    bool   solverInitialized_;
    bool   precInitialized_;
    bool   recompPreconditioner_;

    VectorPtr sol_;

    // grid representation of the state
    Teuchos::RCP<OceanGrid> grid_;

    Teuchos::RCP<Teuchos::ParameterList> solverParams_;

    // Belos flexible GMRES members
    Teuchos::RCP<Teuchos::ParameterList> belosParamList_;
    Teuchos::RCP<Belos::LinearProblem
                 <double, Epetra_MultiVector, Epetra_Operator> > problem_;
    Teuchos::RCP<Belos::BlockGmresSolMgr
                 <double, Epetra_MultiVector, Epetra_Operator> > belosSolver_;

    Teuchos::RCP<Ifpack_Preconditioner> precPtr_;

    // Domain object
    Teuchos::RCP<TRIOS::Domain> domain_;

    // Problem dimensions:
    int N_, M_, L_;

    // Label of current mask
    std::string currentMask_;

    // HDF5 input/output files
    std::string const inputFile_;
    std::string const outputFile_;

    // Flags for loading/saving states from/to hfd5 file
    bool loadState_;
    bool saveState_;
    bool saveMask_;

    //! Enable saving every state in multiple formats
    bool storeEverything_;

    //! Continuation parameter name.
    //! For a list of possible names see THCM::par2int()
    std::string parName_;

    //! Landmask file
    std::string landmaskFile_;

    //! Row map for pressure points P
    Teuchos::RCP<Epetra_Map> mapP_;

    //! Row map for velocity U
    Teuchos::RCP<Epetra_Map> mapU_;

    //! Rows corresponding to singular pressure points
    VectorPtr singRows_;

    //! Index map for surface T values
    Teuchos::RCP<Epetra_BlockMap> tIndexMap_;

    //! Surface temperature importer
    Teuchos::RCP<Epetra_Import> surfaceTimporter_;

    //! Land mask
    Utils::MaskStruct landmask_;

public:

    //! constructor
    Ocean(Teuchos::RCP<Epetra_Comm> Comm,
          Teuchos::RCP<Teuchos::ParameterList> oceanParamList);

    //! destructor
    ~Ocean();

    //! Solve may optionally accept an rhs of VectorPointer type
    void solve(Teuchos::RCP<Epetra_MultiVector> rhs = Teuchos::null);

    //! Calculate explicit residual norm
    double explicitResNorm(VectorPtr rhs);
    void printResidual(VectorPtr rhs);

    void computeRHS();
    void computeJacobian();

    //! All get...() functions that return a VectorPtr accept
    //! a character specifying the access the user has: Copy 'C' or View 'V'
    //! By default this is a copy.
    VectorPtr getSolution(char mode = 'C');
    VectorPtr getState(char mode = 'C');
    VectorPtr getRHS(char mode = 'C');

    //! Binary diagonal for UVTS parts
    Teuchos::RCP<Epetra_Vector> getM(char mode = 'C');

    //! Return pointer to Jacobian
    MatrixPtr getJacobian() {return jac_;}

    //! Return pointer to domain object
    Teuchos::RCP<TRIOS::Domain> getDomain() { return domain_; }

    //! Apply the Jacobian matrix to a vector
    //! out = J*v
    void applyMatrix(Epetra_MultiVector const &v, Epetra_MultiVector &out);

    //! Apply the preconditioner inverse to a vector
    //! out = P^{-1}*v
    void applyPrecon(Epetra_MultiVector const &v, Epetra_MultiVector &out);

    //! Set prec recompute flag
    void recomputePreconditioner() { recompPreconditioner_ = true; }

    //! Build preconditioner
    void buildPreconditioner(bool forceInit = false);

    //! Initialize solver
    void initializeSolver();

    //! Get pointer to preconditioning operator
    PreconPtr getPreconPtr() { return precPtr_; }

    //! The parameter set members wrap the corresponding
    //! Fortran functions. We keep an internal parameter name which
    //! can be overwritten here.
    void setPar(double value);
    void setPar(std::string const &parName, double value);
    void setParName(std::string const &parName) { parName_ = parName; }

    // Return continuation parameter value
    double getPar();
    double getPar(std::string const &parName);
    std::string const getParName() { return parName_; }

    //! Get current landmask
    LandMask getLandMask();

    //! Get global and local (distributed) landmask arrays.
    //! This routine also tests the obtained landmask and it will
    //! therefore also be set in THCM.
    LandMask getLandMask(std::string const & fname);

    //! Check the current mask
    bool isCurrentMask(std::string const &label)
        {
            return (currentMask_ == label);
        }

    //! Set global and local (distributed) landmask arrays
    void setLandMask(LandMask mask, bool global = false);

    //! Apply landmask
    void applyLandMask(LandMask mask, double factor);
    void applyLandMask(VectorPtr x,
                       LandMask mask, double factor);
    void applyLandMask(VectorPtr x,
                       LandMask maskA, LandMask maskB);

    //! Additional monitor for continuation
    bool monitor(){return false;}

    //! Pre and post-convergence processing
    void preProcess();
    void postProcess();

    //! Gather important data to use in continuation summary
    std::string const writeData(bool describe = false);

    //! Get functions for the grid dimensions
    int getNdim() { return N_;}
    int getMdim() { return M_;}
    int getLdim() { return L_;}

    //! Get the degrees of freedom
    int dof() { return _NUN_; }

    //! Get row at the interface, XX is 1-based, i and j are 0-based.
    //! The interface is at surface points.
    int interface_row(int i, int j, int XX)
        { return FIND_ROW2(_NUN_, N_, M_, L_, i, j, L_-1, XX); }

    // ---------------------------------------------------------------
    // The following functions are required for the coupling with other
    // models:

    //! Receive atmosphere data
    //! This routine is called from a CoupledModel object
    void synchronize(std::shared_ptr<AtmospherePar> atmos);

    //! Obtain atmos temp data for debugging
    Teuchos::RCP<Epetra_Vector> getLocalAtmosT();

    //! Obtain atmos humidity for debugging
    Teuchos::RCP<Epetra_Vector> getLocalAtmosQ();

    //! Obtain atmos P data for debugging
    Teuchos::RCP<Epetra_Vector> getLocalAtmosP();

    //! Obtain ocean E data for debugging
    Teuchos::RCP<Epetra_Vector> getLocalOceanE();
    Teuchos::RCP<Epetra_Vector> getE();

    //! Returns the atmos->ocean coupling block in the Jacobian matrix.
    //! I.e. the derivative of our RHS with respect to the atmosphere.
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<AtmospherePar> atmos);

    // Obtain interface (surface) temperature
    Teuchos::RCP<Epetra_Vector> interfaceT();

    // ---------------------------------------------------------------
    // Scale the matrix and rhs of the problem
    void scaleProblem(VectorPtr rhs = Teuchos::null);
    void unscaleProblem(VectorPtr rhs = Teuchos::null);
    VectorPtr getRowScaling();
    VectorPtr getColScaling();

    int loadState() { return loadStateFromFile(inputFile_); }
    int saveState() { return saveStateToFile(outputFile_); }

    void setParameters(Teuchos::RCP<Teuchos::ParameterList> pars);

    void copyFiles();
    void copyFiles(std::string const &fname);
    void copyMask(std::string const &fname);
private:
    // HDF5-based save and load functions for the state and parameters
    int loadStateFromFile(std::string const &name);
    int saveStateToFile(std::string const &name);

    // Write the state of the ocean to traditional fortran out files fort.*
    // Use matlab plot-scripts for visualization
    void printFiles();

    // Initializer members
    void initializeOcean();
    void initializePreconditioner();
    void initializeBelos();

    // This function finds the pressure rows that will give a
    // singular Schur complement.
    // Based on these rows we will try to create a better landmask.
    int analyzeJacobian();

    VectorPtr getVector(char mode, VectorPtr vec);

    void inspectVector(VectorPtr x);
};
#endif
