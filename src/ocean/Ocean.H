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
//!     ... some more (see Continuation)
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

#include "Model.H"

#include <string>

// forward declarations
class Atmosphere;
class OceanGrid;
class SeaIce;
class THCM;

namespace TRIOS
{ class Domain; }

class Epetra_Comm;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_MultiVector;
class Epetra_Operator;
class Epetra_CrsMatrix;


class Ocean : public Model
{
public:

    using PreconPtr = Teuchos::RCP<Ifpack_Preconditioner>;
    using LandMask  = Utils::MaskStruct;

protected:
    Teuchos::ParameterList params_;

    Teuchos::RCP<THCM> thcm_;

    VectorPtr rhs_;
    VectorPtr sst_;
    VectorPtr sss_;
    VectorPtr Qsi_;
    VectorPtr Msi_;
    VectorPtr Gsi_;
    VectorPtr massMatrix_;

    MatrixPtr jac_;
    MatrixPtr frc_;
    VectorPtr diagB_;

    VectorPtr rowScaling_;
    VectorPtr colScaling_;
    VectorPtr rowScalingRecipr_;
    VectorPtr colScalingRecipr_;

    bool   solverInitialized_;
    bool   precInitialized_;
    bool   recompPreconditioner_;
    bool   recompMassMat_;

    VectorPtr sol_;

    // grid representation of the state
    mutable Teuchos::RCP<OceanGrid> grid_;

    // Belos flexible GMRES members
    Teuchos::RCP<Belos::LinearProblem
                 <double, Epetra_MultiVector, Epetra_Operator> > problem_;
    Teuchos::RCP<Belos::BlockGmresSolMgr
                 <double, Epetra_MultiVector, Epetra_Operator> > belosSolver_;

    double effort_;
    mutable int effortCtr_;

    Teuchos::RCP<Ifpack_Preconditioner> precPtr_;

    // Domain object
    Teuchos::RCP<TRIOS::Domain> domain_;

    // Problem dimensions:
    int N_, M_, L_;

    // Label of current mask
    std::string currentMask_;

    // Flags for loading/saving states from/to hfd5 file
    bool loadSalinityFlux_;
    bool saveSalinityFlux_;
    bool loadTemperatureFlux_;
    bool saveTemperatureFlux_;

    //! Use legacy fortran output fort.3 and fort.44
    bool useFort3_, useFort44_;

    //! Integrals of the columns of the Jacobian may indicate errors
    //! in the discretization.
    bool saveColumnIntegral_;

    //! Maximum number of mask fixing cycles
    int maxMaskFixes_;

    //! Select Jacobian analysis
    bool analyzeJacobian_;

    //! Row map for pressure points P
    Teuchos::RCP<Epetra_Map> mapP_;

    //! Row map for velocity U
    Teuchos::RCP<Epetra_Map> mapU_;

    //! Rows corresponding to singular pressure points
    VectorPtr singRows_;

    //! Index maps for surface T,S values
    Teuchos::RCP<Epetra_BlockMap> tIndexMap_, sIndexMap_;

    //! Surface temperature and salinity importers
    Teuchos::RCP<Epetra_Import> surfaceTimporter_, surfaceSimporter_;

    //! Land mask
    Utils::MaskStruct landmask_;

public:
    //! constructor
    Ocean(Teuchos::RCP<Epetra_Comm> Comm);
    Ocean(Teuchos::RCP<Epetra_Comm> Comm, Model::ParameterList oceanParamList);
    Ocean(Teuchos::RCP<Epetra_Comm> Comm, Teuchos::ParameterList& oceanParamList);

    //! destructor
    ~Ocean();

    std::string name() const { return "ocean"; }
    int modelIdent() const { return 0; }
    Teuchos::RCP<Epetra_Comm> Comm() const { return comm_; }

    const Teuchos::ParameterList& getParameters();
    void setParameters(Teuchos::ParameterList& oceanParamList);

    static Teuchos::ParameterList getDefaultInitParameters();
    static Teuchos::ParameterList getDefaultParameters();

    //! Solve may optionally accept an rhs of VectorPointer type
    void solve(Teuchos::RCP<const Epetra_MultiVector> rhs = Teuchos::null);

    //! Calculate explicit residual norm
    double explicitResNorm(VectorPtr rhs);
    void printResidual(VectorPtr rhs);

    //! compute rhs (spatial discretization)
    void computeRHS();

    //! compute derivative of rhs
    void computeJacobian();
    void computeForcing();

    //! compute mass matrix
    void computeMassMat();

    //! All get...() functions that return a VectorPtr accept
    //! a character specifying the access the user has: Copy 'C' or View 'V'
    //! By default this is a copy.
    VectorPtr getSolution(char mode = 'C');
    VectorPtr getState(char mode = 'C');
    VectorPtr getRHS(char mode = 'C');
    VectorPtr getMassMat(char mode = 'C');

    //! Binary diagonal for UVTS parts
    Teuchos::RCP<Epetra_Vector> getM(char mode = 'C');

    //! Return pointer to Jacobian
    MatrixPtr getJacobian() {return jac_;}
    MatrixPtr getForcing() {return frc_;}

    //! Return pointer to domain object
    Teuchos::RCP<TRIOS::Domain> getDomain() { return domain_; }

    //! Apply the Jacobian matrix to a vector
    //! out = J*v
    void applyMatrix(Epetra_MultiVector const &v, Epetra_MultiVector &out);

    //! Apply the preconditioner inverse to a vector
    //! out = P^{-1}*v
    void applyPrecon(Epetra_MultiVector const &v, Epetra_MultiVector &out);

    //! Apply the Jacobian's diagonal inverse to a vector
    //! out = D^{-1}*v
    void applyDiagInv(Epetra_MultiVector const &v, Epetra_MultiVector &out);

    //! Apply the mass matrix
    //! out = B*v
    void applyMassMat(Epetra_MultiVector const &v, Epetra_MultiVector &out);

    //! Set prec recompute flag
    void recomputePreconditioner() { recompPreconditioner_ = true; }

    //! Build preconditioner
    void buildPreconditioner(bool forceInit);
    void buildPreconditioner() { buildPreconditioner(false); }

    //! Initialize solver
    void initializeSolver();

    //! Get pointer to preconditioning operator
    PreconPtr getPreconPtr() { return precPtr_; }

    //! The parameter set members wrap the corresponding
    //! Fortran functions.
    void setPar(std::string const &parName, double value);

    //! Return parameter value
    double getPar(std::string const &parName);

    //! Return total number of continuation parameters
    int npar();

    //! Convert integer parameter index to parameter name
    std::string int2par(int ind) const;

    //! Get current landmask
    Utils::MaskStruct getLandMask()
        {
            return getLandMask(false);
        }

    Utils::MaskStruct getLandMask(bool adjustMask);

    //! Get global and local (distributed) landmask arrays.
    //! This routine also tests the obtained landmask and it will
    //! therefore also be set in THCM.
    Utils::MaskStruct getLandMask(std::string const & fname, bool adjustMask = false);

    //! Check the current mask
    bool isCurrentMask(std::string const &label)
        {
            return (currentMask_ == label);
        }

    //! Set global and local (distributed) landmask arrays
    void setLandMask(Utils::MaskStruct const &mask)
        {
            setLandMask(mask, false);
        }

    void setLandMask(Utils::MaskStruct const &mask, bool global);

    //! Apply landmask
    void applyLandMask(Utils::MaskStruct mask, double factor);
    void applyLandMask(VectorPtr x,
                       Utils::MaskStruct mask, double factor);
    void applyLandMask(VectorPtr x,
                       Utils::MaskStruct maskA, Utils::MaskStruct maskB);

    //! Additional monitor for continuation
    bool monitor(){ return false; }

    //! Pre and post-convergence processing
    void preProcess();
    void postProcess();

    //! Gather important data to use in continuation summary
    std::string writeData(bool describe = false) const;

    //! Get functions for the grid dimensions
    int getNdim() { return N_;}
    int getMdim() { return M_;}
    int getLdim() { return L_;}

    //! Get the minimum and maximum of the Meridional streamfunction
    int getPsiM(double &psiMin, double &psiMax) const;

    //! Get coupling information
    int getCoupledT();
    int getCoupledS();

    //! Get salinity flux correction
    double getSCorr();

    //! Get the degrees of freedom
    int dof();

    //! Get row at the interface, XX is 1-based, i and j are 0-based.
    //! The interface is at surface points.
    int interface_row(int i, int j, int XX);

    //! Get salinity integral condition row
    int getRowIntCon();

    //! The following functions are required for the coupling with
    //! other models. We use an abstract base class Model, from which
    //! all other models are derived. CoupledModel knows only the base
    //! class, hence in Model we downcast to specific derived models.

    //! |      ---Model---        |
    //! |     /     |     \       |
    //! | Ocean Atmosphere SeaIce |

    //! Receive atmosphere data
    void synchronize(std::shared_ptr<Atmosphere> atmos);

    //! Receive sea ice data
    void synchronize(std::shared_ptr<SeaIce> seaice);

    //! Meaningless: dummy implementation
    void synchronize(std::shared_ptr<Ocean> ocean) {}

    //! Obtain atmos temp data for debugging
    Teuchos::RCP<Epetra_Vector> getLocalAtmosT();

    //! Obtain atmos humidity for debugging
    Teuchos::RCP<Epetra_Vector> getLocalAtmosQ();

    //! Obtain atmos P data for debugging
    Teuchos::RCP<Epetra_Vector> getLocalAtmosP();

    //! Obtain ocean E data for debugging
    Teuchos::RCP<Epetra_Vector> getLocalOceanE();

    Teuchos::RCP<Epetra_Vector> getSunO();

    //! Returns the ocean-atmos coupling block in the Jacobian
    //! matrix. I.e. the derivative of our RHS with respect to the
    //! atmosphere. The CouplingBlock class builds a parallel coupling
    //! block from the CRS struct. As the Oceans Jacobian is taken
    //! negative in THCM, this couplingblock should be negated
    //! correspondingly. This means that we add values, i.e., -Ooa,
    //! gamma*eta and gamma, where one would expect opposite signs
    //! based on forcing.F90.
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Atmosphere> atmos);

    //! Return the ocean-seaice coupling block (see previous comment).
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<SeaIce> seaice);

    //! Dummy implementation but required by the base class.
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Ocean> ocean)
        { return std::shared_ptr<Utils::CRSMat>(); }

    // Obtain interface surface temperature
    Teuchos::RCP<Epetra_Vector> interfaceT();

    // Obtain interface surface salinity
    Teuchos::RCP<Epetra_Vector> interfaceS();

    // Obtain evaporation as computed by ocean
    Teuchos::RCP<Epetra_Vector> interfaceE();


    // ---------------------------------------------------------------
    // Scale the matrix and rhs of the problem. These things are
    // confusing and not used.
    void scaleProblem(VectorPtr rhs = Teuchos::null);
    void unscaleProblem(VectorPtr rhs = Teuchos::null);
    VectorPtr getRowScaling();
    VectorPtr getColScaling();

    void copyMask(std::string const &fname);

    void dumpBlocks();

    void integralChecks(Teuchos::RCP<Epetra_Vector> state,
                        double &salt_advection,
                        double &salt_diffusion);

    Teuchos::RCP<Epetra_Vector> getColumnIntegral(Teuchos::RCP<Epetra_CrsMatrix> mat,
                                                  bool useSRES = true);

    Teuchos::RCP<Epetra_Vector> getIntCondCoeff();

    // This function finds the pressure rows that will give a singular
    // Schur complement. Based on these rows we will try to create a
    // better landmask.
    int analyzeJacobian1();

    // Here we integrate columns of the salinity Jacobian over the ocean
    // volume. Nonzero integrals correspond to points in the mask we
    // do not want.
    int analyzeJacobian2();

    // Project pressures modes from a state vector
    void pressureProjection(Teuchos::RCP<Epetra_Vector> vec);

private:
    // HDF5-based save and load functions to load and save components
    // other than the state and parameters.
    void additionalImports(EpetraExt::HDF5 &HDF5, std::string const &filename);
    void additionalExports(EpetraExt::HDF5 &HDF5, std::string const &filename);

    // Write the state of the ocean to traditional fortran out files fort.*
    // Use matlab plot-scripts for visualization
    void printLegacyFiles();

    // Initializer members
    void initializeOcean();
    void initializePreconditioner();
    void initializeBelos();

    // Perform a Newton solve with a small perturbation in the parameter
    Teuchos::RCP<Epetra_Vector> initialState();

    void inspectVector(VectorPtr x);
};
#endif
