#ifndef MODEL_H
#define MODEL_H

#include <fstream>
#include "Utils.H"
#include "TRIOS_Domain.H"

// forward declarations
// namespace Teuchos { template<class T> class RCP; }

namespace TRIOS { class Domain; }

class Epetra_MultiVector;
class Epetra_CrsMatrix;
class Epetra_Vector;

class Ocean;
class Atmosphere;
class SeaIce;

class Model
{

public:

    using VectorPtr      = Teuchos::RCP<Epetra_Vector>;
    using ConstVectorPtr = Teuchos::RCP<const Epetra_Vector>;
    using MatrixPtr      = Teuchos::RCP<Epetra_CrsMatrix>;
    using ParameterList  = Teuchos::RCP<Teuchos::ParameterList>;

    Model(Teuchos::RCP<Epetra_Comm> Comm) : comm_(Comm)
    {}

    //! non-overlapping state vector
    Teuchos::RCP<Epetra_Vector> state_;

    //! Epetra communicator
    Teuchos::RCP<Epetra_Comm> comm_;

    //! input/output flags
    bool loadState_, saveState_, loadMask_, saveMask_;

    //! save/copy frequency
    int saveEvery_;

    //! postprocessing counter
    int ppCtr_;

    //! HDF5 input/output filenames
    std::string inputFile_;
    std::string outputFile_;

    virtual ~Model() {}

    //! compute rhs (spatial discretization)
    virtual void computeRHS() = 0;

    //! compute derivative of rhs
    virtual void computeJacobian() = 0;

    //! compute mass matrix
    virtual void computeMassMat() = 0;

    virtual void applyMatrix(Epetra_MultiVector const &v, Epetra_MultiVector &out) = 0;
    virtual void applyMassMat(Epetra_MultiVector const &v, Epetra_MultiVector &out) = 0;
    virtual void applyPrecon(Epetra_MultiVector const &v, Epetra_MultiVector &out) = 0;

    virtual MatrixPtr getJacobian() = 0;

    virtual VectorPtr getState(char mode) = 0;
    virtual VectorPtr getSolution(char mode) = 0;
    virtual VectorPtr getRHS(char mode) = 0;

    virtual double getPar(std::string const &parName) = 0;
    virtual void setPar(std::string const &parName, double value) = 0;

    //! Return total number of continuation parameters in model
    virtual int npar() = 0;

    //! Convert integer parameter index to parameter name
    virtual std::string int2par(int ind) const = 0;

    virtual Utils::MaskStruct getLandMask() = 0;
    virtual void setLandMask(Utils::MaskStruct const &mask) = 0;

    virtual std::string name() const = 0;
    virtual int modelIdent() const = 0;

    virtual Teuchos::RCP<TRIOS::Domain> getDomain() = 0;

    //! Model's own getBlock to distribute calls among submodels
    template <typename T>
    std::shared_ptr<Utils::CRSMat> getBlock(T model);

    //! getBlock members to compute derivative w.r.t. any other model.
    virtual std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Ocean> ocean)      = 0;
    virtual std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Atmosphere> atmos) = 0;
    virtual std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<SeaIce> seaice)    = 0;

    //! Our own synchronize to distribute synchronizations among submodels
    template <typename T>
    void synchronize(T model);

    //! Derived models need to implement synchronizations with any
    //! other model, including a dummy implementation of synchronize
    //! with itself.
    virtual void synchronize(std::shared_ptr<Ocean> ocean)      = 0;
    virtual void synchronize(std::shared_ptr<Atmosphere> atmos) = 0;
    virtual void synchronize(std::shared_ptr<SeaIce> seaice)    = 0;

    //! degrees of freedom (excluding any auxiliary unknowns)
    virtual int dof() = 0;

    virtual void buildPreconditioner() = 0;

    virtual void preProcess()  = 0;

    virtual void postProcess() = 0;

    //! Plaintext data output
    virtual std::string writeData(bool describe) const = 0;

    //! HDF5-based load function for the state and parameters
    int loadStateFromFile(std::string const &filename);

    //! Additional, model-specific queries for the HDF5 object
    virtual void additionalImports(EpetraExt::HDF5 &HDF5,
                                   std::string const &filename) = 0;

    //! HDF5-based save function for the state and parameters
    int saveStateToFile(std::string const &filename);

    //! Copy outputFile_ to <prepend>outputFile_
    int copyState(std::string const &prepend);

    //! Additional, model-specific writes for the HDF5 object
    virtual void additionalExports(EpetraExt::HDF5 &HDF5,
                                   std::string const &filename) = 0;

    //! Convert global id to coordinates i,j,k,xx and model identification mdl
    void gid2coord(int const &gid, int &mdl,
                   int &i, int &j, int &k, int &xx);

    virtual void initializeState() { state_->PutScalar(0.0); }

    virtual void pressureProjection(VectorPtr vec){}

};

// Implementations
//=============================================================================
inline int Model::loadStateFromFile(std::string const &filename)
{
    INFO("_________________________________________________________");
    if (loadState_)
    {
        INFO("Loading state and parameters from " << filename);
    }
    else
    {
        INFO("Performing only model specific import operations from " << filename);
    }

    // Check whether file exists
    std::ifstream file(filename);
    if (!file)
    {
        WARNING("Can't open " << filename
                << ", continue with trivial state", __FILE__, __LINE__);

        // create trivial state
        state_->PutScalar(0.0);
        return 1;
    }
    else file.close();

    // Create HDF5 object
    EpetraExt::HDF5 HDF5(*comm_);
    Epetra_MultiVector *readState;

    // Open file
    HDF5.Open(filename);

    if (loadState_)
    {
        // Check contents
        if (!HDF5.IsContained("State"))
        {
            ERROR("The group <State> is not contained in hdf5 " << filename,
                  __FILE__, __LINE__);
        }

        // Read the state. To be able to restart with different
        // numbers of procs we do not include the Map in the hdf5. The
        // state as read here will have a linear map and we import it
        // into the current domain decomposition.
        HDF5.Read("State", readState);

        if ( readState->GlobalLength() != getDomain()->GetSolveMap()->NumGlobalElements() )
        {
            WARNING("Loading state from differ #procs", __FILE__, __LINE__);
        }

        // Create importer
        // target map: domain StandardMap
        // source map: state with linear map as read by HDF5.Read
        Teuchos::RCP<Epetra_Import> lin2solve =
            Teuchos::rcp(new Epetra_Import(*(getDomain()->GetSolveMap()),
                                           readState->Map() ));

        // Import state from HDF5 into state_ datamember
        CHECK_ZERO(state_->Import(*((*readState)(0)), *lin2solve, Insert));

        delete readState;

        INFO(" state: ||x|| = " << Utils::norm(state_));

        // Interface between HDF5 and the parameters,
        // put all the <npar> parameters back in atmos.
        std::string parName;
        double parValue;

        // Check contents
        if (!HDF5.IsContained("Parameters"))
        {
            ERROR("The group <Parameters> is not contained in hdf5 " << filename,
                  __FILE__, __LINE__);
        }

        for (int par = 0; par < npar(); ++par)
        {
            parName  = int2par(par);

            // Read continuation parameter and set them in model
            try
            {
                HDF5.Read("Parameters", parName.c_str(), parValue);
            }
            catch (EpetraExt::Exception &e)
            {
                e.Print();
                continue;
            }

            setPar(parName, parValue);
            INFO("   " << parName << " = " << parValue);
        }
    }

    additionalImports(HDF5, filename);

    INFO("_________________________________________________________");
    return 0;
}

//=============================================================================
inline int Model::saveStateToFile(std::string const &filename)
{
    INFO("_________________________________________________________");
    INFO("Create backup of " << outputFile_);
    copyState(".bak");

    INFO("Writing state and parameters to " << filename);
    INFO("   state: ||x|| = " << Utils::norm(state_));

    // Write state, map and continuation parameter
    EpetraExt::HDF5 HDF5(*comm_);
    HDF5.Create(filename);
    HDF5.Write("State", *state_);

    // Interface between HDF5 and the parameters,
    // store all the <npar> parameters in an HDF5 file.
    std::string parName;
    double parValue;
    for (int par = 0; par < npar(); ++par)
    {
        parName  = int2par(par);
        parValue = getPar(parName);
        INFO("   " << parName << " = " << parValue);
        HDF5.Write("Parameters", parName.c_str(), parValue);
    }

    // Write grid information available in domain object
    HDF5.Write("Grid", "n",   getDomain()->GlobalN());
    HDF5.Write("Grid", "m",   getDomain()->GlobalM());
    HDF5.Write("Grid", "l",   getDomain()->GlobalL());
    HDF5.Write("Grid", "nun", getDomain()->Dof());
    HDF5.Write("Grid", "aux", getDomain()->Aux());

    HDF5.Write("Grid", "xmin", getDomain()->Xmin());
    HDF5.Write("Grid", "xmax", getDomain()->Xmax());
    HDF5.Write("Grid", "ymin", getDomain()->Ymin());
    HDF5.Write("Grid", "ymax", getDomain()->Ymax());
    HDF5.Write("Grid", "hdim", getDomain()->Hdim());


    // Write grid arrays
    std::string gridArrays[6] = {"x", "y", "z",
                                 "xu", "yv", "zw"};

    for (int i = 0; i != 6; ++i)
    {
        std::vector<double> array = (*getDomain()->GetGlobalGrid())[i];
        HDF5.Write("Grid", gridArrays[i], H5T_NATIVE_DOUBLE,
                   array.size(), &array[0]);
    }

    additionalExports(HDF5, filename);
    comm_->Barrier();

    INFO("_________________________________________________________");
    return 0;
}

//=============================================================================
inline int Model::copyState(std::string const &append)
{
    if (comm_->MyPID() == 0)
    {
        if (saveState_)
        {
            std::stringstream ss;
            ss << outputFile_ << append;
            INFO("copying " << outputFile_ << " to " << ss.str());
            std::ifstream src(outputFile_.c_str(), std::ios::binary);
            std::ofstream dst(ss.str(), std::ios::binary);
            dst << src.rdbuf();
        }
        else
        {
            WARNING("No use in copying a state when saveState = false",
                    __FILE__, __LINE__);
        }
    }
    return 0;
}

//=============================================================================
inline void Model::gid2coord(int const &gid, int &mdl,
                   int &i, int &j, int &k, int &xx)
{
            mdl = modelIdent();

            int N   = getDomain()->GlobalN();
            int M   = getDomain()->GlobalM();
            int L   = getDomain()->GlobalL();
            int dof = getDomain()->Dof();
            int dim = N*M*L*dof;

            int aux = gid - dim;
            if (aux >= 0) // this is an auxiliary unknown
            {
                xx = dof + aux;
                i  = 0;
                j  = 0;
                k  = 0;
            }
            else
            {
                int tmp = gid;

                xx  = tmp % dof;
                tmp = (tmp - xx) / dof;
                i   = tmp % N;
                tmp = (tmp - i) / N;
                j   = tmp % M;
                k   = (tmp - j) / M;
            }
        }

//=============================================================================
template <typename T>
std::shared_ptr<Utils::CRSMat> Model::getBlock(T model)
{
    auto ocean  = std::dynamic_pointer_cast<Ocean>(model);
    auto atmos  = std::dynamic_pointer_cast<Atmosphere>(model);
    auto seaice = std::dynamic_pointer_cast<SeaIce>(model);

    if (ocean)
        return this->getBlock(ocean);
    else if (atmos)
        return this->getBlock(atmos);
    else if (seaice)
        return this->getBlock(seaice);
    else
    {
        ERROR("Model: downcasting failed", __FILE__, __LINE__);
        return std::shared_ptr<Utils::CRSMat>();
    }
}

//=============================================================================
template <typename T>
void Model::synchronize(T model)
{
    auto ocean  = std::dynamic_pointer_cast<Ocean>(model);
    auto atmos  = std::dynamic_pointer_cast<Atmosphere>(model);
    auto seaice = std::dynamic_pointer_cast<SeaIce>(model);

    if (ocean)
        return this->synchronize(ocean);
    else if (atmos)
        return this->synchronize(atmos);
    else if (seaice)
        return this->synchronize(seaice);
    else
    {
        ERROR("Model: downcasting failed", __FILE__, __LINE__);
    }
}

#endif
