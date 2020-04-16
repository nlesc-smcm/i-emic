#ifndef MODEL_H
#define MODEL_H

#include "Utils.H"

// We include these here, since otherwise all classes that use
// VectorPtr and MatrixPtr will fail to compile.
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// forward declarations
// namespace Teuchos { template<class T> class RCP; }

namespace TRIOS { class Domain; }

class Epetra_MultiVector;

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

    virtual void initializeState();

    virtual void pressureProjection(VectorPtr vec){}

};

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
