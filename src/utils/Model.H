#ifndef MODEL_H
#define MODEL_H

#include "Utils.H"

// forward declarations
namespace Teuchos { template<class T> class RCP; }

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
    
    using VectorPtr = Teuchos::RCP<Epetra_Vector>;
    
    virtual ~Model() {}
    
    virtual void computeRHS() = 0;
    virtual void computeJacobian() = 0;
    
    virtual void applyMatrix(Epetra_MultiVector const &v, Epetra_MultiVector &out) = 0;
    virtual void applyMassMat(Epetra_MultiVector const &v, Epetra_MultiVector &out) = 0;
    virtual void applyPrecon(Epetra_MultiVector const &v, Epetra_MultiVector &out) = 0;
    
    virtual Teuchos::RCP<Epetra_CrsMatrix> getJacobian() = 0;

    virtual VectorPtr getState(char mode) = 0;
    virtual VectorPtr getSolution(char mode) = 0;
    virtual VectorPtr getRHS(char mode) = 0;

    virtual double getPar(std::string const &parName) = 0;

    virtual void setParName(std::string const &parName) = 0;
    virtual void setPar(std::string const &parName, double value) = 0;

    virtual Utils::MaskStruct getLandMask() = 0;
    virtual void setLandMask(Utils::MaskStruct const &mask) = 0;    

    virtual std::string const name() = 0;

    virtual Teuchos::RCP<TRIOS::Domain> getDomain() = 0;

    // Model's own getBlock to distribute calls among submodels
    template <typename T>
    std::shared_ptr<Utils::CRSMat> getBlock(T model);

    // getBlock members to compute derivative w.r.t. any other model. 
    virtual std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Ocean> ocean) = 0;
    virtual std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Atmosphere> atmos) = 0;
    virtual std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<SeaIce> seaice) = 0;

    // Our own synchronize to distribute synchronizations among submodels
    template <typename T>
    void synchronize(T model);
    
    virtual void synchronize(std::shared_ptr<Ocean> ocean) = 0;
    virtual void synchronize(std::shared_ptr<Atmosphere> atmos) = 0;
    virtual void synchronize(std::shared_ptr<SeaIce> seaice) = 0;

    virtual int dof() = 0;
    
    virtual void buildPreconditioner() = 0;

    virtual void preProcess() = 0;
    virtual void postProcess() = 0;

    virtual std::string const writeData(bool describe) = 0;

    
};

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
