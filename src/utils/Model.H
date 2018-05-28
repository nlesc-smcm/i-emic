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
class AtmospherePar;
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
    
    virtual std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Model> model) = 0;

    virtual int dof() = 0;

    virtual void synchronize(std::shared_ptr<Model> model) = 0;

    virtual void buildPreconditioner() = 0;

    virtual void preProcess() = 0;
    virtual void postProcess() = 0;

    virtual std::string const writeData(bool describe) = 0;

};

#endif
