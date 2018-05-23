#ifndef MODEL_H
#define MODEL_H

#include <Teuchos_RCP.hpp>
#include <TRIOS_Domain.H>
#include "Utils.H"

// forward declarations
class Epetra_Vector;

class Model
{
    
public:
    
    using VectorPtr = Teuchos::RCP<Epetra_Vector>;
    
    virtual ~Model() {}
    
    virtual VectorPtr getState(char mode) = 0;
    virtual VectorPtr getSolution(char mode) = 0;
    virtual VectorPtr getRHS(char mode) = 0;
    
    virtual void computeRHS() = 0;

    virtual void setParName(std::string const &parName) = 0;

    virtual Utils::MaskStruct getLandMask() = 0;
    virtual void setLandMask(Utils::MaskStruct const &mask) = 0;    

    virtual std::string const name() = 0;

    virtual Teuchos::RCP<TRIOS::Domain> getDomain() = 0;
    
    virtual std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Model> model) = 0;

    virtual int dof() = 0;
};

#endif
