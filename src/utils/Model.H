#ifndef MODEL_H
#define MODEL_H

#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Vector;

class Model
{
    
public:
    
    using VectorPtr = Teuchos::RCP<Epetra_Vector>;
    
    virtual VectorPtr getState(char mode) = 0;
    virtual VectorPtr getSolution(char mode) = 0;
    virtual VectorPtr getRHS(char mode) = 0;

    virtual void computeRHS() = 0;
    
    
};

#endif
