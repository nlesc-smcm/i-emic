#ifndef OCEANCONT_H
#define OCEANCONT_H

#include "Ocean.H"

//======================================================================

class OceanCont: public Ocean
{
	Teuchos::RCP<Epetra_Vector> oldState_;
	Teuchos::RCP<Epetra_Vector> oldRhs_;

	// THCM continuation parameters and bounds
	int    parIdent_;  // integer identifier
	double par_;       // value of the parameter
	double parStart_;  // its starting value
	double parEnd_;	   // its ending value
	double oldPar_;    // stored value of the parameter
	
public:
	//constructor
	OceanCont(Teuchos::RCP<Epetra_Comm> Comm);

	//destructor
	~OceanCont(){ INFO("OceanCont destructor called...");	}

	// The parameter get and set members wrap the corresponding
	// Fortran functions. Just as in THCM, these functions
	// may do more in the future,
	double  GetPar();
	void    SetPar(double par);

	Teuchos::RCP<Epetra_Vector> GetSolutionCopy();

	void StoreState();
	void RestoreState();
	void StoreRHS();
	void RestoreRHS();

	Teuchos::RCP<Epetra_Vector> GetOldState() { return oldState_; }
	Teuchos::RCP<Epetra_Vector> GetOldRHS()   { return oldRhs_;   }	
};

#endif
