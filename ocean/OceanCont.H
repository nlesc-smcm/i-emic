#ifndef OCEANCONT_H
#define OCEANCONT_H

#include "Ocean.H"

//======================================================================

class OceanCont: public Ocean
{
	Teuchos::RCP<Epetra_Vector> oldState_;
	Teuchos::RCP<Epetra_Vector> oldRhs_;

	// THCM continuation parameters and bounds
	int    thcmParIdent_;  // integer identifier
	double thcmPar_;       // value of the parameter
	double thcmParStart_;  // its starting value
	double thcmParEnd_;	   // its ending value
	
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
