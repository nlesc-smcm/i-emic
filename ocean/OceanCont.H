#ifndef OCEANCONT_H
#define OCEANCONT_H

#include "Ocean.H"

//======================================================================
class OceanCont: public Ocean
{
	Teuchos::RCP<Epetra_Vector> storedState_;
	Teuchos::RCP<Epetra_Vector> storedRhs_;

	// THCM continuation parameters and bounds
	int    parIdent_;   // parameter identifier (integer)
	double parValue_;   // value of the parameter
	double parStart_;   // its starting value
	double parEnd_;	    // its ending value
	double storedPar_;  // stored value of the parameter
		
public:
	//constructor
	OceanCont(Teuchos::RCP<Epetra_Comm> Comm);

	//destructor
	~OceanCont(){ INFO("OceanCont destructor called...");	}

	// The parameter get and set members wrap the corresponding
	// Fortran functions. Just as in THCM, these functions
	// may do more in the future,
	void    SetPar(double value);
	double  GetPar();
	double  GetParDestination() { return parEnd_; }
	
	void Store();
	void Restore();
	
	// All Get...() functions that return a VectorPtr type should accept
	// a character, specifying the access the user has: Copy 'C' / View 'V'
	// By default this should be copy. 
	Teuchos::RCP<Vector> GetStoredState(char mode = 'C');
	Teuchos::RCP<Vector> GetStoredRHS(char mode = 'C');

	double  GetStoredPar() { return storedPar_; };

private:

	void StoreState();
	void RestoreState();
	void StoreRHS();
	void RestoreRHS();
	void StorePar(); 
	void RestorePar(); 
};

#endif
