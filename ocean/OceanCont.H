#ifndef OCEANCONT_H
#define OCEANCONT_H

#include "Ocean.H"

//======================================================================
class OceanCont: public Ocean
{
	// THCM continuation parameters and bounds
	int    parIdent_;    // parameter identifier (integer)
	double parValue_;    // value of the parameter
	double parStart_;    // its starting value
	double parEnd_;	     // its ending value
	double storedPar1_;  // stored value of the parameter
	double storedPar2_;  // stored value of the parameter
		
public:
	//constructor
	OceanCont(Teuchos::RCP<Epetra_Comm> Comm);

	//destructor
	~OceanCont(){ INFO("OceanCont destructor called...");	}

	// The parameter get and set members wrap the corresponding
	// Fortran functions. Just as in THCM, these functions
	// may do more in the future,
	//--> not sure if getpardestination is needed
	void    setPar(double value);
	double  getPar();
	double  getParDestination() { return parEnd_; }
};

#endif
