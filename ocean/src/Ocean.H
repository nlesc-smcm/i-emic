#ifndef OCEAN_H
#define OCEAN_H

#include <Teuchos_RCP.hpp>

#include "THCM.H"

#include "GlobalDefinitions.H"

// For now this is a dummy class

class Ocean
{
	Teuchos::RCP<THCM> thcm_;
public:
	//constructor
	Ocean(Teuchos::RCP<Epetra_Comm> Comm);
	//destructor
	~Ocean();
	void doStuff();
};
#endif
