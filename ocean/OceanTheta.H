#ifndef OCEANTHETA_H
#define OCEANTHETA_H

#include "Ocean.H"

//======================================================================
// class OceanTheta should have all the theta and time stuff and
// inherit the rest
//======================================================================

class OceanTheta: public Ocean
{
public:
	//constructor
	OceanTheta(Teuchos::RCP<Epetra_Comm> Comm);

	//destructor
	~OceanTheta(){ INFO("OceanTheta destructor called...");	}

	void Store();
	void Restore();

	void SetTimestep(double timestep){ timestep_ = timestep; }
	void SetTheta(double theta){ theta_ = theta; }
	
	void ComputeRHS();
	void ComputeJacobian();
	
protected:
	double theta_;
	double timestep_;
	
	Teuchos::RCP<Epetra_Vector> oldState_;
	Teuchos::RCP<Epetra_Vector> stateDot_;
	Teuchos::RCP<Epetra_Vector> oldRhs_;

};

#endif
