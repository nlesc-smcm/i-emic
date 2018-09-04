#ifndef OCEANTHETA_H
#define OCEANTHETA_H

#include "Ocean.H"

//======================================================================
// class OceanTheta has all the theta and time stuff and
// inherits the rest
//======================================================================

class OceanTheta: public Ocean
{
public:
	//constructor
	OceanTheta(Teuchos::RCP<Epetra_Comm> Comm,
			   Teuchos::RCP<Teuchos::ParameterList> params);

	//destructor
	~OceanTheta(){ INFO("OceanTheta destructor called...");	}

	void store();
	void restore();

	void setTimestep(double timestep){ timestep_ = timestep; }
	void setTheta(double theta){ theta_ = theta; }
	
	void computeRHS();
	void computeJacobian();
	
protected:
	double theta_;
	double timestep_;
	
	Teuchos::RCP<Epetra_Vector> oldState_;
	Teuchos::RCP<Epetra_Vector> stateDot_;
	Teuchos::RCP<Epetra_Vector> oldRhs_;

};

#endif
