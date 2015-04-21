#ifndef OCEAN_H
#define OCEAN_H

#include <Teuchos_RCP.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack_Preconditioner.h>

#include "GlobalDefinitions.H"

// forward declarations
class THCM;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Operator;
class Epetra_CrsMatrix;

// Perhaps an abstract base class will be nice
class Ocean
{
public:
	//constructor
	Ocean(Teuchos::RCP<Epetra_Comm> Comm);
	//destructor
	~Ocean()
		{
			INFO("Ocean destructor called...");
		}

	// Randomize state, scale with scaling
	void randomizeState(double scaling);
	void dumpState();
	
	void initializeSolver();
	
	void solve();
	void unscaleProblem();	
	void scaleProblem();

	void computeRHS();
	void computeJacobian();

	double getNormRHS();
	double getNormState();
	
	Teuchos::RCP<Epetra_Vector> getState() { return state_; }
	Teuchos::RCP<Epetra_Vector> getSolution() { return sol_; }
	
	void setState(Teuchos::RCP<Epetra_Vector> state) { state_ = state; }
	
protected:
	Teuchos::RCP<THCM> thcm_;
	Teuchos::RCP<Epetra_Vector> state_;
	Teuchos::RCP<Epetra_Vector> rhs_;
	Teuchos::RCP<Epetra_Vector> massMatrix_;	
	Teuchos::RCP<Epetra_CrsMatrix> jac_;
	Teuchos::RCP<Epetra_Comm> comm_;
	
	Teuchos::RCP<Epetra_Vector> rowScalingRecipr_;
	Teuchos::RCP<Epetra_Vector> colScalingRecipr_;

	bool solverInitialized_;

	Teuchos::RCP<Epetra_Vector> sol_;
	Teuchos::RCP<Teuchos::ParameterList> belosParamList_;
	Teuchos::RCP<Belos::LinearProblem
				 <double, Epetra_MultiVector, Epetra_Operator> > problem_;	
	Teuchos::RCP<Belos::BlockGmresSolMgr
				 <double, Epetra_MultiVector, Epetra_Operator> > belosSolver_;

	Teuchos::RCP<Ifpack_Preconditioner> precPtr_;
};

// class OceanTheta should have all the theta stuff and
// inherit the rest
class OceanTheta: public Ocean
{
public:
	//constructor
	OceanTheta(Teuchos::RCP<Epetra_Comm> Comm);

	//destructor
	~OceanTheta()
		{
			INFO("OceanTheta destructor called...");
		}

	void parkModel();
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
