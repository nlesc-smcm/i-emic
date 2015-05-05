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
class Epetra_Time;

//==============================================================================
// An abstract base class would be nice
//==============================================================================

class Ocean
{
public:
	//constructor
	Ocean(Teuchos::RCP<Epetra_Comm> Comm);
	//destructor
	~Ocean(){ INFO("Ocean destructor called..."); }

	// Randomize state, scale with scaling
	void RandomizeState(double scaling);

	// Write the state of the ocean to traditional FORTRAN out files fort.*
	// Use matlab plot-scripts for visualization 
	void DumpState();

	// 
	void InitializeSolver();

	//
	void Solve(Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null);

	//
	void ScaleProblem();

	//
	void UnscaleProblem();	

	//
	void ComputeRHS();

	//
	void ComputeJacobian();

	//
	double GetNormRHS();

	//
	double GetNormState();

	//
	Teuchos::RCP<Epetra_Vector> GetState() { return state_; }

	//
	Teuchos::RCP<Epetra_Vector> GetSolution() { return sol_; }
	
	void SetState(Teuchos::RCP<Epetra_Vector> state) { state_ = state; }
	void RecomputePreconditioner() { recomputePreconditioner_ = true; }

	void SaveStateToFile(std::string const &name);
	void LoadStateFromFile(std::string const &name);
	
protected:
	Teuchos::RCP<THCM> thcm_;
	Teuchos::RCP<Epetra_Vector> state_;
	Teuchos::RCP<Epetra_Vector> rhs_;
	Teuchos::RCP<Epetra_Vector> massMatrix_;	
	Teuchos::RCP<Epetra_CrsMatrix> jac_;
	Teuchos::RCP<Epetra_Comm> comm_;
	
	Teuchos::RCP<Epetra_Vector> rowScalingRecipr_;
	Teuchos::RCP<Epetra_Vector> colScalingRecipr_;

	Teuchos::RCP<Epetra_Time> timer_;

	bool solverInitialized_;
	bool recomputePreconditioner_;
	bool useScaling_;

	Teuchos::RCP<Epetra_Vector> sol_;
	Teuchos::RCP<Teuchos::ParameterList> belosParamList_;
	Teuchos::RCP<Belos::LinearProblem
				 <double, Epetra_MultiVector, Epetra_Operator> > problem_;	
	Teuchos::RCP<Belos::BlockGmresSolMgr
				 <double, Epetra_MultiVector, Epetra_Operator> > belosSolver_;

	Teuchos::RCP<Ifpack_Preconditioner> precPtr_;

	int belosIters_;
	int recomputeBound_;
};

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
