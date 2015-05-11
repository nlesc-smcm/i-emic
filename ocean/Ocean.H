#ifndef OCEAN_H
#define OCEAN_H

#include <Teuchos_RCP.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack_Preconditioner.h>

#include "GlobalDefinitions.H"
#include "Vector.H"

// forward declarations
class THCM;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Operator;
class Epetra_CrsMatrix;
class Epetra_Time;

//=======================================================================
// An abstract base class would be nice
//=======================================================================

class Ocean
{
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
	
public:

	//constructor
	Ocean(Teuchos::RCP<Epetra_Comm> Comm);
	//destructor
	~Ocean(){ INFO("Ocean destructor called..."); }

	// Solve may optionally accept an rhs of VectorPointer type
	void Solve(Teuchos::RCP<Vector> rhs = Teuchos::null);
	void ComputeRHS();
	void ComputeJacobian();

	double GetNormRHS();
	double GetNormState();

	// All Get...() functions that return a VectorPtr type should accept
	// a character specifying the access the user has: Copy 'C' or View 'V'
	// By default this should be copy. 
	Teuchos::RCP<Vector> GetSolution(char mode = 'C');
	Teuchos::RCP<Vector> GetState(char mode = 'C');
	Teuchos::RCP<Vector> GetRHS(char mode = 'C');
	
	void SetState(Teuchos::RCP<Vector> state) { state_ = state->GetRCP(); }
	void SetRHS(Teuchos::RCP<Vector> rhs)     { rhs_   = rhs->GetRCP();   }

	// --> Trying to prevent thetastepper from breaking
	void SetState(Teuchos::RCP<Epetra_Vector> state) { state_ = state; }
	void SetRHS(Teuchos::RCP<Epetra_Vector> rhs)     { rhs_   = rhs;   }

	void SaveStateToFile(std::string const &name);
	void LoadStateFromFile(std::string const &name);

	// Write the state of the ocean to traditional FORTRAN out files fort.*
	// Use matlab plot-scripts for visualization 
	void DumpState();	
	
private:

	// Randomize state, scale with scaling
	void RandomizeState(double scaling);
	void InitializeSolver();
	void ScaleProblem();
	void UnscaleProblem();		
};
#endif
