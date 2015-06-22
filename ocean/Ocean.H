//!------------------------------------------------------------------------
/* Ocean is a model within the I-EMIC. It is a wrapper for Jonas Thies'
   Trilinos-THCM, using only a fraction of the functionality that is
   available in there.
   
   For now I try to let matrices and solvers be only available at the
   model-level, but I'm doubting whether this is the right approach... In
   the future I will probably need to assemble a big matrix and find an
   appropriate solver. 

   A model is required to have the following member functions:

       void computeJacobian()
	   void computeRHS() 
	   void solve(rhs)	   
	   vector getSolution()
	   vector getState()
	   vector getRHS()

   The vector type that is returned from the get{.} functions will have
   to be known by the interface (e.g. Continuation) that calls this class.
   The calling interface requires the following member functions:

       void update()
	   void scale()
	   int length()	   
	   double dot()
	   double norm()
*/
//!------------------------------------------------------------------------
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

class Ocean
{
protected:

	// 
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
	void solve(Teuchos::RCP<Vector> rhs = Teuchos::null);
	void computeRHS();
	void computeJacobian();
	
	// All get...() functions that return a VectorPtr type should accept
	// a character specifying the access the user has: Copy 'C' or View 'V'
	// By default this should be copy. 
	Teuchos::RCP<Vector> getSolution(char mode = 'C');
	Teuchos::RCP<Vector> getState(char mode = 'C');
	Teuchos::RCP<Vector> getRHS(char mode = 'C');
	
	void setState(Teuchos::RCP<Vector> state) { state_ = state->getRCP(); }
	void setRHS(Teuchos::RCP<Vector> rhs)     { rhs_   = rhs->getRCP();   }

	void saveStateToFile(std::string const &name);
	void loadStateFromFile(std::string const &name);

	// Write the state of the ocean to traditional FORTRAN out files fort.*
	// Use matlab plot-scripts for visualization 
	void dumpState();	
	
private:
	// Randomize state, scale with scaling
	void randomizeState(double scaling);
	void initializeSolver();
	void scaleProblem();
	void unscaleProblem();

	Teuchos::RCP<Vector> getVector(char mode, Teuchos::RCP<Epetra_Vector> vec);
};
#endif
