//!------------------------------------------------------------------------
/* Ocean is a model within the I-EMIC. It is a wrapper for Jonas Thies'
   Trilinos-THCM, using only a fraction of the functionality that is
   available in there.
   
   A model is required to have the following member functions:

	   void computeRHS() 
       void computeJacobian()
	   void solve(rhs)	   
	   	   
	   vector *getSolution()
	   vector *getRHS()
	   vector *getState()	   
	   double getPar()
	   void   setPar()
	   
	   void preProcess()
	   void postProcess()

   The vector that is returned from the get{.} functions will have
   to be known by the interface (e.g. Continuation) that calls this class.

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
#include "SuperVector.H"

// forward declarations
class THCM;

namespace TRIOS
{ class Domain; }

class Epetra_Comm;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Operator;
class Epetra_CrsMatrix;

class Ocean
{
public:
	using VectorPtr = Teuchos::RCP<SuperVector>;

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
	
	bool solverInitialized_;
	bool adaptivePrecCompute_;
	bool recomputePreconditioner_;
	double gmresIters_;
	bool useScaling_;
	
	Teuchos::RCP<Epetra_Vector> sol_;
	
	Teuchos::RCP<Teuchos::ParameterList> belosParamList_;
	Teuchos::RCP<Belos::LinearProblem
				 <double, Epetra_MultiVector, Epetra_Operator> > problem_;	
	Teuchos::RCP<Belos::BlockGmresSolMgr
				 <double, Epetra_MultiVector, Epetra_Operator> > belosSolver_;
	
	Teuchos::RCP<Ifpack_Preconditioner> precPtr_;
	
	int recomputeBound_;

	// Domain object
	Teuchos::RCP<TRIOS::Domain> domain_;

	// Problem dimensions:
	int N_, M_, L_;

 	// HDF5 input/output files
	std::string const inputFile_;
	std::string const outputFile_;

	// Flag for loading pre-existing state
	bool useExistingState_;
	
	//------------------------------------------------------------------
	// Datamembers in charge of the continuation parameter
	// --> only parValue should be in here...
	int    parIdent_;    // parameter identifier (integer)
	double parValue_;    // value of the parameter
	double parStart_;    // its starting value
	double parEnd_;	     // its ending value
	double storedPar1_;  // stored value of the parameter
	double storedPar2_;  // stored value of the parameter
	
public:

	//constructor
	Ocean(Teuchos::RCP<Epetra_Comm> Comm,
		  Teuchos::RCP<Teuchos::ParameterList> oceanParamList);

	//destructor
	~Ocean();

	// ----------------------------------------------------------------------+
	// The following functions are required by the continuation interface:
	//
	// Solve may optionally accept an rhs of VectorPointer type
	void solve(Teuchos::RCP<SuperVector> rhs = Teuchos::null);
	void computeRHS();
	void computeJacobian();
	
	// All get...() functions that return a VectorPtr type should accept
	// a character specifying the access the user has: Copy 'C' or View 'V'
	// By default this should be copy. 
	Teuchos::RCP<SuperVector> getSolution(char mode = 'C');
	Teuchos::RCP<SuperVector> getState(char mode = 'C');
	Teuchos::RCP<SuperVector> getRHS(char mode = 'C');
	
	void setState(Teuchos::RCP<SuperVector> state)
		{ state_ = state->getOceanVector(); }
	
	void setRHS(Teuchos::RCP<SuperVector> rhs)
		{ rhs_   = rhs->getOceanVector(); }

	// The parameter get and set members wrap the corresponding
	// Fortran functions. Just as in THCM, these functions
	// may do more in the future,
	//--> not sure if using getpardestination is ok from a design
	//    perspective. The continuation routine should have that
	//    information, not the model right?
	void    setPar(double value);

	// Return continuation parameter or related information.
	// Different modes   'V': return parameter value (default)
	//                   'D': return parameter destination value
	//                   'S': return parameter start value
	double  getPar(char mode = 'V');
	
	double  getParDestination() { return parEnd_; }

	// ! Pre and post-convergence processing
	void preProcess();
	void postProcess();
	
	// Get functions for the grid dimensions
	int getNdim() { return N_;}
	int getMdim() { return M_;}
	int getLdim() { return L_;}

	// ----------------------------------------------------------------------+
	// The following functions are required for the coupling with other
	// models:

	// Obtain an external atmosphere
	void setAtmosphere(std::vector<double> const &atmos);

	// Obtain the Jacobian from THCM
	Teuchos::RCP<Epetra_CrsMatrix> getJacobian();

    // Get the diagonal elements of the atmosphere to ocean
	// contribution in the Jacobian
	std::shared_ptr<std::vector<double> > getAtmosBlock();
	
	// Obtain ocean/land surface temperature
	std::shared_ptr<std::vector<double> > getSurfaceT();
	
	// Get the rows in the matrix associated with surface temperature
	std::shared_ptr<std::vector<int> > getSurfaceTRows();

	// Get the ocean's landmask;
	std::shared_ptr<std::vector<int> > getLandMask();
	
private:
	// HDF5-based save and load functions for the state and parameters
	void saveStateToFile(std::string const &name);
	void loadStateFromFile(std::string const &name);

	// Write the state of the ocean to traditional FORTRAN out files fort.*
	// Use matlab plot-scripts for visualization 
	void writeFortFiles();

	// Randomize state, scale with scaling
	void randomizeState(double scaling);

	void initializeSolver();

	// Scale the matrix and rhs of the problem
	void scaleProblem();
	void unscaleProblem();

	Teuchos::RCP<SuperVector> getVector(char mode, Teuchos::RCP<Epetra_Vector> vec);
};
#endif








