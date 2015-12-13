//!---------------------------------------------------------------------
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
//!---------------------------------------------------------------------
#ifndef OCEAN_H
#define OCEAN_H

#include <Teuchos_RCP.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack_Preconditioner.h>

#include "GlobalDefinitions.H"
#include "SuperVector.H"
#include "IDRSolver.H"

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

	// 'F' Flexible GMRES (Belos/Trilinos)
	// 'I' IDR 
	char solverType_;
	
	bool   solverInitialized_;
	bool   precInitialized_;
	bool   recompPreconditioner_;
	
	Teuchos::RCP<Epetra_Vector> sol_;

	Teuchos::RCP<Teuchos::ParameterList> solverParams_;
	
	// Belos flexible GMRES members
	Teuchos::RCP<Teuchos::ParameterList> belosParamList_;
	Teuchos::RCP<Belos::LinearProblem
				 <double, Epetra_MultiVector, Epetra_Operator> > problem_;	
	Teuchos::RCP<Belos::BlockGmresSolMgr
				 <double, Epetra_MultiVector, Epetra_Operator> > belosSolver_;

	// IDR solver
	IDRSolver<Ocean, Teuchos::RCP<SuperVector> > idrSolver_;

	Teuchos::RCP<Ifpack_Preconditioner> precPtr_;
	
	// Domain object
	Teuchos::RCP<TRIOS::Domain> domain_;

	// Problem dimensions:
	int N_, M_, L_;

 	// HDF5 input/output files
	std::string const inputFile_;
	std::string const outputFile_;

	// Flags for loading/saving states from/to hfd5 file
	bool loadState_;
	bool saveState_;
	
	//------------------------------------------------------------------
	// Continuation parameter name.
	// For a list of possible names see THCM::par2int()
	std::string parName_;
	
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
	void solve(VectorPtr rhs = Teuchos::null);

	// Clear the search space of the IDR solver
	void clearSearchSpace() { idrSolver_.clearSearchSpace();}

	// Calculate explicit residual norm
	double explicitResNorm(VectorPtr rhs);
	
	void computeRHS();
	void computeJacobian();
	
	// All get...() functions that return a VectorPtr type should accept
	// a character specifying the access the user has: Copy 'C' or View 'V'
	// By default this should be copy. 
	Teuchos::RCP<SuperVector> getSolution(char mode = 'C');
	Teuchos::RCP<SuperVector> getState(char mode = 'C');
	Teuchos::RCP<SuperVector> getRHS(char mode = 'C');

	// Apply the Jacobian to a vector and return the result u = J*v
	Teuchos::RCP<SuperVector> applyMatrix(SuperVector const &v);
	void applyMatrix(SuperVector const &v, SuperVector &out);

	// Apply the preconditioner inverse to a vector and return
	// the result u = P^{-1}*v
	Teuchos::RCP<SuperVector> applyPrecon(SuperVector const &v);
	void applyPrecon(SuperVector const &v, SuperVector &out);

	void recomputePreconditioner() { recompPreconditioner_ = true; }

	// // THESE ARE OBSOLETE!!! there is a minor dependency in the continuation
	// storage stuff. Since we save to HDF5 we can get rid of that
	void setState(Teuchos::RCP<SuperVector> state)
		{ state_ = state->getOceanVector(); }
	
	void setRHS(Teuchos::RCP<SuperVector> rhs)
		{ rhs_   = rhs->getOceanVector(); }

	// The parameter set members wrap the corresponding
	// Fortran functions. We keep an internal parameter name which
	// can be overruled here.
	void setPar(double value);
	void setPar(std::string parName, double value);

	// Return continuation parameter value
	double getPar();
	double getPar(std::string const &parName);

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
	
	// Scale the matrix and rhs of the problem
	void scaleProblem(VectorPtr rhs = Teuchos::null);
	void unscaleProblem(VectorPtr rhs = Teuchos::null);

	int loadState() { return loadStateFromFile(inputFile_); }
	int saveState() { return saveStateToFile(outputFile_); }
	
private:
	// HDF5-based save and load functions for the state and parameters
	int loadStateFromFile(std::string const &name);
	int saveStateToFile(std::string const &name);

	// Write the state of the ocean to traditional FORTRAN out files fort.*
	// Use matlab plot-scripts for visualization 
	void writeFortFiles();

	// Initializer members
	void initializeOcean();
	void initializeSolver();
	void initializePreconditioner();
	void initializeBelos();
	void initializeIDR();

	// Build preconditioner
	void buildPreconditioner();

	Teuchos::RCP<SuperVector> getVector(char mode, Teuchos::RCP<Epetra_Vector> vec);
};
#endif








