//!
//! Ocean is a model within the I-EMIC. It is a wrapper for Jonas Thies'
//!   Trilinos-THCM, using only a fraction of the functionality that is
//!   available in there. 
//!   
//!   A model is required to have the following member functions:
//!
//!     VectorPtr getState()
//!     VectorPtr getRHS()
//!     VectorPtr getSolution()
//!
//!     double getPar()
//!     void   setPar()
//!
//!     void preProcess()
//!     void postProcess()
//!
//!     void computeRHS()
//!     void computeJacobian()
//!     void solve()
//!   
//!   The VectorPtr type that is returned from the get{.} functions will have
//!   to be known by the interface (e.g. Continuation) that calls this class.


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

#include <string>

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
	// We need both a distributed and a global version of the land mask, so
	// we use a struct.
	struct maskStruct
	{
		Teuchos::RCP<Epetra_IntVector> local;
		std::shared_ptr<std::vector<int> > global;
		std::string label;
	};
	
public:
	
	using VectorPtr = Teuchos::RCP<SuperVector>;
	using MatrixPtr = Teuchos::RCP<Epetra_CrsMatrix>;
	using PreconPtr = Teuchos::RCP<Ifpack_Preconditioner>;
	using LandMask  = maskStruct;

protected:
	
	// 
	Teuchos::RCP<THCM> thcm_;
	
	Teuchos::RCP<Epetra_Vector> state_;
	Teuchos::RCP<Epetra_Vector> rhs_;
	Teuchos::RCP<Epetra_Vector> massMatrix_;	
	
	Teuchos::RCP<Epetra_CrsMatrix> jac_;
	Teuchos::RCP<Epetra_Comm> comm_;
	
	Teuchos::RCP<Epetra_Vector> rowScaling_;
	Teuchos::RCP<Epetra_Vector> colScaling_;
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

	// Gathererd landmask, including borders and lids
	std::shared_ptr<std::vector<int> > landmask_;

	// Gathered surface mask, without borders
	std::shared_ptr<std::vector<int> > surfmask_;

	// Label of current mask
	std::string currentMask_;

 	// HDF5 input/output files
	std::string const inputFile_;
	std::string const outputFile_;

	// Flags for loading/saving states from/to hfd5 file
	bool loadState_;
	bool saveState_;
	bool saveMask_;

	// Enable saving every state in multiple formats
	bool storeEverything_;
	
	// Continuation parameter name.
	// For a list of possible names see THCM::par2int()
	std::string parName_;

	// Landmask file
	std::string landmaskFile_;

	// Row map for pressure points P
	Teuchos::RCP<Epetra_Map> mapP_;

	// Row map for velocity U
	Teuchos::RCP<Epetra_Map> mapU_;

	// Rows corresponding to singular pressure points
	Teuchos::RCP<Epetra_Vector> singRows_;

public:

	//! constructor
	Ocean(Teuchos::RCP<Epetra_Comm> Comm,
		  Teuchos::RCP<Teuchos::ParameterList> oceanParamList);

	//! destructor
	~Ocean();

	//! Solve may optionally accept an rhs of VectorPointer type
	void solve(VectorPtr rhs = Teuchos::null);

	//! Calculate explicit residual norm
	double explicitResNorm(VectorPtr rhs);
	void printResidual(VectorPtr rhs);
	
	void computeRHS();
	void computeJacobian();
	
	//! All get...() functions that return a VectorPtr type should accept
	//! a character specifying the access the user has: Copy 'C' or View 'V'
	//! By default this should be copy. 
	Teuchos::RCP<SuperVector> getSolution(char mode = 'C');
	Teuchos::RCP<SuperVector> getState(char mode = 'C');
	Teuchos::RCP<SuperVector> getRHS(char mode = 'C');

	//! Return pointer to Jacobian
	Teuchos::RCP<Epetra_CrsMatrix> getJacobian() {return jac_;}
	
	//! Apply the Jacobian to a vector and return the result u = J*v
	Teuchos::RCP<SuperVector> applyMatrix(SuperVector const &v);
	void applyMatrix(SuperVector const &v, SuperVector &out);
	
	//! Apply the preconditioner inverse to a vector and return
	//! the result u = P^{-1}*v
	Teuchos::RCP<SuperVector> applyPrecon(SuperVector const &v);
	void applyPrecon(SuperVector const &v, SuperVector &out);

	//! Set prec recompute flag
	void recomputePreconditioner() { recompPreconditioner_ = true; }

 	//! Build preconditioner
	void buildPreconditioner(bool forceInit = false);

	//! Initialize solver
	void initializeSolver();

	//! Get pointer to preconditioning operator
	PreconPtr getPreconPtr() { return precPtr_; }

	//! The parameter set members wrap the corresponding
	//! Fortran functions. We keep an internal parameter name which
	//! can be overwritten here.
	void setPar(double value);
	void setPar(std::string const &parName, double value);
	void setParName(std::string const &parName) { parName_ = parName; }

	// Return continuation parameter value
	double getPar();
	double getPar(std::string const &parName);

	//! Get global and local (distributed) landmask arrays.
	//! This routine also tests the obtained landmask and it will
	//! therefore also be set in THCM.
	LandMask getLandMask(std::string const & fname);

	//! Check the current mask
	bool isCurrentMask(std::string const &label)
		{
			return (currentMask_ == label);
		}

	//! Set global and local (distributed) landmask arrays
	void setLandMask(LandMask mask, bool global = false);

	//! Apply landmask 
	void applyLandMask(LandMask mask, double factor);	
	void applyLandMask(Teuchos::RCP<Epetra_Vector> x,
					   LandMask mask, double factor);	
	void applyLandMask(Teuchos::RCP<Epetra_Vector> x,
					   LandMask maskA, LandMask maskB);

	//! Additional monitor for continuation
	bool monitor(){return false;}
		
	//! Pre and post-convergence processing
	void preProcess();
	void postProcess();
	
	//! Get functions for the grid dimensions
	int getNdim() { return N_;}
	int getMdim() { return M_;}
	int getLdim() { return L_;}

	// ---------------------------------------------------------------
	// The following functions are required for the coupling with other
	// models:

	// Obtain an external atmosphere
	void setAtmosphere(std::vector<double> const &atmos);

    // Get the atmosphere to ocean
	// contribution in the Jacobian, see forcing.F90.
	void getAtmosBlock(std::vector<double> &values,
					   std::vector<int> &row_inds);
	
	std::shared_ptr<std::vector<double> >
	getLandTemperature(std::shared_ptr<std::vector<double> > tatm);
	
	// Obtain ocean/land surface temperature
	std::shared_ptr<std::vector<double> > getSurfaceT();
	
	// Get the rows in the matrix associated with surface temperature
	std::shared_ptr<std::vector<int> > getSurfaceTRows();

	// Get the ocean's current landmask (gathered on all procs);
	std::shared_ptr<std::vector<int> > getLandMask()    { return landmask_; }

	// Get the ocean's current surface landmask (gathered on all procs);
	std::shared_ptr<std::vector<int> > getSurfaceMask() { return surfmask_; }
	
	// ---------------------------------------------------------------
	// Scale the matrix and rhs of the problem
	void scaleProblem(VectorPtr rhs = Teuchos::null);
	void unscaleProblem(VectorPtr rhs = Teuchos::null);
	Teuchos::RCP<Epetra_Vector> getRowScaling();
	Teuchos::RCP<Epetra_Vector> getColScaling();

	int loadState() { return loadStateFromFile(inputFile_); }
	int saveState() { return saveStateToFile(outputFile_); }
 
	void setParameters(Teuchos::RCP<Teuchos::ParameterList> pars);
	
 	void copyFiles();
 	void copyFiles(std::string const &fname);
	void copyMask(std::string const &fname);
private:
	// HDF5-based save and load functions for the state and parameters
	int loadStateFromFile(std::string const &name);
	int saveStateToFile(std::string const &name);

	// Write the state of the ocean to traditional fortran out files fort.*
	// Use matlab plot-scripts for visualization 
	void printFiles();

	// Initializer members
	void initializeOcean();
	void initializePreconditioner();
	void initializeBelos();
	void initializeIDR();

	// This function finds the pressure rows that will give a
	// singular Schur complement.
	// Based on these rows we will try to create a better landmask.
	int analyzeJacobian();

	
	Teuchos::RCP<SuperVector>
	getVector(char mode, Teuchos::RCP<Epetra_Vector> vec);

	void inspectVector(Teuchos::RCP<Epetra_Vector> x);
};
#endif
