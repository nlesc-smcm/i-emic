#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

#include "Ocean.H"
#include "Atmosphere.H"
#include "SuperVector.H"

#include <vector>
#include <memory>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

// This class assembles multiple sub-models into a CoupledModel. Which models and couplings are used should be supplied in the accompanying xml file (coupledmodel_params.xml). For the state we designed a special container called a SuperVector, see it's class header for more documentation.

// Most routines in CoupledModel simply call their corresponding routines in the sub-models. For instance, a call to computeJacobian() translates to a call to computeJacobian() in the sub-models.

// The frequency at which the states of the sub-models are synchronized depends on the solving scheme. Moreover we use hashes to make sure nothing is computed and synchronized too often

// The Jacobians in the sub-models form the diagonal blocks of the CoupledModel's Jacobian. If we require each sub-model to maintain its own solve() routine we can use a block Gauss Seidel iteration to solve the coupled problem. This circumvents the assembly of a big matrix. Another option is to solve the problem decoupled and only synchronize when postProcess() is called.

// Forward declarations
class Epetra_Comm;

class CoupledModel
{
public:
	using VectorPtr = std::shared_ptr<SuperVector>;
	
private:
	// Adjacency graph describing the couplings between the models
	// NOT IMPLEMENTED YET
    // Graph couplings_;

	// Trilinos-based parallel implicit ocean model (Trilinos-THCM) with a
	// Belos linear solver performing a parallel FGMRES 
	Teuchos::RCP<Ocean>  ocean_;
	
	// STL std::vector-based implicit atmosphere model with a banded Lapack solver
	std::shared_ptr<Atmosphere> atmos_;

	// Combined state, right hand side and solution
	std::shared_ptr<SuperVector> stateView_;
	std::shared_ptr<SuperVector> rhsView_;
	std::shared_ptr<SuperVector> solView_;

	// Use existing (hdf5) state/par
	bool useExistingState_;
	
	// Solving scheme
	//   'D': decoupled
	//   'G': backward block Gauss Seidel
	char   solvingScheme_;

	int    iterGS_;      // Max iterations in GS scheme
	double toleranceGS_; // Tolerance in GS scheme

	// Use hashing to decide to act based on whether a state has changed
	// Hashes should be remembered for each task they get used for.
	bool   useHash_;  // flag to turn hashing on/off
	std::size_t syncHash_; // hash to base synchronize() decision on
	std::size_t rhsHash_;  // hash to base computeRHS() decision on	
	std::size_t jacHash_;  // hash to base computeJacobian() decision on

	IDRSolver<CoupledModel, std::shared_ptr<SuperVector> > idrSolver_;
	bool idrInitialized_;
	int idrSolveCtr_;
	int clearSPFreq_;
	
	// Trilinos' MPI-like communicator
	Teuchos::RCP<Epetra_Comm> comm_;

	// Atmosphere -> Ocean block in Jacobian
	std::shared_ptr<std::vector<double> > B_;	
	std::shared_ptr<std::vector<int> > rowsB_;

	// Ocean -> Atmosphere block in Jacobian
	std::shared_ptr<std::vector<double> > C_;	
	std::shared_ptr<std::vector<int> > rowsC_;
	
public:
	CoupledModel(Teuchos::RCP<Ocean> ocean,
				 std::shared_ptr<Atmosphere> atmos,
		         Teuchos::RCP<Teuchos::ParameterList> params);
	
	~CoupledModel(){}

	void computeJacobian();
	void computeRHS();

	// Solve the models
	void solve(std::shared_ptr<SuperVector> rhs = std::shared_ptr<SuperVector>());

	void initializeIDR();

	void applyMatrix(SuperVector const &v, SuperVector &out);
	void applyPrecon(SuperVector const &v, SuperVector &out);
	
	std::shared_ptr<SuperVector> applyMatrix(SuperVector const &v);
	std::shared_ptr<SuperVector> applyPrecon(SuperVector const &v);

	std::shared_ptr<SuperVector> getSolution(char mode = 'C');
	std::shared_ptr<SuperVector> getState(char mode = 'C');
	std::shared_ptr<SuperVector> getRHS(char mode = 'C');

	void setState(std::shared_ptr<SuperVector> state);
	void setRHS(std::shared_ptr<SuperVector> rhs);

	void   setPar(double value);
	double getPar();
	double getParDestination();

	// A controlling scheme may allow some pre-processing, for instance at the
	// start of a Newton process.
	void preProcess();
	
	// Similarly we can supply some post-processing, for instance when
	// a Newton process has converged.
	void postProcess();

private:
	// Solve the system with a few symmetric block GS iterations
	void blockGSSolve(std::shared_ptr<SuperVector> rhs);
	
	// For internal use: both modes for the underlying
	// getSolution methods can be called here.
	std::shared_ptr<SuperVector> getSolution(char mode1, char mode2);

	// Compute the residual ||b-Jx||
	double computeResidual(std::shared_ptr<SuperVector> rhs);
	
	// Synchronize the states between the models that are needed to communicate
	void synchronize();

	// Synchronize the landmask between the ocean and the atmosphere
	void synchronizeLandmask();

	// Create a hash of the extended state + parameter in an XOR / rotating manner
	std::size_t getHash();

	// Output, same as in Atmosphere but templated
	template<typename T>
	void write(std::vector<T> &vector, const std::string &filename);

	// print profile
	void printProfile(ProfileType profile);
	
	// A few tests are collected here
	void test();
};

#endif
