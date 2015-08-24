#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

#include "Ocean.H"
#include "Atmosphere.H"
#include "SuperVector.H"

#include <vector>
#include <memory>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

// This is a very important class and it should be documented well...

// We will use hashes to make sure nothing is computed too often

// typedef std::vector<std::vector<bool> > Graph;

// Forward declarations
class Epetra_Comm;

class CoupledModel
{
public:
	using VectorPtr = std::shared_ptr<SuperVector>;
	
private:
	// Adjacency graph describing the couplings between the models
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
	
	// Order of Neumann expansion in elimination based solve
	int kNeumann_;

	// Solving scheme
	//   'D': decoupled
    //   'E': elimination based solve using the Schur complement and one or
    //        more terms in the Neumann expansion of the Schur inverse
	//   'S': backward block SOR
	char solvingScheme_;

	// Use hashing to decide to act based on whether a state has changed
	// Hashes should be remembered for each task they get used for.
	bool   useHash_;  // flag to turn hashing on/off
	std::size_t syncHash_; // hash to base synchronize() decision on
	std::size_t rhsHash_;  // hash to base computeRHS() decision on	
	std::size_t jacHash_;  // hash to base computeJacobian() decision on

	int iterGS_; // Iterations in GS scheme
	
	double relaxSOR_;  // Relaxation in SOR scheme
	int iterSOR_;      // Iterations in SOR scheme
	
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

	std::shared_ptr<SuperVector> getSolution(char mode = 'C');
	std::shared_ptr<SuperVector> getState(char mode = 'C');
	std::shared_ptr<SuperVector> getRHS(char mode = 'C');

	void setState(std::shared_ptr<SuperVector> state);
	void setRHS(std::shared_ptr<SuperVector> rhs);

	void   setPar(double value);
	double getPar();
	double getParDestination();

	// If a step in a controlling scheme has converged some post-processing
	// is required, which is done here.
	void postConvergence();

private:
	// Solve the system with a few symmetric block GS iterations
	void blockGSSolve(std::shared_ptr<SuperVector> rhs);
	
	// Solve the coupled system in a few backward block SOR iterations.
	void blockSORSolve(std::shared_ptr<SuperVector> rhs);
	
	// Solve the system using an elimination based approach.
	void eliminationSolve(std::shared_ptr<SuperVector> rhs);
	
	// For internal use: both modes for the underlying
	// getSolution methods can be called here.
	std::shared_ptr<SuperVector> getSolution(char mode1, char mode2);
	
	// Synchronize the states between the models that are needed to communicate
	void synchronize();

	// Create a hash of the extended state + parameter in an XOR / rotating manner
	std::size_t getHash();

	// Output, same as in Atmosphere but templated^^
	template<typename T>
	void write(std::vector<T> &vector, const std::string &filename);

	// A few tests are collected here
	void test();
};

#endif
