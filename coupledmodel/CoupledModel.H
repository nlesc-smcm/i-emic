#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

#include "Ocean.H"
#include "Atmosphere.H"
#include "SuperVector.H"

#include <vector>
#include <memory>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

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

	// Trilinos-based parallel implicit ocean model (Trilinos-THCM)
	Teuchos::RCP<Ocean>  ocean_;
	
	// Plain std::vector-based atmosphere model
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
	char solvingScheme_;
	
	// Trilinos' MPI-like communicator
	Teuchos::RCP<Epetra_Comm> comm_;

	// Atmosphere -> Ocean block in Jacobian
	std::shared_ptr<std::vector<double> > B_;	
	std::shared_ptr<std::vector<int> > rowsB_;

	// Ocean -> Atmosphere block in Jacobian
	std::shared_ptr<std::vector<double> > C_;	
	std::shared_ptr<std::vector<int> > rowsC_;
	
public:
	CoupledModel(Teuchos::RCP<Epetra_Comm> comm,
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
	// for internal use in the inexact solver
	std::shared_ptr<SuperVector> getSolution(char mode1, char mode2);

	// Synchronize the states between the models that are needed to communicate
	// Relaxation might not be a good idea, set it to 1
	void synchronize(double relaxation = 1.0);

	// A few tests are collected here
	void test();
};

#endif
