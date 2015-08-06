#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

#include "Ocean.H"
#include "Atmosphere.H"
#include "SuperVector.H"

#include <vector>
#include <memory>

#include <Teuchos_RCP.hpp>

typedef std::vector<std::vector<bool> > Graph;

// Forward declarations
class Epetra_Comm;

class CoupledModel
{
	// Adjacency graph describing the couplings between the models
	Graph couplings_;

	// Trilinos-based parallel implicit ocean model (Trilinos-THCM)
	Teuchos::RCP<Ocean>  ocean_;
	
	// Plain std::vector-based atmosphere model
	std::shared_ptr<Atmosphere> atmos_;

	// Combined state, right hand side and solution
	std::shared_ptr<SuperVector> stateView_;
	std::shared_ptr<SuperVector> rhsView_;
	std::shared_ptr<SuperVector> solView_;

	// number of models;
	int nm_;

	// Order of Neumann expansion in elimination based solve
	int kNeumann_;

	// Trilinos' MPI-like communicator
	Teuchos::RCP<Epetra_Comm> comm_;

	// Atmosphere -> Ocean block in Jacobian
	std::shared_ptr<std::vector<double> > B_;	
	std::shared_ptr<std::vector<int> > rowsB_;

	// Ocean -> Atmosphere block in Jacobian
	std::shared_ptr<std::vector<double> > C_;	
	std::shared_ptr<std::vector<int> > rowsC_;
	
public:
	CoupledModel(Graph &couplings,
				 Teuchos::RCP<Epetra_Comm> comm);
	~CoupledModel(){}

	void computeJacobian();
	void computeRHS();

	// Solve the models
	//    mode = 'D': decoupled
    //    mode = 'S': coupled using the Schur complement and one or more terms
	//                in the Neumann expansion of the inverse
	void solve(std::shared_ptr<SuperVector> rhs = std::shared_ptr<SuperVector>(),
			   char mode = 'S');

	std::shared_ptr<SuperVector> getSolution(char mode = 'C');
	std::shared_ptr<SuperVector> getState(char mode = 'C');
	std::shared_ptr<SuperVector> getRHS(char mode = 'C');

	void setState(std::shared_ptr<SuperVector> state);
	void setRHS(std::shared_ptr<SuperVector> rhs);

	void   setPar(double value);
	double getPar();
	double getParDestination();

	void dumpState();

private:
	// for internal use in the inexact solver
	std::shared_ptr<SuperVector> getSolution(char mode1, char mode2);
	void synchronize(double relaxation = 1.0);
	void test();
};

#endif
