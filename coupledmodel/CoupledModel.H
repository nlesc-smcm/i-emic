#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

#include "Ocean.H"
#include "Atmosphere.H"
#include "Vector.H"

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
	std::shared_ptr<Atmosphere> atmosphere_;

	// Combined state, right hand side and solution
	std::shared_ptr<Vector> stateView_;
	std::shared_ptr<Vector> rhsView_;
	std::shared_ptr<Vector> solView_;

	// number of models;
	int nm_; 

	// Trilinos' MPI-like communicator
	Teuchos::RCP<Epetra_Comm> comm_;
	
public:
	CoupledModel(Graph &couplings,
				 Teuchos::RCP<Epetra_Comm> comm);
	~CoupledModel(){}

	void computeJacobian();
	void computeRHS();
	void solve(std::shared_ptr<Vector> rhs = nullptr);

	std::shared_ptr<Vector> getSolution(char mode = 'C');
	std::shared_ptr<Vector> getState(char mode = 'C');
	std::shared_ptr<Vector> getRHS(char mode = 'C');

	void setState(std::shared_ptr<Vector> state);
	void setRHS(std::shared_ptr<Vector> rhs);

	void   setPar(double value);
	double getPar();
	double getParDestination();

	void dumpState();

private:
	void synchronize();
	void test();
	
};

#endif
