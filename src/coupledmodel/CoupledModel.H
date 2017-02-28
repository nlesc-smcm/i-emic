#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

#include "Ocean.H"
#include "AtmospherePar.H"
#include "Combined_MultiVec.H"
#include "IDRSolver.H"
#include "GMRESSolver.H"

#include "CouplingBlock.H"

#include <vector>
#include <memory>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

//! This class assembles multiple sub-models into a CoupledModel. Which models and couplings are used should be supplied in the accompanying xml file (coupledmodel_params.xml).

//! Most routines in CoupledModel simply call their corresponding routines in the sub-models. For instance, a call to computeJacobian() translates to a call to computeJacobian() in the sub-models.

//! The frequency at which the states of the sub-models are synchronized depends on the solving scheme. 

//! Forward declarations
class Epetra_Comm;
class Epetra_IntVector;
class Epetra_Vector;

class CoupledModel
{
public:
	using VectorPtr = std::shared_ptr<Combined_MultiVec>;
	using LandMask  = Ocean::LandMask;
	
private:
	//! Adjacency graph describing the couplings between the models
	//! NOT IMPLEMENTED YET
    //! Graph couplings_;

	//! Trilinos-based parallel implicit ocean model (Trilinos-THCM) 
	std::shared_ptr<Ocean>  ocean_;

	//! Parallel implicit atmosphere model
	std::shared_ptr<AtmospherePar> atmos_;

	//! Combined state vector
	std::shared_ptr<Combined_MultiVec> stateView_;

	//! Combined solution vector
	std::shared_ptr<Combined_MultiVec> solView_;

	//! Combined rhs vectir
	std::shared_ptr<Combined_MultiVec> rhsView_;

	//! Determines the continuation parameter
	//! For a list see THCM::par2int()
	std::string parName_;
	
	//! Solving scheme
	//!   'D': decoupled
	//!   'B': backward block Gauss Seidel with sub-model solves
	//!   'F': FGMRES (Belos)
	char solvingScheme_;

	//! Max iterations in GS preconditioning scheme
	int iterGS_;

	//! Tolerance in GS preconditioning scheme
	double toleranceGS_;

	//! keep track of syncs
	int syncCtr_;

	//! Trilinos MPI-like communicator
	Teuchos::RCP<Epetra_Comm> comm_;
	
public:
	CoupledModel(std::shared_ptr<Ocean> ocean,
				 std::shared_ptr<AtmospherePar> atmos,
		         Teuchos::RCP<Teuchos::ParameterList> params);
	
	~CoupledModel() { INFO("CoupledModel destructor"); }

	//! Compute Jacobian matrix
	void computeJacobian();

	//! Compute RHS
	void computeRHS();

	//! Solve Jx=b
	void solve(std::shared_ptr<Combined_MultiVec> rhs);

	//! Initialize FGMRES (Belos) solver
	void initializeFGMRES();

	//! Apply the Jacobian matrix: out=J*v
	//!   coupled   mode 'C'
	//!   decoupled mode 'D'
	//! 
	void applyMatrix(Combined_MultiVec const &v, Combined_MultiVec &out, char mode = 'C');

	//! Apply the preconditioning: out=inv(P)*v
	//!   coupled   mode 'C'
	//!   decoupled mode 'D'
	void applyPrecon(Combined_MultiVec const &v, Combined_MultiVec &out, char mode = 'C');

	//! Apply the mass matrix: out=B*v
	void applyMassMatrix(Combined_MultiVec const &v, Combined_MultiVec &out);

	//! Get an RCP to the solution vector
	std::shared_ptr<Combined_MultiVec> getSolution(char mode = 'C');

	//! Get an RCP to the state vector
	std::shared_ptr<Combined_MultiVec> getState(char mode = 'C');

	//! Get an RCP to the rhs vector
	std::shared_ptr<Combined_MultiVec> getRHS(char mode = 'C');

	//! Set the continuation parameter
	void setPar(double value);
	
	//! Get the continuation parameter
	double getPar();

	//! Pre-processing.
	//! A controlling scheme might allow some pre-processing,
	//! for instance at the start of a Newton process.
	void preProcess();

	//! Post-processing.
	//! Similarly we can supply some post-processing, for instance when
	//! a Newton process has converged.
	void postProcess();

	//! Additional monitor for continuation
	bool monitor() { return false; }

private:
	//! Solve the system with a few symmetric block GS iterations
	void blockGSSolve(std::shared_ptr<Combined_MultiVec> rhs);
	
	//! Solve the system using FGMRES
	void FGMRESSolve(std::shared_ptr<Combined_MultiVec> rhs);
	
	//! For internal use: both modes for the underlying
	//! getSolution methods can be called here.
	std::shared_ptr<Combined_MultiVec> getSolution(char mode1, char mode2);

	//! Compute the residual ||b-Jx||
	double computeResidual(std::shared_ptr<Combined_MultiVec> rhs);

	//! coupled   mode 'C'
	//! decoupled mode 'D'
	double computeResidual(Combined_MultiVec const &rhs,
						   Combined_MultiVec const &x,
						   char mode = 'C');
		
	//! Synchronize the states between the models that are needed to communicate
	void synchronize();

};

#endif
