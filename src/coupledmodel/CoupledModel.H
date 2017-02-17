#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

#include "Ocean.H"
#include "Atmosphere.H"
#include "Combined_MultiVec.H"
#include "IDRSolver.H"
#include "GMRESSolver.H"

#include "CouplingBlock.H"

#include <vector>
#include <memory>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

//! This class assembles multiple sub-models into a CoupledModel. Which models and couplings are used should be supplied in the accompanying xml file (coupledmodel_params.xml). For the state we designed a special container called a SuperVector, see it's class header for more documentation.

//! Most routines in CoupledModel simply call their corresponding routines in the sub-models. For instance, a call to computeJacobian() translates to a call to computeJacobian() in the sub-models.

//! The frequency at which the states of the sub-models are synchronized depends on the solving scheme. Moreover we use hashes to make sure nothing is computed and synchronized too often

//! Forward declarations
class Epetra_Comm;

class CoupledModel
{
public:
	using VectorPtr = std::shared_ptr<Combined_MultiVec>;
	
private:
	//! Adjacency graph describing the couplings between the models
	//! NOT IMPLEMENTED YET
    //! Graph couplings_;

	//! Trilinos-based parallel implicit ocean model (Trilinos-THCM) 
	Teuchos::RCP<Ocean>  ocean_;

	CouplingBlock C12_, C21_;
	
	//! STL std::vector-based implicit atmosphere model with a banded Lapack solver
	std::shared_ptr<Atmosphere> atmos_;

	//! Combined state, right hand side and solution
	std::shared_ptr<Combined_MultiVec> stateView_;
	std::shared_ptr<Combined_MultiVec> solView_;
	std::shared_ptr<Combined_MultiVec> rhsView_;

	//! Let the models load existing states during initialization
	//! saving is done from the postprocess operation
	//! The respective in- and output files need to be set in the
	//! separate model parameters
	bool useExistingState_;

	//! Determines the continuation parameter
	//! For a list see THCM::par2int()
	std::string parName_;
	
	//! Solving scheme
	//!   'D': decoupled
	//!   'B': backward block Gauss Seidel with sub-model solves
	//!   'F': FGMRES (Belos)
	char solvingScheme_;

	bool useScaling_;  //! scale the problem before a solve
	
	int    iterGS_;      //! Max iterations in GS scheme
	double toleranceGS_; //! Tolerance in GS scheme

	//! Use hashing to decide to act based on whether a state has changed
	//! Hashes should be remembered for each task they get used for.
	bool useHash_;  //! flag to turn hashing on/off
	std::size_t syncHash_; //! hash to base synchronize() decision on
	std::size_t rhsHash_;  //! hash to base computeRHS() decision on	
	std::size_t jacHash_;  //! hash to base computeJacobian() decision on

	//! Synchronization counter
	int syncCtr_;

	//! Preconditioning strategy based on changes in state
	int buildPrecEvery_;

	//! Trilinos' MPI-like communicator
	Teuchos::RCP<Epetra_Comm> comm_;

	//! Atmosphere -> Ocean block in Jacobian
	std::shared_ptr<std::vector<double> > B_;	
	std::shared_ptr<std::vector<int> > rowsB_;

	//! Ocean -> Atmosphere block in Jacobian
	std::shared_ptr<std::vector<double> > C_;	
	std::shared_ptr<std::vector<int> > rowsC_;
	
public:
	CoupledModel(Teuchos::RCP<Ocean> ocean,
				 std::shared_ptr<Atmosphere> atmos,
		         Teuchos::RCP<Teuchos::ParameterList> params);
	
	~CoupledModel(){INFO("CoupledModel destructor");}

	void computeJacobian();
	void computeRHS();

	//! Solve the models
	void solve(std::shared_ptr<Combined_MultiVec> rhs = std::shared_ptr<Combined_MultiVec>());
	
	void initializeFGMRES();
	
	//! coupled   mode 'C'
	//! decoupled mode 'D'
	void applyMatrix(Combined_MultiVec const &v, Combined_MultiVec &out, char mode = 'C');
	void applyPrecon(Combined_MultiVec const &v, Combined_MultiVec &out, char mode = 'C');
	void applyMassMatrix(Combined_MultiVec const &v, Combined_MultiVec &out);

	//! these are deprecated
	std::shared_ptr<Combined_MultiVec> applyMatrix(Combined_MultiVec const &v);
	std::shared_ptr<Combined_MultiVec> applyPrecon(Combined_MultiVec const &v);

	std::shared_ptr<Combined_MultiVec> getSolution(char mode = 'C');
	std::shared_ptr<Combined_MultiVec> getState(char mode = 'C');
	std::shared_ptr<Combined_MultiVec> getRHS(char mode = 'C');

	//! Getters and setters for the continuation parameter
	//!  Used by Continuation class
	void   setPar(double value);
	double getPar();

	//! A controlling scheme might allow some pre-processing,
	//! for instance at the start of a Newton process.
	void preProcess();
	
	//! Similarly we can supply some post-processing, for instance when
	//! a Newton process has converged.
	void postProcess();

	//! Print the Jacobian
	void printJacobian(std::string const filename) const;

	void printResidual(VectorPtr b){}

	//! Additional monitor for continuation
	bool monitor(){return false;}

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

	//! Create a hash of the extended state + parameter in an XOR / rotating manner
	std::size_t getHash();

	//! Output, same as in Atmosphere but templated
	template<typename T>
	void write(std::vector<T> &vector, const std::string &filename);

	//! print profile
	void printProfile(ProfileType profile);
	
	//! A few tests are collected here
	void test();
};

#endif
