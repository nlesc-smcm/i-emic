#ifndef ATMOSPHEREPAR_H
#define ATMOSPHEREPAR_H

//! This class takes care of the parallelization of the atmosphere model.
//! For the domain decomposition we make use of the TRIOS_Domain stuff.
//! This als ensures that, in a coupled case, we can move parallel objects
//! around with maps that make sense.

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_CrsGraph.h>
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include "Atmosphere.H"
#include "TRIOS_Domain.H"
#include "GlobalDefinitions.H"

class Ocean;

class AtmospherePar
{
public:
	using ParameterList = Teuchos::RCP<Teuchos::ParameterList>;
	using VectorPtr = Teuchos::RCP<Epetra_Vector>;
	
private:	
	//! serial atmosphere class, working on local subdomains
	std::shared_ptr<Atmosphere> atmos_;

	//! parameterlist
	ParameterList params_;

	//! communicator object
	Teuchos::RCP<Epetra_Comm> comm_;

	//! grid size
	int n_,m_,l_, dim_;

	//! periodic bdc
	bool periodic_;

	//! degrees of freedom
	int dof_;
	
	//! grid dimensions
	double xmin_,xmax_,ymin_,ymax_;

	//! TRIOS domain object for parallelization strategy
	Teuchos::RCP<TRIOS::Domain> domain_;

	//! standard map, without ghost nodes (non-overlapping).
	Teuchos::RCP<Epetra_Map> standardMap_;

	//! assembly map, with ghost nodes (overlapping).
	Teuchos::RCP<Epetra_Map> assemblyMap_;

	//! depth-averaged, single unknown for ocean surface temperature
	//! standard map, without ghost nodes (non-overlapping).
	Teuchos::RCP<Epetra_Map> standardSurfaceMap_;

	//! depth-averaged, single unknown for ocean surface temperature
	//! assembly map, with ghost nodes (overlapping).
	Teuchos::RCP<Epetra_Map> assemblySurfaceMap_;

	//! parallel atmosphere state (non-overlapping)
	Teuchos::RCP<Epetra_Vector> state_;

	//! parallel atmosphere state (overlapping)
	Teuchos::RCP<Epetra_Vector> localState_;

 	//! parallel atmosphere rhs (non-overlapping)
	Teuchos::RCP<Epetra_Vector> rhs_;

	//! parallel atmosphere rhs (overlapping)
	Teuchos::RCP<Epetra_Vector> localRHS_;

	//! parallel atmosphere sst (non-overlapping)
	Teuchos::RCP<Epetra_Vector> sst_;

	//! parallel atmosphere sst (overlapping)
	Teuchos::RCP<Epetra_Vector> localSST_;
	
	//! parallel atmosphere sol (non-overlapping)
	Teuchos::RCP<Epetra_Vector> sol_;

	//! parallel atmosphere sol (overlapping)
	Teuchos::RCP<Epetra_Vector> localSol_;

	//! input file (HDF5)
	std::string const inputFile_;

	//! output file (HDF5)
	std::string const outputFile_;

	//! input/output flags
	bool loadState_, saveState_;

	//! Jacobian matrix
	Teuchos::RCP<Epetra_CrsMatrix> jac_;

	//! matrix graph, necessary to build Jacobian matrix
	Teuchos::RCP<Epetra_CrsGraph> matrixGraph_;

	//! global surface mask
	std::shared_ptr<std::vector<int> > surfmask_;

public:
	//! Constructor based on parameterlist
	AtmospherePar(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);
	
	~AtmospherePar() {INFO("AtmospherePar destructor");}

	void computeRHS();

	void computeJacobian();

	void applyMatrix(Epetra_MultiVector const &in,
					 Epetra_MultiVector &out);

	void applyPrecon(Epetra_MultiVector &in,
					 Epetra_MultiVector &out);

	void solve(Teuchos::RCP<Epetra_MultiVector> const &b);

	void setState(Teuchos::RCP<Epetra_Vector> input) { state_ = input; }

	void setOceanTemperature(Teuchos::RCP<Epetra_Vector> in);

	void setLandMask(Utils::MaskStruct const &mask);

	//! set continuation parameter
	void setPar(double value) {atmos_->setPar(value); }

	//! set parameter
	void setPar(std::string const &parName, double value) {atmos_->setPar(parName, value); }

	//! Set continuation parameter name
	void setParName(std::string const &parName) { atmos_->setParName(parName); }

	//! Get continuation parameter
	double getPar() { return atmos_->getPar(); }

	//! Get continuation parameter
	double getPar(std::string const &parName) { return atmos_->getPar(parName); }

	void preProcess() {}
	void postProcess() {}

	//! Let the atmos object use idealized values,
	//! mainly for testing purposes
	void idealized();

	//! Obtain problem size
	int dim() { return dim_; }

	//! Obtain degrees of freedom
	int dof() { return dof_; }
	
	Teuchos::RCP<Epetra_Vector> getState(char mode = 'C') { return getVector(mode, state_); }
	Teuchos::RCP<Epetra_Vector> getSolution(char mode = 'C') { return getVector(mode, sol_); } 
	Teuchos::RCP<Epetra_Vector> getRHS(char mode = 'C')   { return getVector(mode, rhs_); }
	Teuchos::RCP<Epetra_Vector> getVector(char mode, Teuchos::RCP<Epetra_Vector> vec);
	Teuchos::RCP<Epetra_CrsMatrix> getJacobian() { return jac_; }

	//! Return pointer to domain object
	Teuchos::RCP<TRIOS::Domain> getDomain() { return domain_; }

	Teuchos::RCP<Epetra_Vector> getT();

	Teuchos::RCP<Epetra_BlockMap> getStandardSurfaceMap() {return standardSurfaceMap_;}
	Teuchos::RCP<Epetra_BlockMap> getAssemblySurfaceMap() {return assemblySurfaceMap_;}

	//! Get ocean dependency
	std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Ocean> ocean);

	//! Get row at the interface, XX is 1-based, i and j 0-based
	int interface_row(int i, int j, int XX)
		{ return FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_, i, j, l_-1, ATMOS_TT_); }

private:
	//! solve with right hand side b on assembly subdomains
	//! using dgbtrf and dgbtrs
	void solveSubDomain(Teuchos::RCP<Epetra_MultiVector> const &b);

	//! create dependency graph for the Jacobian matrix
	//! function borrowed from THCM,
	//! --> could be factorized within the project when bored
	void createMatrixGraph();

	//! insert a dependency into the graph
	//! function borrowed from THCM,
	//! --> could be factorized within the project when bored
	void insert_graph_entry(int* indices, int& pos,
							int i, int j, int k, int xx,
							int N, int M, int L) const;
};
	
#endif
