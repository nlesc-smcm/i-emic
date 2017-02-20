#ifndef ATMOSPHEREPAR_H
#define ATMOSPHEREPAR_H

//! This class takes care of the parallelization of the atmosphere model.
//! For the domain decomposition we make use of the TRIOS_Domain stuff.
//! This als ensures that, in a coupled case, we can move parallel objects
//! around with maps that make sense.

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Atmosphere.H"
#include "TRIOS_Domain.H"
#include "GlobalDefinitions.H"


class AtmospherePar
{
public:
	using ParameterList = Teuchos::RCP<Teuchos::ParameterList>;
	
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

	//! parallel atmosphere state (non-overlapping)
	Teuchos::RCP<Epetra_Vector> rhs_;

	//! parallel atmosphere state (overlapping)
	Teuchos::RCP<Epetra_Vector> localRHS_;

	//! parallel atmosphere state (non-overlapping)
	Teuchos::RCP<Epetra_Vector> sst_;

	//! parallel atmosphere state (overlapping)
	Teuchos::RCP<Epetra_Vector> localSST_;

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

public:
	//! Constructor based on parameterlist
	AtmospherePar(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);
	
	~AtmospherePar() {INFO("AtmospherePar destructor");}

	void computeRHS();

	void computeJacobian();

	void applyMatrix(Teuchos::RCP<Epetra_Vector> const &in,
					 Teuchos::RCP<Epetra_Vector> &out);

	void applyPrecon(Teuchos::RCP<Epetra_Vector> const &in,
					 Teuchos::RCP<Epetra_Vector> &out);

	void solve(Teuchos::RCP<Epetra_Vector> rhs);	

	void setState(Teuchos::RCP<Epetra_Vector> input) { state_ = input; }

	void setOceanTemperature(Teuchos::RCP<Epetra_Vector> in);

	// Let the atmos object use idealized values,
	// mainly for testing purposes
	void idealized();

	// Obtain problem size
	int  dim() { return dim_; }

	Teuchos::RCP<Epetra_Vector> getState(char mode = 'C') { return getVector(mode, state_); }
	Teuchos::RCP<Epetra_Vector> getRHS(char mode = 'C')   { return getVector(mode, rhs_); }
	Teuchos::RCP<Epetra_Vector> getVector(char mode, Teuchos::RCP<Epetra_Vector> vec);


	Teuchos::RCP<Epetra_BlockMap> getStandardSurfaceMap() {return standardSurfaceMap_;}
	Teuchos::RCP<Epetra_BlockMap> getAssemblySurfaceMap() {return assemblySurfaceMap_;}

private:
	// A few functions borrowed from THCM, --> could be factorized when bored
	void createMatrixGraph();
	void insert_graph_entry(int* indices, int& pos,
							int i, int j, int k, int xx,
							int N, int M, int L) const;
};
	
#endif
