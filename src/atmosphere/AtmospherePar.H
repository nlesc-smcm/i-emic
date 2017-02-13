#ifndef ATMOSPHEREPAR_H
#define ATMOSPHEREPAR_H

//! For the parallelization strategy we make use of the TRIOS routines

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
	//! serial atmosphere class, working on subdomains
	std::shared_ptr<Atmosphere> atmos_;

	//! parameterlist
	ParameterList params_;

	//! communicator object
	Teuchos::RCP<Epetra_Comm> comm_;

	//! grid size
	int n_,m_,l_;

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

	//! parallel atmosphere state (non-overlapping)
	Teuchos::RCP<Epetra_Vector> state_;

	//! parallel atmosphere state (overlapping)
	Teuchos::RCP<Epetra_Vector> localState_;

	//! parallel atmosphere state (non-overlapping)
	Teuchos::RCP<Epetra_Vector> rhs_;

	//! parallel atmosphere state (overlapping)
	Teuchos::RCP<Epetra_Vector> localRHS_;

	//! input file (HDF5)
	std::string const inputFile_;

	//! output file (HDF5)
	std::string const outputFile_;

	//! input/output flags
	bool loadState_, saveState_;

public:
	//! Constructor based on parameterlist
	AtmospherePar(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);
	
	~AtmospherePar() {INFO("AtmospherePar destructor");}

	void computeRHS();

	void setState(Teuchos::RCP<Epetra_Vector> input);

	// Let the atmos object use idealized values,
	// mainly for testing purposes
	void idealized() { atmos_->idealized(); }
		
};
	
#endif
