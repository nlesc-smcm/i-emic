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

public:
	//! Constructor based on paramterlist
	AtmospherePar(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);
	
	~AtmospherePar() {INFO("AtmospherePar destructor");}

};
	
#endif
