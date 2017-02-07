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


public:
	//! Constructor based on paramterlist
	AtmospherePar(ParameterList params);
	
	~AtmospherePar(){INFO("AtmospherePar destructor");}

	
	
}
	
#endif
