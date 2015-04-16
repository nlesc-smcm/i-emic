#include <Epetra_Comm.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Ocean.H"
#include "THCM.H"
#include "THCMdefs.H"

#include "GlobalDefinitions.H"

using Teuchos::RCP;
using Teuchos::rcp;


Ocean::Ocean(RCP<Epetra_Comm> Comm)
{
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");

	INFO("Entering Ocean constructor...");   
	// ---------------------------------------------------------------
	// Setup THCM parameters:
	// ---------------------------------------------------------------
	RCP<Teuchos::ParameterList> globalParamList =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("../ocean/parameters/thcm_params.xml",
								globalParamList.ptr());
	Teuchos::ParameterList &thcmList =
		globalParamList->sublist("THCM");
	DEBUG(*globalParamList);
	DEBUG(thcmList);
    //-------------------------------------------------------------------
	// Create THCM object
	//-------------------------------------------------------------------
    thcm_ = rcp(new THCM(thcmList, Comm));
	INFO("Leaving Ocean constructor...");
}

Ocean::~Ocean()
{
	INFO("Ocean destructor called...");
}

void Ocean::doStuff()
{
	RCP<Epetra_Vector> soln = THCM::Instance().getSolution();
	INFO("STUFF GEDAAN")
}
