//=====================================================================
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

//=====================================================================
#include "OceanCont.H"
#include "Vector.H"
#include "THCM.H"
#include "THCMdefs.H"
#include "GlobalDefinitions.H"

//=====================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//=====================================================================
// Fortran: get access to parameter get and set functions
extern "C"
{
	_SUBROUTINE_(getparcs)(int*,double*);
	_SUBROUTINE_(setparcs)(int*,double*); 
}

//=====================================================================
// OceanCont
//=====================================================================
OceanCont::OceanCont(Teuchos::RCP<Epetra_Comm> Comm)
	:
	Ocean(Comm)
{
	// THCM continuation parameter and its bounds
	parIdent_ = 19;
	parValue_ = 0.0;
	parStart_ = 0.0;
	parEnd_   = 1.0;

	// Put the initial value into the model.
	SetPar(parStart_);
}

//====================================================================
double OceanCont::GetPar()
{
	double thcmPar;
	FNAME(getparcs)(&parIdent_, &thcmPar);
	if (thcmPar != parValue_)
	{
		INFO("OceanCont::GetPar(): Parameter synchronization");
		INFO("              thcm: " << thcmPar);
		INFO("         OceanCont: " << parValue_);
		INFO("     fixing this...");
		FNAME(setparcs)(&parIdent_, &parValue_);
		FNAME(getparcs)(&parIdent_, &thcmPar);
		INFO("              thcm: " << thcmPar);
	}					  
	return parValue_;
}

//====================================================================
void OceanCont::SetPar(double value)
{
	parValue_ = value;
	FNAME(setparcs)(&parIdent_, &value);
}
