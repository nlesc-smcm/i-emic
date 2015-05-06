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
	// Initialize a few datamembers
	oldState_ = rcp(new Epetra_Vector(*state_));
	oldRhs_   = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));

	// THCM continuation parameter and its bounds
	thcmParIdent_ = 19;
	thcmPar_      = 0.0;
	thcmParStart_ = 0.0;
	thcmParEnd_   = 1.0;

	// Put the initial value into the model.
	SetPar(thcmParStart_);
}

//=====================================================================
void OceanCont::StoreState()
{
	DEBUG("Storing the state of OceanCont");
	
	if (oldState_ == Teuchos::null or
		!oldState_->Map().SameAs(state_->Map()))
	{
		oldState_ =
			rcp(new Epetra_Vector(state_->Map()));
	}
	
	*oldState_ = *state_;
}

//=====================================================================
void OceanCont::StoreRHS()
{
	DEBUG("Storing the RHS of OceanCont");
	if (oldRhs_ == Teuchos::null or
		!oldRhs_->Map().SameAs(rhs_->Map()))
	{
		oldRhs_ =
			rcp(new Epetra_Vector(rhs_->Map()));
	}
	
	*oldRhs_   = *rhs_;
}

//=====================================================================
void OceanCont::RestoreState()
{
	DEBUG("Restoring the state of OceanCont");

	if (state_ == Teuchos::null or
		!state_->Map().SameAs(oldState_->Map()))
	{
		state_ =
			rcp(new Epetra_Vector(oldState_->Map()));
	}
	
	*state_ = *oldState_;
}

//====================================================================
void OceanCont::RestoreRHS()
{
	DEBUG("Restoring the RHS of OceanCont");

	if (rhs_ == Teuchos::null or
		!rhs_->Map().SameAs(oldRhs_->Map()))
	{
		rhs_ =
			rcp(new Epetra_Vector(oldRhs_->Map()));
	}
	
	*rhs_   = *oldRhs_;
}

//====================================================================
Teuchos::RCP<Epetra_Vector> OceanCont::GetSolutionCopy()
{
	Teuchos::RCP<Epetra_Vector> copySol =
		Teuchos::rcp(new Epetra_Vector(*sol_));
	return copySol;
}

//====================================================================
double OceanCont::GetPar()
{
	FNAME(getparcs)(&thcmParIdent_, &thcmPar_);
	return thcmPar_;
}

//====================================================================
void OceanCont::SetPar(double value)
{
	FNAME(setparcs)(&thcmParIdent_, &value);
}
