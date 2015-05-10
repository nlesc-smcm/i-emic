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
	// Initialize a few datamembers
	storedState_ = rcp(new Epetra_Vector(*state_));
	storedRhs_   = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));

	// THCM continuation parameter and its bounds
	parIdent_ = 19;
	parValue_ = 0.0;
	parStart_ = 0.0;
	parEnd_   = 1.0;

	// Put the initial value into the model.
	SetPar(parStart_);
}

//=====================================================================
void OceanCont::Store()
{
	StoreState();
	StoreRHS();
	StorePar();
}

//=====================================================================
void OceanCont::Restore()
{
	RestoreState();
	RestoreRHS();
	RestorePar();
}

//=====================================================================
void OceanCont::StoreState()
{
	DEBUG("Storing the state of OceanCont");
	
	if (storedState_ == Teuchos::null or
		!storedState_->Map().SameAs(state_->Map()))
	{
		storedState_ =
			rcp(new Epetra_Vector(state_->Map()));
	}	
	*storedState_ = *state_;
}

//=====================================================================
void OceanCont::StoreRHS()
{
	DEBUG("Storing the RHS of OceanCont");
	if (storedRhs_ == Teuchos::null or
		!storedRhs_->Map().SameAs(rhs_->Map()))
	{
		storedRhs_ =
			rcp(new Epetra_Vector(rhs_->Map()));
	}
	
	*storedRhs_   = *rhs_;
}

//=====================================================================
void OceanCont::StorePar()
{
	storedPar_ = parValue_  ;
}

//=====================================================================
void OceanCont::RestoreState()
{
	DEBUG("Restoring the state of OceanCont");

	if (state_ == Teuchos::null or
		!state_->Map().SameAs(storedState_->Map()))
	{
		state_ =
			rcp(new Epetra_Vector(storedState_->Map()));
	}
	
	*state_ = *storedState_;
}

//====================================================================
void OceanCont::RestoreRHS()
{
	DEBUG("Restoring the RHS of OceanCont");

	if (rhs_ == Teuchos::null or
		!rhs_->Map().SameAs(storedRhs_->Map()))
	{
		rhs_ =
			rcp(new Epetra_Vector(storedRhs_->Map()));
	}
	
	*rhs_   = *storedRhs_;
}

//=====================================================================
void OceanCont::RestorePar()
{
	parValue_  = storedPar_ ;
	
	// Let THCM know that we want to restore the parameter
	SetPar(parValue_);
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

//====================================================================
Teuchos::RCP<Vector> OceanCont::GetStoredState(char mode)
{
	if (mode == 'C')
	{
		RCP<Epetra_Vector> copyStoredState =
			rcp(new Epetra_Vector(*storedState_));
		RCP<Vector> storedStatePtr =
			rcp(new Vector(copyStoredState));
		return storedStatePtr;
	}
	else if (mode == 'V')
	{
		RCP<Vector> storedStatePtr =
			rcp(new Vector(storedState_));
		return storedStatePtr;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//====================================================================
Teuchos::RCP<Vector> OceanCont::GetStoredRHS(char mode)
{
	if (mode == 'C')
	{
		RCP<Epetra_Vector> copyStoredRhs =
			rcp(new Epetra_Vector(*storedRhs_));
		RCP<Vector> storedRhsPtr =
			rcp(new Vector(copyStoredRhs));
		return storedRhsPtr;
	}
	else if (mode == 'V')
	{
		RCP<Vector> storedRhsPtr =
			rcp(new Vector(storedRhs_));
		return storedRhsPtr;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	} 
}
