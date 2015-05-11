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
	storedState1_ = rcp(new Epetra_Vector(*state_));
	storedState2_ = rcp(new Epetra_Vector(*state_));
	storedRhs_    = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));

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
	
	// First we put storedState1 -> storedState2
	if (storedState2_ == Teuchos::null or
		!storedState2_->Map().SameAs(storedState1_->Map()))
	{
		storedState2_ =
			rcp(new Epetra_Vector(storedState1_->Map()));
	}	
	*storedState2_ = *storedState1_;
	
	// Then we put state -> storedState1
	if (storedState1_ == Teuchos::null or
		!storedState1_->Map().SameAs(state_->Map()))
	{
		storedState1_ =
			rcp(new Epetra_Vector(state_->Map()));
	}	
	*storedState1_ = *state_;
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
	storedPar2_ = storedPar1_;
	storedPar1_ = parValue_  ;
}

//=====================================================================
void OceanCont::RestoreState()
{
	DEBUG("Restoring the state of OceanCont");

	// First restore state
	if (state_ == Teuchos::null or
		!state_->Map().SameAs(storedState1_->Map()))
	{
		state_ =
			rcp(new Epetra_Vector(storedState1_->Map()));
	}	
	*state_ = *storedState1_;

	// Then restore storedState1
	if (storedState1_ == Teuchos::null or
		!storedState1_->Map().SameAs(storedState2_->Map()))
	{
		storedState1_ =
			rcp(new Epetra_Vector(storedState2_->Map()));
	}	
	*storedState1_ = *storedState2_;
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
	parValue_  = storedPar1_ ;
	// Let THCM know that we restore the parameter
	SetPar(parValue_);
	storedPar1_ = storedPar2_;
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
			rcp(new Epetra_Vector(*storedState1_));
		RCP<Vector> storedStatePtr =
			rcp(new Vector(copyStoredState));
		return storedStatePtr;
	}
	else if (mode == 'V')
	{
		RCP<Vector> storedStatePtr =
			rcp(new Vector(storedState1_));
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
