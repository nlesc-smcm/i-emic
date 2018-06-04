//=====================================================================
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

//=====================================================================
#include "OceanTheta.H"
#include "THCM.H"
#include "THCMdefs.H"
#include "GlobalDefinitions.H"

//=====================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//=====================================================================
// OceanTheta
//=====================================================================
OceanTheta::OceanTheta(Teuchos::RCP<Epetra_Comm> Comm,
					   Teuchos::RCP<Teuchos::ParameterList> params)
	:
	Ocean(Comm, params),
	theta_(1.0),
	timestep_(1.0e-03)
{

	// Initialize a few datamembers
	oldState_ = rcp(new Epetra_Vector(*state_));
	stateDot_ = rcp(new Epetra_Vector(*state_));
	oldRhs_   = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
}

//=====================================================================
void OceanTheta::store()
{
	DEBUG("Storing the model");
	
	if (oldState_ == Teuchos::null or
		!oldState_->Map().SameAs(state_->Map()))
	{
		oldState_ =
			rcp(new Epetra_Vector(state_->Map()));
	}
	
	*oldState_ = *state_;

	if (oldRhs_ == Teuchos::null or
		!oldRhs_->Map().SameAs(rhs_->Map()))
	{
		oldRhs_ =
			rcp(new Epetra_Vector(rhs_->Map()));
	}
	
	*oldRhs_   = *rhs_;
}

//=====================================================================
void OceanTheta::restore()
{
	DEBUG("Restoring the model");

	if (state_ == Teuchos::null or
		!state_->Map().SameAs(oldState_->Map()))
	{
		state_ =
			rcp(new Epetra_Vector(oldState_->Map()));
	}
	
	*state_ = *oldState_;
	
	if (rhs_ == Teuchos::null or
		!rhs_->Map().SameAs(oldRhs_->Map()))
	{
		rhs_ =
			rcp(new Epetra_Vector(oldRhs_->Map()));
	}
	
	*rhs_   = *oldRhs_;
}

//=====================================================================
void OceanTheta::computeRHS()
{
    // THCM::Instance().fixMixing(0);
	THCM::Instance().evaluate(*state_, rhs_, false);

    // Calculate mass matrix
	THCM::Instance().evaluateB();

    // Get the mass matrix from THCM 
	diagB_ = THCM::Instance().DiagB();

    // Create timestep scaling vector
    // Epetra_Vector leftScale(*diagB_);

    // int numMyElements = leftScale.Map().NumMyElements();
    // for (int i = 0; i != numMyElements; ++i)
    //     if (std::abs((*diagB_)[i]) < 1e-12)
    //         leftScale[i] = 1.0;
    //     else
    //         leftScale[i] = timestep_;

    // Calculate dx = (xnew - xold)
    double s = 1.0 / timestep_ / theta_;
	stateDot_->Update(s , *state_, -s ,
					  *oldState_, 0.0);
	
	// Obtain number of local elements
	int numMyElements = stateDot_->Map().NumMyElements();

    // Get a list of the global element IDs owned by the calling proc
	int *myGlobalElements = stateDot_->Map().MyGlobalElements();

	// Scale xdot with the values in the mass matrix
	double value;
	for (int i = 0; i != numMyElements; ++i)
	{
		value = (*diagB_)[i] * (*stateDot_)[i];
		stateDot_->ReplaceGlobalValues(1, &value, myGlobalElements + i);
	}

    // The final theta timestepping rhs is given by
	// -1 * ( (B dx/dt)/theta + F(x) + (theta-1)/theta * F(x_old))
    double f = (theta_-1)/theta_;
    rhs_->Update(f, *oldRhs_, 1.0);
	rhs_->Update(1.0, *stateDot_, 1.0);
	rhs_->Scale(-1.0);
}

//=====================================================================
void OceanTheta::computeJacobian()
{
    // Check theta
	if (theta_ <= 0 || theta_ > 1)
	{
		WARNING("Ocean: Incorrect theta: " << theta_,
				__FILE__, __LINE__);
	}	

    // First compute the Jacobian in THCM using the current state
    // THCM::Instance().fixMixing(0);
	THCM::Instance().evaluate(*state_, Teuchos::null, true);

    // Get the plain Jacobian from THCM
	jac_ = THCM::Instance().getJacobian();

    // Get the mass matrix from THCM (which is actually just a
	//    vector with diagonal elements)
	// Wrap it in a non-owning RCP
	diagB_ = THCM::Instance().DiagB();

    // Epetra_Vector leftScale(*diagB_);
    // int numMyElements = leftScale.Map().NumMyElements();
    // for (int i = 0; i != numMyElements; ++i)
    //     if (std::abs((*diagB_)[i]) < 1e-12)
    //         leftScale[i] = 1.0;
    //     else
    //         leftScale[i] = timestep_ * theta_;

    // CHECK_ZERO( jac_->LeftScale(leftScale) );
      
    // Get the number of local elements
	int numMyElements =	jac_->Map().NumMyElements();

    // Get a list of the global element IDs owned by the calling proc
	int *myGlobalElements = jac_->Map().MyGlobalElements();

    // Add to the Jacobian the values B[i]. The result is B / (theta * dt) + J.
	double value;
	for (int i = 0; i != numMyElements; ++i)
	{
		value = (*diagB_)[i] / timestep_ / theta_;
		jac_->SumIntoGlobalValues(myGlobalElements[i], 1,
								  &value, myGlobalElements + i);
	}
	jac_->FillComplete();
}
