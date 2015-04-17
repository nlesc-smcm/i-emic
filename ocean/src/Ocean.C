#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>

#include "Ocean.H"
#include "THCM.H"
#include "THCMdefs.H"

#include "GlobalDefinitions.H"

using Teuchos::RCP;
using Teuchos::rcp;

Ocean::Ocean(RCP<Epetra_Comm> Comm)
	:
	solverInitialized_(false)
{
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");

	DEBUG("Entering Ocean constructor...");   
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
	//-------------------------------------------------------------------
	// Obtain solution vector from THCM
	//  THCM is implemented as a Singleton, which allows only a single
	//  instance at a time. The Ocean class can access THCM with a call
	//  to THCM::Instance()
	//-------------------------------------------------------------------
	state_ = THCM::Instance().getSolution();
	INFO("  Obtained solution from THCM");
	//-------------------------------------------------------------------
	// Randomize state vector 
	//-------------------------------------------------------------------
	double randScale = 1.0e-3;
	randomizeState(randScale);
	INFO("  Randomized solution and scaled with a factor " << randScale);
	//-------------------------------------------------------------------
	// Obtain Jacobian from THCM
    //-------------------------------------------------------------------
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	jac_ = THCM::Instance().getJacobian();
	INFO("  Obtained jacobian from THCM");
    //-------------------------------------------------------------------
	// Initialize a few datamembers
	//-------------------------------------------------------------------
	rhs_ = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	DEBUG("Leaving Ocean constructor...");	
}

void Ocean::randomizeState(double scaling)
{
	DEBUG("Entering Ocean::randomizeState()...");
	state_->Random();
	state_->Scale(scaling);
	INFO("Initialized solution vector");
	DEBUG("Leaving  Ocean::randomizeState()...");
}

void Ocean::initializeSolver()
{
	DEBUG("Entering Ocean::initializeSolver()...");
	problem_ =
		rcp(new Belos::LinearProblem<double, Epetra_MultiVector,
			Epetra_Operator>(jac_, dir_, rhs_) );
	
	solverInitialized_ = true;
	DEBUG("Leaving Ocean::initializeSolver()...");
}

void Ocean::solve()
{
	DEBUG("Entering Ocean::solve()");
	if (!solverInitialized_)
		initializeSolver();
	DEBUG("Leaving Ocean::solve()");
	
}

void Ocean::computeJacobian()
{
	DEBUG("Entering Ocean::computeJacobian()...");
	// Compute the Jacobian in THCM using the current state
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	// Get the Jacobian from THCM
	jac_ = THCM::Instance().getJacobian();	
	DEBUG("Leaving  Ocean::computeJacobian()...");
}

OceanTheta::OceanTheta(Teuchos::RCP<Epetra_Comm> Comm)
	:
	Ocean(Comm),
	theta_(0.0),
	timestep_(0.001)
{
	DEBUG("Entering OceanTheta constructor");
	//-------------------------------------------------------------------
	// Initialize a few datamembers
	//-------------------------------------------------------------------
	oldState_ = rcp(new Epetra_Vector(jac_->OperatorDomainMap()));
	stateDot_ = rcp(new Epetra_Vector(jac_->OperatorDomainMap()));
	oldRhs_   = rcp(new Epetra_Vector(jac_->OperatorRangeMap()));
	DEBUG("Leaving OceanTheta constructor");
}

void OceanTheta::parkModel()
{
	DEBUG("Entering Ocean::parkModel()");
	INFO("Parking the model"); 
	*oldState_ = *state_;
	*oldRhs_   = *rhs_;
	DEBUG("Leaving Ocean::parkModel()");
}

void OceanTheta::computeJacobian()
{
	DEBUG("Entering OceanTheta::computeJacobian()...");
	// Check theta
	if (theta_ < 0 || theta_ > 1)
	{
		INFO("Incorrect theta: " << theta_);
	}	
	// First compute the Jacobian in THCM using the current state
	THCM::Instance().evaluate(*state_, Teuchos::null, true);
	// Get the plain Jacobian from THCM
	jac_ = THCM::Instance().getJacobian();
	// Scale it with theta
	jac_->Scale(theta_);
	// Get the mass matrix from THCM (which is actually just a
	//    vector with diagonal elements)
	// Wrap it in a non-owning RCP
	massMatrix_ = rcp(&THCM::Instance().DiagB(), false);

	// Get the number of local elements
	int numMyElements     =	jac_->Map().NumMyElements();
	// Get a list of the global element IDs owned by the calling proc
	int *myGlobalElements = jac_->Map().MyGlobalElements();

    // Add to the Jacobian the values B[i]/dt
	double value;
	for (int i = 0; i != numMyElements; ++i)
	{
		value = (*massMatrix_)[i] / timestep_;
		jac_->SumIntoGlobalValues(myGlobalElements[i], 1,
								  &value, myGlobalElements + i);
	}
	jac_->FillComplete();
	DEBUG("Leaving  OceanTheta::computeJacobian()...");
}

void Ocean::computeRHS()
{
	DEBUG("Entering Ocean::computeRHS()");
	THCM::Instance().evaluate(*state_, rhs_, false);
	DEBUG("Leaving Ocean::computeRHS()");
}

void OceanTheta::computeRHS()
{
	DEBUG("Entering OceanTheta::computeRHS()");
	THCM::Instance().evaluate(*state_, rhs_, false);
	// Calculate mass matrix
	THCM::Instance().evaluateB();
	// Get the mass matrix from THCM 
	massMatrix_ = rcp(&THCM::Instance().DiagB(), false);
	// Calculate d/dt x = (xnew - xold)/dt
	stateDot_->Update(1.0 / timestep_, *state_, -1.0 / timestep_,
					  *oldState_, 0.0);
	// Obtain number of local elements
	int numMyElements     = stateDot_->Map().NumMyElements();
	// Get a list of the global element IDs owned by the calling proc
	int *myGlobalElements = stateDot_->Map().MyGlobalElements();
	// Scale xdot with the values in the mass matrix
	double value;
	for (int i = 0; i != numMyElements; ++i)
	{
		value = (*massMatrix_)[i] * (*stateDot_)[i];
		stateDot_->ReplaceGlobalValues(1, &value, myGlobalElements + i);
	}
	// The final theta timestepping rhs is given by
	// B d/dt x + theta*F(x) + (theta-1) * F(x_old)
	rhs_->Update(theta_ - 1.0, *oldRhs_, theta_);
	rhs_->Update(1.0, *stateDot_, 1.0);
		
	DEBUG("Leaving Ocean::computeRHSTheta()");
}


