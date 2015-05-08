//======================================================================
#ifndef CONTINUATION_H
#define CONTINUATION_H

//======================================================================
#include "ContinuationDecl.H"

//======================================================================
//Constructor
template<typename ModelPtr, typename VectorPtr>
Continuation<ModelPtr, VectorPtr>::Continuation(ModelPtr model)
	:
	ds_(1.0),              // continuation step-size
	epsilon_(1.0e-8),      // variation used for numerical finite difference
	initialTangent_(true), // ensure correct computation of initial tangent
	newtonIterations_(10)  // set the number of iterations in the corrector 
{	
	model_  = model;                  // obtain the model
	state_  = model_->GetState('V');  // obtain the state
	
	par_ = model_->GetPar(); // get the current parameter value 
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::CreateInitialTangent()
{
	//-----------------------------------------------------------------
	// Initial tangent: 
	// 1) Take the derivative of the RHS w.r.t. the continuation par
	//    (dFdPar) using a finite difference.
	// 2) Solve J*statedot = -dFdPar.
	//-----------------------------------------------------------------

	// 1) Compute dFdPar_
	ComputeDFDPar();
	
	// 2) Solve J*stateDot_ = -dFdPar_
	model_->ComputeJacobian();
	dFdPar_->Scale(-1.0);
	model_->Solve(dFdPar_);

	// To obtain the solution from this solve without the risk
	// of altering it we get a copy.
	stateDot_ = model_->GetSolution('C');

	// Set the initial parameter tangent to 1
	parDot_ = 1.0;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::CreateTangent()
{
	// ------------------------------------------------------------
	// Calculate stateDot_:
	// --> we require the existence of a stored state and parameter
	// 1) Initialize with a copy of the stored state in the model 
	// 2) Calculate finite difference for the tangent based on the
	//    stored state and the current state
	// 3) Calculate finite difference for the parameter tangent
	// ------------------------------------------------------------

	// 1) Get a copy
	stateDot_ = model_->GetStoredState('C');

	// Obtain view of the current state
	VectorPtr currentState = model_->GetState('V');

	// 2) Compute stateDot = (state1 - state0)/ds
	stateDot_->Update(1.0 / ds_, *currentState, -1.0 / ds_);

	// 3) Compute parDot = (par1 - par0)/ds
	double storedPar  = model_->GetStoredPar();
	double currentPar = model_->GetPar();
	parDot_ = (currentPar - storedPar) / ds_;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::EulerPredictor()
{
	// At the end of this function the model will be
	// in 'predicted' state.
	
	// Obtain view of the current state
	VectorPtr currentState = model_->GetState('V');

	// Apply predictor to the state in the model
	// Compute: state = state0 + ds * statedot
	// Note that at this point state0 and state are equal.
	currentState->Update(ds_, *stateDot_, 1.0);

	// Compute  par = par0 + ds * pardot
	// and put it in the model. Make sure that our par_ and the value
	// in the model are the same.
	par_ = model_->GetStoredPar() + ds_ * parDot_;
	model_->SetPar(par_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::ComputeDFDPar()
{
	// Let the model compute an RHS F(par) with its current state
	// --> Note that this is usually a predicted state
	model_->ComputeRHS();
	
	// Get a copy of this RHS, store it in our rhs_ member
	rhs_ = model_->GetRHS('C');

	// Calculate new RHS
	model_->SetPar(par_ + epsilon_);  // increment parameter --> par + eps
	model_->ComputeRHS();             // compute new RHS     --> F(par+eps)
	model_->SetPar(par_);             // restore parameter   --> par

	// Get a copy of the new RHS F(par+eps), put it in dFdPar_
	dFdPar_ = model_->GetRHS('C');

	// Update dFdPar_ such that it contains (F(par+eps) - F(par))/eps
	dFdPar_->Update(-1.0 / epsilon_, *rhs_, 1.0 / epsilon_);	
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::NewtonCorrector()
{

}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Step()
{
	// Create tangents.
	// --> initialTangent_ is set to false in NewtonCorrector()
	if (initialTangent_)
	{
		INFO("Create initial tangent...");
		CreateInitialTangent();
		INFO("Create initial tangent... done");
	}
	else
	{
		INFO("Create tangent...");
		CreateTangent();
		INFO("Create tangent... done");
	}
	
	// After the stored data in the model is used for creating the
	// tangents, we are allowed to perform a new Store(),
	// overwriting the state, rhs and par with their current values.
	model_->Store();

	// We scale the tangent of the state with its length and norm.
	// --> Not sure if length is needed, but smaller is safer...
	int length = stateDot_->Length();
	double nrm = stateDot_->Norm();
	stateDot_->Scale(1.0 / (length * nrm));
	
	INFO("Apply Euler predictor...");
	EulerPredictor();
	INFO("Apply Euler predictor... done");

	INFO("Apply Newton corrector...");
	NewtonCorrector();
	INFO("Apply Newton corrector... done");
}

#endif
