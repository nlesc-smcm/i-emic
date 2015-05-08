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
	initialTangent_(true)  // ensure correct computation of initial tangent
{	
	model_  = model;                  // obtain the model
	state_  = model_->GetState('V');  // obtain the state
	
	// get the current parameter value from the model
	par_ = model_->GetPar();
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::CreateInitialTangent()
{
	//-----------------------------------------------------------------
	// Initial tangent: 
	// Taking the derivative of the RHS w.r.t. the continuation param
	// using a finite difference.
	// Then we solve J*soldot = -dFdPar.
	//-----------------------------------------------------------------
	model_->ComputeRHS();   // compute RHS --> F(par)
	model_->Store();        // let the model store the computed RHS

	// Obtain a view of F(par) as we are not going to change it.
	VectorPtr oldRhs = model_->GetStoredRHS('V'); 

	model_->SetPar(par_ + epsilon_); // increment parameter --> par + eps
	model_->ComputeRHS();            // compute new RHS     --> F(par+eps)
	model_->SetPar(par_);            // restore parameter   --> par

	// Initialize dFdPar_ with a copy of the newly computed RHS: F(par+eps)
	dFdPar_ = model_->GetRHS('C');

	// Update dFdPar_ such that it contains (F(par+eps) - F(par))/eps
	dFdPar_->Update(-1.0 / epsilon_, *oldRhs, 1.0 / epsilon_);
	
	// Solve J*solDot_ = -dFdPar_
	model_->ComputeJacobian();
	dFdPar_->Scale(-1.0);
	model_->Solve(dFdPar_);

	// To obtain the solution from this solve without the risk
	// of altering it we get a copy.
	stateDot_ = model_->GetSolution('C');

	// The inital parameter tangent
	parDot_ = 1.0;

	// Restore the RHS in the model to F(par)
	// --> Note that this does not affect affect dFdPar_
	model_->Restore(); 
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::CreateTangent()
{	
	// Finite difference for the tangent based on the old
	// state and the new state...
	// Calculate stateDot_, initialize with copy of the
	// old state in the model as we are going to overwrite
	// the old state later on.
	stateDot_ = model_->GetStoredState('C');

	// Obtain view of the current state
	VectorPtr currentState = model_->GetState('V');

	// Compute stateDot = (state1 - state0)/ds
	stateDot_->Update(1.0 / ds_, *currentState, -1.0 / ds_);

	// Compute parDot = (par1 - par0)/ds
	double oldPar = model_->GetStoredPar();
	double curPar = model_->GetPar();
	parDot_ = (curPar - oldPar) / ds_;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Step()
{
	// Create tangents.
	if (initialTangent_)
	{
		INFO("Create initial tangent...");
		CreateInitialTangent();
		INFO("Create initial tangent... done");
		initialTangent_ = false; // make sure this is done once
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
	
	// ------------------- Euler predictor ------------------------
	INFO("Compute Euler predictor...");
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
	INFO("Compute Euler predictor... done");
	// -----------------------------------------------------------

	// ------------------- Newton corrector-----------------------
	
	// -----------------------------------------------------------
	
}

#endif
