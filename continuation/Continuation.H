//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
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
	ds_(2.0e-2),               // initial continuation step-size
	dsMin_(1.0e-5),            // minimum continuation step-size
	dsMax_(1.0e01),            // maximum continuation step-size
	dsScale1_(1.20),           // conservative scaling for adaptive step
	dsScale2_(2.00),           // drastic scaling for adaptive step
	epsilon_(1.0e-8),          // variation used for numerical finite difference
	backTracking_(true),       // perform backtracking when Newton fails
	numBackTrackingSteps_(10), // initialize number of backtracking steps
	iter_(0),                  // initialize Newton iteration counter
	newtonIterations_(10),     // set the number of Newton iterations
	minAllowedIterations_(4),  // if less then we inrease the step size
	maxAllowedIterations_(5),  // if more then we decrease the step size
	newtonTolerance_(1.0e-3)   // set the tolerance in the Newton corrector
{
	model_  = model;                  // obtain the model
	// --> Right now it is somewhat arbitrary which VectorPtr contains a Copy
	//     and which is a View. This should be made more explicit. Perhaps
	//     in the naming of the members...
	state_  = model_->GetState('V');  // obtain a view of the state
	par_    = model_->GetPar();       // get the current parameter value
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Run()
{
	INFO("Continuation run... ");
	// set all the appropriate initial flags and counters
	step_ = 0;
	initialTangent_ = true;
	destinationReached_ = false;

	while (!destinationReached_)
	{
		// Perform a continuation step
		Step();

		// Get an intermediate solution
		model_->DumpState();

		// Detect special points:
		//  if necessary converge on them (-->not implemented yet)
		Detect();

		++step_;
	}
	INFO("Continuation run... done");
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Step()
{
	INFO("-----------------------------------------");
	INFO("Continuation status: ");
	INFO("            ds:  " << ds_);
	INFO("           par:  " << par_);
	INFO("          step:  " << step_);
	INFO("-----------------------------------------");
	
	// Create tangents.
	if (initialTangent_)
	{
		// make sure this is done only during a startup
		initialTangent_ = false;
		INFO("Create initial tangent...");
		CreateInitialTangent();
		INFO("Create initial tangent... done");
	}
	else
		CreateTangent();


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
	int status = NewtonCorrector();
	INFO("Apply Newton corrector... done");

	// If necessary reset the model, otherwise perform a normal
	// step adjustment based on the number of Newton iterations.
	// A reset will also set initialTangent = true, since we
	// lose the information to compute stateDot with a
	// finite difference.
	if (status == 1)
		Reset();
	else
		AdjustStep();
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

	// 1) Compute dFdPar_, force compute of RHS
	ComputeDFDPar('F');

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

	// 2) Compute stateDot = (state1 - state0)/ds
	// --> Note that our state_ member contains the
	//     current state of the model
	stateDot_->Update(1.0 / ds_, *state_, -1.0 / ds_);

	// 3) Compute parDot = (par1 - par0)/ds
	double storedPar  = model_->GetStoredPar();
	double currentPar = model_->GetPar();
	parDot_ = (currentPar - storedPar) / ds_;
}


//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::ComputeDFDPar(char mode)
{
	// Let the model compute an RHS F(par) with its current state
	// --> Note that this is usually a predicted state
	// --> The RHS for this particular state may already have been
	//     computed.
	//
	if (mode == 'F')
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
void Continuation<ModelPtr, VectorPtr>::EulerPredictor()
{
	// At the end of this function the model will be
	// in 'predicted' state.

	// Apply predictor to the state in the model
	// Compute: state = state0 + ds * statedot
	// --> 1) Note that at this point state0 and state are equal.
	// --> 2) Note that our state_ member contains a view of the
	//        current state of the model.
	state_->Update(ds_, *stateDot_, 1.0);

	// Compute  par = par0 + ds * pardot
	// and put it in the model. Make sure that our par_ and the value
	// in the model are the same.
	par_ = model_->GetStoredPar() + ds_ * parDot_;
	model_->SetPar(par_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
int Continuation<ModelPtr, VectorPtr>::NewtonCorrector()
{
	int backTrackingStatus;  // status of backtracking
	char mode;	             // mode for ComputeDFDPar()
	for (iter_ = 0; iter_ != newtonIterations_; ++iter_)
	{
		// Taking the derivative of the RHS w.r.t. the continuation
		// parameter using a finite difference. In the first iteration
		// the computation of the RHS is required.
		mode = (iter_ == 0) ? 'F' : 'A'; 
		ComputeDFDPar(mode);

		// Obtain the upper part (R) of the continuation RHS.
		// A copy of F(par) is obtained in ComputeDFDPar(), so
		// we use that one to get -F(par).
		// --> Note that we cannot use the current rhs in the model
		//     since at this point it corresponds to F(par+eps)
		VectorPtr R = rhs_;
		R->Scale(-1.0);

		// Obtain the lower part (r) of the continuation RHS,
		// Here we compute
		// r = ds - stateDot'*(state - storedState)
		//           - parDot*(par - storedPar)
		VectorPtr storedState = model_->GetStoredState('V');
		double      storedPar = model_->GetStoredPar();
		double r = ds_ - stateDot_->Dot(*state_)
			           - stateDot_->Dot(*storedState)
			           - parDot_ * (par_ - storedPar);

		// At this point the model contains the predicted state and
		// parameter. The Jacobian will be computed based on the
		// predicted data.
		model_->ComputeJacobian();

		// Now we will perform 2 solves to solve the bordered system:
		//  - The vector containing the first solution will have to make
		//    way for the second solution so we get a copy.
		//  - For the vector containing the second solution we can use
		//    a view.
		model_->Solve(dFdPar_);
		VectorPtr y1 = model_->GetSolution('C');
		model_->Solve(R);
		VectorPtr z1 = model_->GetSolution('V');

		// Determine the directions
		// First for the parameter:
		parDir_ = (r - stateDot_->Dot(*z1)) / (parDot_ - stateDot_->Dot(*y1));

		// Then for the state:
		//  - We perform an update on y1 and
		//  - keep the new direction there
		y1->Update(1.0, *z1, -1.0 * parDir_);
		stateDir_ = y1;

		// Update the state and the parameter in the model
		state_->Update(1.0, *stateDir_, 1.0);
		par_ = par_ + parDir_; // update our parameter
		model_->SetPar(par_);  // set it in the model

		// test new RHS, if it satisfies the tolerance we exit the loop
		normRHS_ = model_->GetNormRHS();
		model_->ComputeRHS();
		normRHStest_ = model_->GetNormRHS();
		INFO("   NewtonCorrector    iter: " << iter_ );
		INFO("   NewtonCorrector    norm: " << normRHStest_ );
		if (normRHStest_ < newtonTolerance_)
			break;

		// if necessary perform backtracking
		if (backTracking_ and (normRHS_ < normRHStest_) )
			backTrackingStatus = RunBackTracking();
	}
	if (iter_ == newtonIterations_)
	{
		WARNING("Continuation: Newton failed!!", __FILE__, __LINE__);
		return 1;
	}
	return 0;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
int Continuation<ModelPtr, VectorPtr>::RunBackTracking()
{
	// Initialize reduction with -1/2
	double reduction = -1.0 / 2;

	for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_)
	{
		if (normRHStest_ < normRHS_)
			break;

		// Apply reduction to the state vector and the parameter
		state_->Update(reduction, *stateDir_, 1.0); // update the state
		par_ = par_ - reduction * parDir_;  // update our parameter
		model_->SetPar(par_);		        // set it in the model

		// Compute new RHS and obtain its norm
		model_->ComputeRHS();
		normRHStest_ = model_->GetNormRHS();

		INFO("    NewtonCorrector:  backtracking: " );
		INFO("       step:      "   << backTrack_   );
		INFO("       reduction: "   << reduction    );
		INFO("       norm:      "   << normRHStest_ );

		// Update reduction
		reduction /= 2.0;
	}
	if (backTrack_ == numBackTrackingSteps_)
	{
		WARNING("Continuation: --> BACKTRACKING FAILED", __FILE__, __LINE__);
		return 1;
	}
	return 0;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Detect(char mode)
{
	// Just to be on the safe side
	par_ = model_->GetPar();
	double f1;
	f1 = par_ - model_->GetParDestination();

	// This is a naive initial implementation
	if (f1 > 0)
		destinationReached_ = true;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::AdjustStep()
{
	if (iter_ < minAllowedIterations_)
		ds_ = MIN(ds_ * dsScale1_, dsMax_);
	else if (iter_ > maxAllowedIterations_)
		ds_ = MAX(ds_ / dsScale1_, dsMin_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Reset()
{
	INFO("Continuation: RESET!");
	ds_ = MAX(ds_ / dsScale2_, dsMin_);
	INFO("Continuation:  Restoring model...");
	model_->Restore();
	INFO("Continuation:  Restoring model... done");
	initialTangent_ = true;
}
#endif
