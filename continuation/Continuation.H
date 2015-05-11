//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef CONTINUATION_H
#define CONTINUATION_H

//======================================================================
#include "ContinuationDecl.H"
#include "GlobalDefinitions.H"

//======================================================================
//Constructor
template<typename ModelPtr, typename VectorPtr>
Continuation<ModelPtr, VectorPtr>::Continuation(ModelPtr model)
	:
	ds_(1.0e-2),               // initial continuation step-size
	dsMin_(1.0e-5),            // minimum continuation step-size
	dsMax_(1.0e01),            // maximum continuation step-size
	scale1_(1.00),             // conservative scaling for adaptive step
	scale2_(2.00),             // drastic scaling for resetting adaptive step
	epsilon_(1.0e-8),          // variation used for numerical finite difference
	backTracking_(true),       // perform backtracking when Newton fails
	forceIncreaseParam_(true), // only allow increasing parameter
	numBackTrackingSteps_(10), // number of backtracking steps
	newtonIterations_(15),     // set the max for the Newton iterations
	minAllowedIterations_(4),  // if less then we increase the step size
	maxAllowedIterations_(5),  // if more then we decrease the step size
	newtonTolerance_(1.0e-4),  // set the tolerance in the Newton corrector
	parDotSign_(1)             // sign of parDot, used in Detect()
{
	model_  = model;           // obtain the model
	// --> Right now it is somewhat arbitrary which VectorPtr contains a Copy
	//     and which is a View. This should be made more explicit. Perhaps
	//     in the naming of the members...
	state_  = model_->GetState('V');  // obtain a view of the state

	rhsView_ = model_->GetRHS('V');   // view of the rhs in the model
	
	par_    = model_->GetPar();       // get the current parameter value

	// initialize Storage struct
	storage_.ds0    = ds_;
	storage_.par0   = par_;
	storage_.state0 = model_->GetState('C');
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Run()
{
	INFO("Continuation run... ");
	// set all the appropriate initial flags and counters
	step_ = 0;
	destinationReached_ = false;

	// Create the first tangent
	INFO("Create initial tangent...");
	CreateInitialTangent();
	INFO("Create initial tangent... done");

	Store(); // Store the current state and parameter
	
	while (!destinationReached_)
	{
		Step(); 	// Perform a continuation step

		// model_->DumpState();	// Get an intermediate solution

		// Detect special points:
		//  if necessary converge on them (-->not implemented yet)
		// Detect();

		++step_;
	}
	INFO("Continuation run... done");
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Step()
{
	Info(); 	                      // Print continuation info

	EulerPredictor(); 	              // Apply Euler predictor

	int status = NewtonCorrector();   // Apply Newton corrector

	// If necessary reset the step, otherwise perform a normal
	// calculation of the tangent and step adjustment
	if (status)
		Reset();
	else
	{		
		// Create new tangents based on new and stored state
		CreateTangent();			
	
		// step adjustment based on the number of Newton iterations.
		AdjustStep();

		Store();    // Store the current state and parameter			
	}
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

	// We scale the tangent of the state with its length and norm.
	// --> Not sure if this is needed, but smaller is safer?
	int length = stateDot_->Length();
	double nrm = stateDot_->Norm();	
	INFO("   normalizing with " << length << " * " << nrm);
	INFO("   before: " << nrm);
	stateDot_->Scale(1.0 / (length * nrm));
	nrm = stateDot_->Norm();	
	INFO("    after: " << nrm);

	// Set the initial parameter tangent to 1
	parDot_ = 1.0;

	INFO("  norm^2  = " << nrm * nrm << " + " << parDot_ * parDot_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::CreateTangent()
{
	// ------------------------------------------------------------
	// Calculate stateDot_:
	// --> we require the existence of a stored state and parameter
	//  - Initialize with a copy of the state in the model
	//  - Calculate finite difference for the tangent based on the
	//     stored state and the current state
	//  - Calculate finite difference for the parameter tangent
	// ------------------------------------------------------------

	// Get a copy of the current state
	stateDot_ = model_->GetState('C');

	// Get the previous state and par from our storage
	VectorPtr state0 = storage_.state0;
	double par0      = storage_.par0;

	// Compute stateDot = (state1 - state0)/ds
	// --> Note that our state_ member contains a view of the
	//     current state in the model
	stateDot_->Update(-1.0 / ds_, *state0, 1.0 / ds_);

	// Normalize stateDot
	double nrm = stateDot_->Norm();
	INFO("   normalizing with " << nrm);
	INFO("   before: " << nrm);	
	stateDot_->Scale(1.0 / nrm);
	nrm = stateDot_->Norm();
	INFO("    after " << nrm);	
	
    // Compute parDot = (par1 - par0)/ds
	par_    = model_->GetPar();             // update our own par_
	parDot_ = (par_ - par0) / ds_;
	
	INFO("  norm^2  = " << nrm * nrm << " + " << parDot_ * parDot_);
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
	INFO("Apply Euler predictor...");
	// At the end of this function the model will be
	// in a 'predicted' state.

	// Apply predictor to the state in the model
	// Compute: state = state0 + ds * statedot
	// --> 1) Note that at this point state0 and state are equal.
	// --> 2) Note that our state_ member contains a view of the
	//        current state of the model.
	state_->Update(ds_, *stateDot_, 1.0);

	// Compute  par = par0 + ds * pardot
	// --> Note that at this point par0 and par are equal.
	par_ = par_ + ds_ * parDot_;
	INFO("                      old par: " << storage_.par0);
	INFO("                predicted par: " << par_);
	INFO("               norm old state: " << storage_.state0->Norm());
	INFO("         norm predicted state: " << state_->Norm());

	// Make sure the model has the same par
	model_->SetPar(par_);
	INFO("Apply Euler predictor... done");
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
int Continuation<ModelPtr, VectorPtr>::NewtonCorrector()
{
	INFO("Apply Newton corrector...");		
	int backTrackingStatus;  // status of backtracking
	char mode;	             // mode for ComputeDFDPar()
	
	VectorPtr stateDir;      // direction for the state
	double parDir;           // direction for the parameter

	
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
		// r = ds - stateDot'*(state1 - state0)
		//           - parDot*(par2   - par0)
		VectorPtr state0 = storage_.state0;
		double    par0   = storage_.par0;
		double    r      = ds_ - stateDot_->Dot(*state_)
			                   - stateDot_->Dot(*state0)
			                   - parDot_ * (par_ - par0);

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
		parDir = (r - stateDot_->Dot(*z1)) / (parDot_ - stateDot_->Dot(*y1));

		// Then for the state:
		//  - We perform an update on y1 and
		//  - let that be the new direction
		y1->Update(1.0, *z1, -1.0 * parDir);
		stateDir = y1;

		// Update the state and the parameter in the model
		state_->Update(1.0, *stateDir, 1.0);
		par_ = par_ + parDir;  // update our parameter
		model_->SetPar(par_);  // set it in the model

		// test new RHS, if it satisfies the tolerance we exit the loop
		normRHS_ = rhsView_->Norm();
		model_->ComputeRHS();
		normRHStest_ = rhsView_->Norm();
		INFO("   NewtonCorrector      iter: " << iter_ );
		INFO("   NewtonCorrector  old norm: " << normRHS_ );
		INFO("   NewtonCorrector  new norm: " << normRHStest_ );
		if (normRHStest_ < newtonTolerance_)
			break;

		// if we don't see a decrease we run the backtracking routine
		if (backTracking_ and (normRHS_ < normRHStest_) )
			backTrackingStatus = RunBackTracking(stateDir, parDir);
	}
	INFO("  NewtonCorrector:          corrected par: " << par_);
	INFO("  NewtonCorrector:  norm  corrected state: " << state_->Norm());

	if (iter_ == newtonIterations_)
	{
		WARNING("Continuation: Newton failed!!", __FILE__, __LINE__);
		return 1;
	}

	if (par_ < storage_.par0 && forceIncreaseParam_)
	{
		WARNING("Non-increasing parameter detected", __FILE__, __LINE__);
		return 2;
	}
	
	INFO("Apply Newton corrector... done");		
	return 0;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
int Continuation<ModelPtr, VectorPtr>::RunBackTracking(VectorPtr stateDir,
													   double parDir)
{
	// Initialize reduction with -1/2
	double reduction = -1.0 / 2;
	// double reduction = -.1;

	for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_)
	{
		if (normRHStest_ < normRHS_)
			break;

		// Apply reduction to the state vector and the parameter
		state_->Update(reduction, *stateDir, 1.0); // update the state
		par_ = par_ + reduction * parDir;          // update our parameter
		model_->SetPar(par_);		               // set it in the model

		// Compute new RHS and obtain its norm
		model_->ComputeRHS();
		normRHStest_ = rhsView_->Norm();

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

	// Converge on sign switch
	if (SGN(parDot_) != parDotSign_)
	{
		INFO("Continuation::Detect():  sign switch detected");
		parDotSign_ = SGN(parDot_);
	}
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::AdjustStep()
{
	if (iter_ < minAllowedIterations_)
		ds_ = MIN(ds_ * scale1_, dsMax_);
	else if (iter_ > maxAllowedIterations_)
		ds_ = MAX(ds_ / scale1_, dsMin_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Reset()
{
	INFO("Continuation: RESET!");
	INFO("Continuation:  Restoring model...");
	Restore();
	if (ds_ > dsMin_ && 0)
	{
		ds_ = ds_ / scale2_;
		stateDot_->Scale(1.0 / scale2_);
		parDot_ = parDot_ / scale2_;
	}
	INFO("Continuation:  Restoring model... done");
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Store()
{
	// We keep two previous states, state0 and state00
	storage_.state00 = storage_.state0;
	storage_.state0  = model_->GetState('C');
	storage_.par00   = storage_.par0;
	storage_.par0    = model_->GetPar();
	storage_.ds00    = storage_.ds0;
	storage_.ds0     = ds_;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Restore()
{
 	INFO("Continuation::Restore()");
	INFO("    norm state before set: " << state_->Norm());
	INFO("       norm stored state0: " << storage_.state0->Norm());
	INFO("      norm stored state00: " << storage_.state00->Norm());
	model_->SetState(storage_.state0);

	// make sure we have a view of the state
	state_ = model_->GetState('V');
	INFO("     norm state after set: " << state_->Norm());

	INFO("    par before set: " << model_->GetPar());
	model_->SetPar(storage_.par0);
	INFO("     par after set: " << model_->GetPar());
	
	par_ = storage_.par0;
	// ds_  = storage_.ds0;

	storage_.state0  = storage_.state00;
	storage_.state00 = model_->GetState('C');
	
	storage_.par0    = storage_.par00;
	// storage_.ds0     = storage_.ds00;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr>
void Continuation<ModelPtr, VectorPtr>::Info()
{
	INFO("-----------------------------------------");
	INFO("Continuation status: ");
	INFO("            ds:  " << ds_);
	INFO("           par:  " << par_);
	INFO("        parDot:  " << parDot_);
	INFO("          step:  " << step_);
	INFO("    norm state:  " << state_->Norm());
	INFO("-----------------------------------------");
}

#endif
