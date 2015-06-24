//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef CONTINUATION_H
#define CONTINUATION_H

//======================================================================
#include "ContinuationDecl.H"
#include "GlobalDefinitions.H"
#include <math.h>
//======================================================================
//Constructor
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
Continuation<ModelPtr, VectorPtr, ParameterList>::
Continuation(ModelPtr model, ParameterList pars)
	:

	ds_(pars->get("ds",0.02)),
	dsMin_(pars->get("dsMin", 1.0e-5)),
	dsMax_(pars->get("dsMax", 1.0e01)),
	maxSteps_(pars->get("maxSteps", -1)),

	scale1_(pars->get("scale1", 1.00)),
	scale2_(pars->get("scale2", 2.00)),

	epsilon_(pars->get("epsilon", 1.0e-8)),

	backTracking_(pars->get("backTracking", true)),
	numBackTrackingSteps_(pars->get("numBackTrackingSteps", 10)),

	maxNewtonIterations_(pars->get("maxNewtonIterations", 15)),
	minAllowedIterations_(pars->get("minAllowedIterations", 4)),
	maxAllowedIterations_(pars->get("maxAllowedIterations", 5)),
	newtonTolerance_(pars->get("newtonTolerance", 1.0e-3))
{
	model_  = model;

	// Get views of the important components of the model
	stateView_ = model_->getState('V');
	rhsView_   = model_->getRHS('V');
	solView_   = model_->getSolution('V');
	
 	par_       = model_->getPar();       

	// print some info on the model
	modelInfo();
	
	// initialize Storage struct
	storage_.ds0    = ds_;
	storage_.par0   = par_;
	storage_.state0 = model_->getState('C');

	// scaling
	zeta_ = 1.0 / stateView_->length();	

	newtonIter_ = 0;	
	backTrack_  = 0;	
	parDotSign_ = 1;	
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
run()
{
	// set all the appropriate initial flags and counters
	step_ = 0;
	destinationReached_ = false;

	// Create the first tangent
	createInitialTangent();

	store(); // Store the current state and parameter
	
	while (!destinationReached_ && step_ != maxSteps_)
	{
		step(); 	// Perform a continuation step

		model_->dumpState();	// Get an intermediate solution

		// Detect special points:
		//  if necessary converge on them (-->not implemented yet)
		detect();

		++step_;
	}
	INFO("---------Finished continuation run------------");
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
step()
{
	info(); 	                      // Print continuation info

	eulerPredictor(); 	              // Apply Euler predictor

	int status = newtonCorrector();   // Apply Newton corrector

	// If necessary reset the step, otherwise perform a normal
	// calculation of the tangent and step adjustment
	if (status)
	{
		reset();
	}
	else
	{		
		// Create new tangents based on result from newtonCorrector
		createTangent('E');			
	
		// step adjustment based on the number of Newton iterations.
		adjustStep();

		store();    // Store the current state and parameter			
	}
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
createInitialTangent()
{
	//-----------------------------------------------------------------
	// Initial tangent:
	// 1) Take the derivative of the RHS w.r.t. the continuation par
	//    (dFdPar) using a finite difference.
	// 2) Solve J*statedot = -dFdPar.
	//-----------------------------------------------------------------

	// 1) Compute dFdPar_, force compute of RHS
	computeDFDPar('F');

	// 2) Solve J*stateDot_ = -dFdPar_
	model_->computeJacobian();
	dFdPar_->scale(-1.0);
	model_->solve(dFdPar_);

	// To obtain the solution from this solve without the risk
	// of altering it we get a copy.
	stateDot_ = model_->getSolution('C');

	// We scale the tangent of the state and parameter such they are
	// normalized.
	normalize(stateDot_, zeta_, parDot_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
computeDFDPar(char mode)
{
	// Let the model compute an RHS F(par) with its current state
	// - Note that this is usually a predicted state
	// - The RHS for this particular state may already have been
	//   computed.
	//
	if (mode == 'F')
		model_->computeRHS();
	
	// Get a copy of this RHS, store it in our rhsCopy_ member
	rhsCopy_ = model_->getRHS('C');

	// Calculate new RHS
	model_->setPar(par_ + epsilon_);  // increment parameter --> par + eps
	model_->computeRHS();             // compute new RHS     --> F(par+eps)
	model_->setPar(par_);             // restore parameter   --> par

	// Get a copy of the new RHS F(par+eps), put it in dFdPar_
	dFdPar_ = model_->getRHS('C');

	// Update dFdPar_ such that it contains (F(par+eps) - F(par))/eps
	dFdPar_->update(-1.0 / epsilon_, *rhsCopy_, 1.0 / epsilon_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
createTangent(char mode)
{
	if (mode == 'S') // Secant
	{
		// ------------------------------------------------------------
		// Calculate stateDot_ with the Secant approach:
		//  - Initialize with a copy of the state in the model
		//  - Calculate finite difference for the tangent based on the
		//     stored state and the current state
		//  - Calculate finite difference for the parameter tangent
		// ------------------------------------------------------------
		
		// Get a copy of the current state
		stateDot_ = model_->getState('C');
		
		// Get the previous state and par from our storage
		VectorPtr state0 = storage_.state0;
		double par0      = storage_.par0;
		
		// Compute stateDot = (state1 - state0)/ds
		stateDot_->update(-1.0 / ds_, *state0, 1.0 / ds_);
		
		// Compute parDot = (par1 - par0)/ds
		par_    = model_->getPar();      // update our par_
		parDot_ = (par_ - par0) / ds_;

		// Normalize the combined tangent
		// --> Not sure about zeta_ here
		normalize(stateDot_, zeta_, parDot_);
	}
	else if (mode == 'E' && newtonIter_ != 0) // Euler 
	{
		// ------------------------------------------------------------
		// Calculate stateDot_ with the Euler approach. Here we require
		// that the corrector has been run at least once, so:
		//   - At this point stateDot_ contains y1 from J*y1 = dFdPar
		//   - This means we can flip its sign to create a new tangent
		// ------------------------------------------------------------
		stateDot_->scale(-1.0);

		// We still need to make sure everything is normalized
		normalize(stateDot_, zeta_, parDot_);
	}
	else
		WARNING("(Continuation::createTangent) Invalid mode!",
				__FILE__, __LINE__);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
normalize(VectorPtr stateDot, double scale, double parDot)
{
	// -----------------------------------------------------------------
	// We normalize the tangent of the branch (stateDot, parDot)
	// 1) First we scale the tangent of the state
	// 2) Then we compute the norm of the combined vector
	// 3) Finally we normalize both stateDot and parDot
	// -----------------------------------------------------------------

	// 1) scale stateDot
	stateDot_->scale(zeta_);

	// 2) compute norm of combined vector
	double nrm      = stateDot_->norm();
	double normComb = sqrt(nrm * nrm + 1);

	// 3) scale the components
	stateDot_->scale(1.0 / normComb);
	parDot_ = 1.0 / normComb;

	DEBUG(" Normalize():");
	DEBUG("  norm  = " << stateDot_->norm() * stateDot_->norm()
		                  + parDot_ * parDot_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
eulerPredictor()
{	
	// At the end of this function the model will be
	// in a 'predicted' state.

	// Apply predictor to the state in the model
	// Compute: state = state0 + ds * statedot
	//  - Note that at this point state0 and state are equal.
	stateView_->update(ds_, *stateDot_, 1.0);

	// Compute  par = par0 + ds * pardot
	// - Note that at this point par0 and par are equal.
	par_ = par_ + ds_ * parDot_;
	
	DEBUG("                      old par: " << storage_.par0);
	DEBUG("                predicted par: " << par_);
	DEBUG("               norm old state: " << storage_.state0->norm());
	DEBUG("         norm predicted state: " << stateView_->norm());

	// Make sure the model has the same par
	model_->setPar(par_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
int Continuation<ModelPtr, VectorPtr, ParameterList>::
newtonCorrector()
{
	int backTrackingStatus;  // status of backtracking
	char mode;	             // mode for ComputeDFDPar()
	
	VectorPtr stateDir;      // direction for the state
	VectorPtr y1;            // solution of solve with dFdPar
	VectorPtr z1;            // solution of solve with F
	double parDir;           // direction for the parameter
	newtonIter_ = 0;	
	while (newtonIter_ < maxNewtonIterations_)
	{
		// Taking the derivative of the RHS w.r.t. the continuation
		// parameter using a finite difference. In the first iteration
		// the computation of the RHS is required.
		mode = (newtonIter_ == 0) ? 'F' : 'A'; 
		computeDFDPar(mode);

		// Obtain the upper part (R) of the continuation RHS.
		// A copy of F(par) is obtained in ComputeDFDPar(), so
		// we use that one to get -F(par).
		// > Note that we cannot use the current rhs in the model
		//   since at this point it corresponds to F(par+eps)
		VectorPtr R = rhsCopy_;
		R->scale(-1.0);
		
		// Obtain the lower part (r) of the continuation RHS,
		// Here we compute
		// r = ds - stateDot'*(state1 - state0)*zeta
		//           - parDot*(par2   - par0)
		// > zeta is an additional scaling for the state tangent
		//   and makes sure that we have a proper parDot after
		//   normalization.
		VectorPtr state0 = storage_.state0;
		double    par0   = storage_.par0;
		double    r      =  ds_ - stateDot_->dot(*stateView_) * zeta_
			                - stateDot_->dot(*state0) * zeta_
			                - parDot_ * (par_ - par0);

		// At this point the model contains the predicted state and
		// parameter. The Jacobian will be computed based on the
		// predicted data.
		model_->computeJacobian();

		// Now we will perform 2 solves to solve the bordered system:
		// In both cases we obtain copies of the solution. Both copies
		// wil have their use either here or in the computation of the
		// next tangent.
		model_->solve(dFdPar_);
		y1 = model_->getSolution('C');
		INFO("... about to do the standard newton solve ...");
		model_->solve(R);
		z1 = model_->getSolution('C');

		// Determine the directions
		// First for the parameter:
		parDir = (r - zeta_ * stateDot_->dot(*z1))
			   / (parDot_ - zeta_ * stateDot_->dot(*y1));

		// --> Perhaps we can optimize the next bit (copies and views)		
		// Then for the state:
		//  - We perform an update on z1 and
		//  - let that be the new direction
		z1->update(-1.0 * parDir, *y1, 1.0);
		stateDir = z1;

		// Update the state and the parameter in the model
		stateView_->update(1.0, *stateDir, 1.0);
		par_ = par_ + parDir;  // update our parameter
		model_->setPar(par_);  // set it in the model
		++newtonIter_;
		
		// test new RHS, if it satisfies the tolerance we exit the loop
		normRHS_ = rhsView_->norm();
		model_->computeRHS();
		normRHStest_ = rhsView_->norm();
		INFO("Contin.: corrector    iter: " << newtonIter_ );
		INFO("Contin.: corrector  norm R: " << normRHStest_ );
				
		if (normRHStest_ < newtonTolerance_)
		{
			break;
		}
		
		// if we don't see a decrease we run the backtracking routine
		if (backTracking_ and (normRHS_ < normRHStest_) )
			backTrackingStatus = runBackTracking(stateDir, parDir);
	}
	
	// Store the last obtained y1 in stateDot, such that it can be used by
	// CreateTangent()
	stateDot_ = y1;
	
	if (newtonIter_ == maxNewtonIterations_)
	{
		WARNING("Continuation: Newton failed", __FILE__, __LINE__);
		return 1;
	}

	INFO("Contin.: corrector converged in " << newtonIter_ << " steps");
	return 0;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
int Continuation<ModelPtr, VectorPtr, ParameterList>::
runBackTracking(VectorPtr stateDir, double parDir)
{
	// Initialize reduction with -1/2
	double reduction = -1.0 / 2;
	// double reduction = -.1;

	for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_)
	{
		if (normRHStest_ < normRHS_)
			break;

		// Apply reduction to the model state and its parameter
		stateView_->update(reduction, *stateDir, 1.0); // update the state
		par_ = par_ + reduction * parDir;          // update our parameter
		model_->setPar(par_);		               // set it in the model

		// Compute new RHS and obtain its norm
		model_->computeRHS();
		normRHStest_ = rhsView_->norm();

		INFO("Contin.: backtracking " << backTrack_ <<
			 " norm: " << normRHStest_ << " < " << normRHS_ << " ? ");

		reduction /= 2.0;	// Update reduction
	}
	if (backTrack_ == numBackTrackingSteps_)
	{
		WARNING("Continuation: backtracking failed ", __FILE__, __LINE__);
		return 1;
	}
	return 0;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
detect(char mode)
{
	// Just to be on the safe side
	par_ = model_->getPar();
	double f1;
	f1 = par_ - model_->getParDestination();

	// This is a naive initial implementation
	if (f1 > 0)
		destinationReached_ = true;

	// Converge on sign switch
	if (SGN(parDot_) != parDotSign_)
	{
		INFO("Continuation::detect():  sign switch detected");
		parDotSign_ = SGN(parDot_);
	}
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
adjustStep()
{
	if (newtonIter_ < minAllowedIterations_)
		ds_ = MIN(ds_ * scale1_, dsMax_);
	else if (newtonIter_ > maxAllowedIterations_)
		ds_ = MAX(ds_ / scale1_, dsMin_);
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
reset()
{
	INFO("***********************************************");
	INFO("Continuation: RESETTING!");
	INFO("Continuation:  Step = " << step_);
	INFO("Continuation:  Restoring model...");
	restore();
	ds_ = MAX(ds_ / scale2_, dsMin_);
	INFO("Continuation:  Restoring model... done");
	INFO("***********************************************");
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
store()
{
	// We keep two previous states, state0 and state00
	storage_.state00 = storage_.state0;
	storage_.state0  = model_->getState('C');
	storage_.par00   = storage_.par0;
	storage_.par0    = model_->getPar();
	storage_.ds00    = storage_.ds0;
	storage_.ds0     = ds_;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
restore()
{
 	INFO("Continuation::Restore()");
	INFO("               norm state: " << stateView_->norm());
	INFO("       norm stored state0: " << storage_.state0->norm());
	INFO("      norm stored state00: " << storage_.state00->norm());
	model_->setState(storage_.state0);

	// make sure we have a view of the state
	stateView_ = model_->getState('V');
	model_->setPar(storage_.par0);
	
	par_ = storage_.par0;
	ds_  = storage_.ds0;

	storage_.state0  = storage_.state00;
	storage_.state00 = model_->getState('C');
	
	storage_.par0    = storage_.par00;
	storage_.ds0     = storage_.ds00;
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
info()
{
	INFO("-----------------------------------------");
	INFO("Continuation summary: ");
	INFO("               ds:  " << ds_);
	INFO("              par:  " << par_);
	INFO("           parDot:  " << parDot_);
	if (maxSteps_ >= 0)
	{
		INFO("             step:  " << step_ << "/" << maxSteps_);
	}
	else
	{
		INFO("             step:  " << step_);
	}
	INFO("       norm state:  " << stateView_->norm());
	INFO("     newton iters:  " << newtonIter_);
	INFO("-----------------------------------------");
}

//======================================================================
template<typename ModelPtr, typename VectorPtr, typename ParameterList>
void Continuation<ModelPtr, VectorPtr, ParameterList>::
modelInfo()
{
	INFO("-----------------------------------------");
	INFO("Model view: ");
	INFO("     norm state :  " << stateView_->norm());
	INFO("     norm rhs   :  " << rhsView_->norm());
	INFO("     norm sol   :  " << solView_->norm());
	INFO("-----------------------------------------");
}

#endif
