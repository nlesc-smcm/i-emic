//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef CONTINUATION_H
#define CONTINUATION_H

//======================================================================
#include "ContinuationDecl.H"
#include "GlobalDefinitions.H"
#include <math.h> // pow(), sqrt()
#include <iomanip>
//======================================================================
//Constructor
template<typename Model, typename ParameterList>
Continuation<Model, ParameterList>::
Continuation(Model model, ParameterList pars)
	:
	ds_                    (pars->get("initial step size", 0.02)),
	dsMin_                 (pars->get("minimum step size", 1.0e-5)),
	dsMax_                 (pars->get("maximum step size", 1.0e01)),
	parDest_               (pars->get("parameter destination", 1.0)),
	scale1_                (pars->get("increase step size", 2.00)),
	scale2_                (pars->get("decrease step size", 2.00)),
	epsilon_               (pars->get("epsilon increment", 1.0e-8)),
	backTracking_          (pars->get("enable backtracking", true)),
	numBackTrackingSteps_  (pars->get("backtracking steps", 10)),
	maxSteps_              (pars->get("maximum number of steps", -1)),
	maxNewtonIterations_   (pars->get("maximum Newton iterations", 15)),
	minAllowedIterations_  (pars->get("minimum desired Newton iterations", 4)),
	maxAllowedIterations_  (pars->get("maximum desired Newton iterations", 5)),
	newtonTolerance_       (pars->get("Newton tolerance", 1.0e-3)),
	rejectFailedNewton_    (pars->get("reject failed iteration", true)),
	cycleTolerance_        (pars->get("enable tolerance cycling", false)),
	usePracticalTol_       (pars->get("enable practical tolerance", false))	
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

	//-->experimental
	double zetapower = pars->get("zetapower", 1.0);
	zeta_ = pow(zeta_, zetapower );

	newtonIter_ = 0;	
	backTrack_  = 0;	
	parDotSign_ = 1;
}

//======================================================================
//Destructor
template<typename Model, typename ParameterList>
Continuation<Model, ParameterList>::
~Continuation()
{}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
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
		// Let the model do some administrative work before a step
		model_->preProcess();
		
		step();  // Perform a continuation step
		
		// Let the model do some administrative work after a step 
		model_->postProcess();	   
		
		// Detect special points:
		//  if necessary converge on them (-->not implemented yet)
		detect();
		
		++step_;
	}
	INFO("---------Finished continuation run------------");
	modelInfo();
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
step()
{
	computeTolerance();             // Calculate practical tolerance 
	info(); 	                    // Print continuation info
	eulerPredictor(); 	            // Apply Euler predictor
	int status = newtonCorrector(); // Apply Newton corrector

	// If necessary reset the step, otherwise perform a normal
	// calculation of the tangent and step adjustment
	if (status)
		reset();
	else
	{
		// Put the parameter and norm of the state in the history
		parHist_.push_back(par_);
		stateNormHist_.push_back(stateView_->norm());

		// Inspect the history for weird behaviour
		analyzeHist();
		
		// Create new tangents based on result from newtonCorrector
		createTangent('E');			
	
		// step adjustment based on the number of Newton iterations.
		adjustStep();

		store();    // Store the current state (-->superfluous?)
	}
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
computeTolerance()
{
	// Scale the requested tolerance with sqrt(n)
	if (usePracticalTol_)
		practicalTol_ = newtonTolerance_ * sqrt(stateView_->length());
	else
		practicalTol_ = newtonTolerance_;

	// If requested we cycle between several multiples of practicalTol_
	if (cycleTolerance_)
	{
		int multiples[5] = {1,5,10,10,5};
		practicalTol_ *= multiples[step_ % 5];
	}
}
	
//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
createInitialTangent()
{
  INFO("Continuation: create initial tangent...");
	//-----------------------------------------------------------------
	// Initial tangent:
	// 1) Take the derivative of the RHS w.r.t. the continuation par
	//    (dFdPar) using a finite difference.
	// 2) Solve J*statedot = -dFdPar.
	//-----------------------------------------------------------------

	// 1) Compute dFdPar_, force computation of RHS
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
	INFO("Continuation: create initial tangent... done");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
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
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
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
		Vector state0 = storage_.state0;
		double par0   = storage_.par0;
		
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
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
normalize(Vector stateDot, double scale, double parDot)
{
	// -----------------------------------------------------------------
	// We normalize the tangent of the branch (stateDot, parDot)
	// 1) First we scale the tangent of the state
	// 2) Then we compute the norm of the combined vector
	// 3) Finally we normalize (stateDot, parDot)
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
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
eulerPredictor()
{
	INFO("Continuation: predictor");
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
template<typename Model, typename ParameterList>
int Continuation<Model, ParameterList>::
newtonCorrector()
{
  INFO("Continuation: Newton corrector...");
	char   mode;	      // mode for ComputeDFDPar()
	Vector stateDir;      // direction for the state
	Vector y1;            // solution of solve with dFdPar
	Vector z1;            // solution of solve with F
	double parDir;        // direction for the parameter

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
		Vector R = rhsCopy_;
		R->scale(-1.0);
		
		// Obtain the lower part (r) of the continuation RHS,
		// Here we compute
		// r = ds - stateDot'*(state1 - state0)*zeta
		//           - parDot*(par2   - par0)
		// > zeta is an additional scaling for the state tangent
		//   and makes sure that we have a proper parDot after
		//   normalization.
		Vector state0 = storage_.state0;
		double par0   = storage_.par0;
		double r      =  ds_ - stateDot_->dot(*stateView_) * zeta_
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
		normRHStest_ = rhsView_->norm('V');

		INFO("----------------------------------------------------------");
		INFO("           Newton corrector    iter: " << newtonIter_ );
		INFO("                          tolerance: " << practicalTol_);
		INFO("                              ||R||: " << normRHStest_);
		INFO("                             ||dx||: " << stateDir->norm());
		INFO("----------------------------------------------------------");		
		
		if (normRHStest_ < practicalTol_) break;
		
		// if we don't see a decrease we run the backtracking routine
		if (backTracking_ and (normRHS_ < normRHStest_) )
		{
			runBackTracking(stateDir, parDir);
			normRHStest_ = rhsView_->norm(); 
			if (normRHStest_ < practicalTol_) break;					
		}
	}
	
	// Store the last obtained y1 in stateDot, such that it can be used by
	// CreateTangent()
	stateDot_ = y1;
	
	// Do some administrative work
	TRACK_ITERATIONS("Continuation: Newton iterations...", newtonIter_);
	
	if (newtonIter_ == maxNewtonIterations_)
	  {
	    WARNING("Continuation: Newton failed after " << newtonIter_ << " steps", 
		    __FILE__, __LINE__);
	    if (rejectFailedNewton_)
	      return 1;
	    else
	      INFO("Continuation: Proceeding with unconverged result");
	  }
	else
	  INFO("Continuation: corrector converged in " << newtonIter_ << " steps");
	
	INFO("Continuation: Newton corrector... done");
	return 0;
}

//======================================================================
template<typename Model, typename ParameterList>
int Continuation<Model, ParameterList>::
runBackTracking(Vector stateDir, double parDir)
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

		INFO("    backtracking step " << backTrack_ <<
			 ", norm: " << normRHStest_ << " < " << normRHS_ << " ? ");

		reduction /= 2.0;	// Update reduction
	}

	TRACK_ITERATIONS("Continuation: backtracking steps...", backTrack_);
	
	if (backTrack_ == numBackTrackingSteps_)
	{
		WARNING("Continuation: backtracking failed ", __FILE__, __LINE__);
		return 1;
	}
	return 0;
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
detect(char mode)
{
	// Just to be on the safe side
	par_ = model_->getPar();
	double f1;
	f1 = par_ - parDest_;
	
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
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
adjustStep()
{
	// we may need something like
	//   if (newtonIter_ < 1.5 * work)
	//       increase ds
	//       work = newtonIter_
	//   else
	//       decrease ds       
	         
	if (newtonIter_ < minAllowedIterations_)
	{
		INFO("      adjusting ds... ");
		INFO("      old ds: " << ds_);
		ds_ = std::min(ds_ * scale1_, dsMax_);
		INFO("      new ds: " << ds_);
	}
	else if (newtonIter_ > maxAllowedIterations_)
	{
		INFO("      adjusting ds... ");
		INFO("      old ds: " << ds_);
		ds_ = std::max(ds_ / scale2_, dsMin_);
		INFO("      new ds: " << ds_);
	}
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
analyzeHist()
{
	if (parHist_.size() > 5)
	{
		// look at the difference between the past 2 steps
		double test = std::abs(*(parHist_.end()-1) - *(parHist_.end()-3));
		if (test < 1e-8)
		{
			INFO("************************************************");
			INFO("Parameter appears to stagnate...");
			INFO("  possible fixes:  adjust zeta");
			INFO("                   adjust ds");
			INFO("************************************************");
		}
	}
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
reset()
{
	INFO("***********************************************");
	INFO("Continuation: RESETTING!");
	INFO("              Step = " << step_);
	INFO("              Restoring model...");
	restore();
	ds_ = MAX(ds_ / scale2_, dsMin_);
	INFO("              Restoring model... done");
	INFO("***********************************************");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
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
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
restore()
{
 	INFO("   Continuation::Restore()");
	INFO("                 norm state: " << stateView_->norm());
	INFO("         norm stored state0: " << storage_.state0->norm());
	INFO("        norm stored state00: " << storage_.state00->norm());
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
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
test()
{
	model_->computeRHS();
	testCopyView();
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
testCopyView()
{
	INFO("******************************************************");
	INFO("Continuation: Test Copy / View behaviour of the model.");
	Vector view1 = model_->getState('V');
	Vector copy1 = model_->getState('C');
	Vector view2 = model_->getState('V');
	Vector copy2 = model_->getState('C');
	INFO("view1      norm: " << view1->norm());
	INFO("view2      norm: " << view2->norm());
	INFO("copy1      norm: " << copy1->norm());
	INFO("copy2      norm: " << copy2->norm());
	
	INFO("Testing view");
	INFO("Scaling view1 with 2 ---------------");
	view1->scale(2.0);
	INFO("view1      norm: " << view1->norm());
	INFO("view2      norm: " << view2->norm());
	if (view1->norm() != view2->norm())
	{
		INFO("Test failed: ||view1|| != ||view2||");
	}
	else
	{
		INFO("Test passed: ||view1|| == ||view2||");
	}
	INFO("Restore view: scale with 1/2");
	view1->scale(0.5);
	INFO("Testing copy");
	INFO("Scaling copy1 with 3 ---------------");
	copy1->scale(3.0);
	INFO("copy1      norm: " << copy1->norm());
	INFO("copy2      norm: " << copy2->norm());
	if (copy1->norm() == copy2->norm())
	{
		INFO("Test failed: ||copy1|| == ||copy2||");
	}
	else
	{
		INFO("Test passed: ||copy1|| != ||copy2||");
	}
	INFO("******************************************************");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
info()
{
	INFO(std::setprecision(8));
	INFO("-----------------------------------------");
	INFO("Continuation summary ");
	INFO("        step size:  " << ds_);
	INFO("  parameter value:  " << par_);
	INFO("      destination:  " << parDest_);
	INFO("        d/ds(par):  " << parDot_);

	if (maxSteps_ >= 0)
	{
		INFO("             step:  " << step_ << "/" << maxSteps_);
	}
	else
	{
		INFO("             step:  " << step_);
	}
	
	INFO("            ||x||:  " << stateView_->norm());
	INFO("      ||d/ds(x)||:  " << stateDot_->norm());
	INFO("     newton iters:  " << newtonIter_);
	INFO("-----------------------------------------");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
modelInfo()
{
	INFO(std::setprecision(10));
	INFO("-----------------------------------------");
	INFO("Model view: ");
	INFO("     norm state :  " << stateView_->norm());
	INFO("     norm rhs   :  " << rhsView_->norm());
	INFO("     norm sol   :  " << solView_->norm());
	INFO("    parameter   :  " << model_->getPar());
	INFO("-----------------------------------------");
}

#endif
