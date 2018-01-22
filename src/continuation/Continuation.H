//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef CONTINUATION_H
#define CONTINUATION_H

//======================================================================
#include "ContinuationDecl.H"
#include "GlobalDefinitions.H"
#include "Utils.H"

#include <math.h> // pow(), sqrt()
#include <ctime>
#include <iomanip>

//======================================================================
//Constructor
template<typename Model, typename ParameterList>
Continuation<Model, ParameterList>::
Continuation(Model model, ParameterList pars)
    :
    model_(model),
    dsInit_                (pars->get("initial step size", 0.01)),
    dsMin_                 (pars->get("minimum step size", 1.0e-5)),
    dsMax_                 (pars->get("maximum step size", 1.0e-1)),
    scale1_                (pars->get("increase step size", 1.25)),
    scale2_                (pars->get("decrease step size", 1.25)),
    epsilon_               (pars->get("epsilon increment", 1.0e-4)),
    backTracking_          (pars->get("enable backtracking", false)),
    numBackTrackingSteps_  (pars->get("backtracking steps", 10)),
    backTrackIncrease_     (pars->get("backtracking increase", 10.0)),
    maxSteps_              (pars->get("maximum number of steps", -1)),
    maxNewtonIterations_   (pars->get("maximum Newton iterations", 15)),
    minNewtonIterations_   (pars->get("minimum Newton iterations", 1)),
    minAllowedIterations_  (pars->get("minimum desired Newton iterations", 4)),
    maxAllowedIterations_  (pars->get("maximum desired Newton iterations", 5)),
    newtonTolerance_       (pars->get("Newton tolerance", 1.0e-4)),
    destinationTolerance_  (pars->get("destination tolerance", 1.0e-4)),
    userDetect_            (pars->get("enable custom monitor", false)),
    tanScaling_            (pars->get("state tangent scaling", 1.0e0)),
    normalizeStrategy_     (pars->get("normalize strategy", 'N')),

    rejectFailedNewton_    (pars->get("reject failed iteration", true)),
    giveUpAtdsMin_         (pars->get("give up at minimum step size", false)),
    cycleTolerance_        (pars->get("enable tolerance cycling", false)),
    usePracticalTol_       (pars->get("enable practical tolerance", false)),
    newtChordHybr_         (pars->get("enable Newton Chord hybrid solve", false)),
    tangentType_           (pars->get("tangent type", 'S')),
    residualTest_          (pars->get("corrector residual test", 'D')),
    initialTangent_        (pars->get("initial tangent type", 'E')),
    printImportantVectors_ (pars->get("print important vectors", false)),
    predictorBound_        (pars->get("predictor bound", 1e5)),
    signMonitor_           (0)
{
    // Set the step size
    ds_ = dsInit_;

    // Get the parameter destinations
    // We allow for convergence with destination tolerance at multiple destinations
    double dest;
    std::stringstream destID;
    for (int i = 0; i != 999; ++i)
    {
        destID << "destination " << i;
        dest = pars->get(destID.str().c_str(), -999.0);
        if (std::abs(dest + 999.) < 1e-7)
            break;

        destinations_.push_back(dest);
        destID.str("");
        destID.clear();
    }
    if (destinations_.empty())
        ERROR("No destinations given in continuation xml!",
              __FILE__, __LINE__);

    INFO(std::endl << "Continuation parameter list: "
         << std::endl << *pars << std::endl);
}

//======================================================================
//Destructor
template<typename Model, typename ParameterList>
Continuation<Model, ParameterList>::
~Continuation()
{
    INFO("Continuation destructor");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
initialize()
{
    // (Re)Initialize step size
    ds_ = dsInit_;

    // Compute rhs of the model
    model_->computeRHS();

    // Get views of the important components of the model
    stateView_   = model_->getState('V');
    rhsView_     = model_->getRHS('V');
    solView_     = model_->getSolution('V');
    par_         = model_->getPar();
    startingPar_ = model_->getPar();
    
    // Save the residual
    if (printImportantVectors_)
        Utils::save(rhsView_, "residual");

    // print some info on the model
    modelInfo();

    // initialize Storage struct
    storage_.ds0    = ds_;
    storage_.par0   = par_;
    storage_.state0 = model_->getState('C');

    // initializations for detect()
    signMonitor_ = 0;
    secant_      = false;
    activeDest_  = -999; // hopefully this is never an actual dest^^

    // scaling
    if (normalizeStrategy_ == 'O')
        zeta_ = 1.0 / stateView_->GlobalLength();
    else if (normalizeStrategy_ == 'N')
        zeta_ = tanScaling_ / stateView_->GlobalLength();
    else
    {
        WARNING(" undefined normalization strategy!",__FILE__, __LINE__);
    }

    newtonIter_    = 0;
    sumNewtonIter_ = 0;
    backTrack_     = 0;
    parDotSign_    = 1;

    // initializations for main loop
    step_               = 0;
    resetCounter_       = 0;
    destinationReached_ = false;
    abortFlag_          = false;
    fixStepSize_        = false;

    if (userDetect_) { INFO("Continuation: custom monitor enabled");}
    else { INFO("Continuation: custom monitor disabled");}
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
run()
{
    INFO("Continuation::run initialize...");
    initialize();    // core initializations
    INFO("Continuation::run initialize... done");
    
    TIMER_START("Continuation: run");

    createInitialTangent(); // Create the first tangent

    detect();     // Detect special points
    userDetect(); // Use additional targets provided by model class

    int status = 0;    

    while (!destinationReached_ && (step_ != maxSteps_) && !abortFlag_)
    {
        ++step_;          // Increment step counter
        info();           // Print continuation summary
        store();          // Store the current state and parameter
        status = step();  // Perform a continuation step
        
        if (status) // If step failed we reset and continue
        {
            reset();                       // Reset the model
            continue;                      // Skip the rest
        }

        detect();     // Detect special points
        userDetect(); // Use additional targets provided by model class
        adjustStep(); // step size adjustment
    }
    info();
    modelInfo();
    TIMER_STOP("Continuation: run");
    INFO("---------Finished continuation run--------------");
}

//======================================================================
template<typename Model, typename ParameterList>
int Continuation<Model, ParameterList>::
step()
{
    TIMER_START("Continuation: step");
    // Let the model do some administrative work at the beginning of a step
    model_->preProcess();

    computeTolerance();         // Calculate practical tolerance

    int status = 0;
    status = eulerPredictor();  // Apply Euler predictor

    // If necessary reset the step, otherwise perform a normal
    // calculation of the tangent and step adjustment
    if (status)   // Failure
    {
        TIMER_STOP("Continuation: step");
        return 1; // Consider this step failed
    }

    status = newtonCorrector(); // Apply Newton corrector
    
    if (status)   // Failure
    {
        TIMER_STOP("Continuation: step");
        return 1; // Consider this step failed
    }

    if (printImportantVectors_)
    {
        // Save the residual and last solution from linear solve to binaries
        Utils::save(rhsView_, "residual");
        Utils::save(model_->getSolution('V'), "solution");
    }

    // Put the parameter and norm of the state in the history
    parHist_.push_back(par_);
    stateNormHist_.push_back(Utils::norm(stateView_));

    // Inspect the history for weird behaviour
    analyzeHist();

    // Create new tangents based on result from newtonCorrector
    createTangent(tangentType_);

    // Let the model do some administrative work at the end of a succesful step
    model_->postProcess();

    // Write continuation data to file
    bool describe = (step_ == 1) ? true : false;
    writeData(describe);

    TIMER_STOP("Continuation: step");
    return 0; // Exiting normally
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
computeTolerance()
{
    // Scale the requested tolerance with sqrt(n)
    if (usePracticalTol_)
        practicalTol_ = newtonTolerance_ * sqrt(stateView_->GlobalLength());
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

    if (printImportantVectors_) Utils::save(dFdPar_, "dFdpar");

    if ( (initialTangent_ == 'E') || (initialTangent_ == 'S'))   // 2a) Solve J*dxdpar = -dFdPar_
    {
        model_->preProcess();
        model_->computeJacobian();
        dFdPar_->Scale(-1.0);
        model_->solve(dFdPar_);
        stateDot_ = model_->getSolution('C');
    }
    else if (initialTangent_ == 'A') // 2b) Do stateDot_ = -dFdPar_
    {
        stateDot_ = model_->getSolution('C');
        dFdPar_->Scale(-1.0);
        *stateDot_ = *dFdPar_;
    }
    else
    {
        WARNING(" initialTangent invalid!" , __FILE__, __LINE__);
    }

    if (printImportantVectors_) // Print initial tangent
    {
        std::stringstream fname;
        fname << "init_tan_" << std::setprecision(4) << par_;
        Utils::print(stateDot_, fname.str());
    }

    // We scale the tangent of the state and parameter such they are
    // normalized.
    normalize();

    INFO("                 ||state|| = " << Utils::norm(stateView_));
    INFO("              ||stateDot|| = " << Utils::norm(stateDot_));
    INFO("                  parDot   = " << parDot_);
    INFO("                      ds   = " << ds_);
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
    if (mode == 'F')
        model_->computeRHS();

    INFO("    creating dFdPar...");
    INFO("       |                       F(x,l) = " << Utils::norm(rhsView_));
    // Get a copy of this RHS, store it in our rhsCopy_ member
    rhsCopy_ = model_->getRHS('C');

    // Calculate new RHS
    model_->setPar(par_ + epsilon_);  // increment parameter --> par + eps
    model_->computeRHS();             // compute new RHS     --> F(par+eps)
    INFO("       |                            l = " << par_);
    INFO("       |                          eps = " << epsilon_);
    INFO("       |                  F(x, l+eps) = " << Utils::norm(rhsView_));
    model_->setPar(par_);             // restore parameter   --> par

    // Get a copy of the new RHS F(par+eps), put it in dFdPar_
    dFdPar_ = model_->getRHS('C');

    // Update dFdPar_ such that it contains (F(par+eps) - F(par))/eps
    dFdPar_->Update(-1.0 / epsilon_, *rhsCopy_, 1.0 / epsilon_);
    INFO("       |  (F(x,l+eps) - F(x,l)) / eps = " << Utils::norm(dFdPar_));
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
createTangent(char mode)
{
    INFO("Continuation: create tangent " << mode);
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

        // Get the previous state, par and ds from our storage
        // These are used to compute a secant tangent
        double    par0   = storage_.par0;
        double    ds0    = storage_.ds0;

        // Compute stateDot = (state1 - state0)/ds0
        stateDot_->Update(-1.0 / ds0, *storage_.state0, 1.0 / ds0);
        
        // Compute parDot = (par1 - par0)/ds
        par_    = model_->getPar();      // update our par_
        parDot_ = (par_ - par0) / ds0;
        
    }
    else if (mode == 'E')  // Euler
    {
        // Calculate stateDot_ with the Euler approach
        if (newtChordHybr_)
        {
            // 1) Compute dFdPar_, force computation of RHS
            computeDFDPar('F');

            // 2) Solve J*stateDot_ = -dFdPar_
            model_->computeJacobian();
            dFdPar_->Scale(-1.0);
            model_->solve(dFdPar_);

            // To obtain the solution from this solve without the risk
            // of altering it we get a copy.
            stateDot_ = model_->getSolution('C');

        }
        else if (newtonIter_ != 0)
        {
            // ------------------------------------------------------------
            //  If the corrector has been run we have
            //  stateDot = y from J*y = dFdPar.
            //  This means we can flip its sign to create a new tangent.
            // ------------------------------------------------------------
            stateDot_->Scale(-1.0);
        }
        else
        {
            WARNING(" undefined behaviour!",__FILE__, __LINE__);
        }

        // Adjust (normalize) the Euler tangent
        normalize();
    }
    else
    {
        WARNING(" invalid mode!",__FILE__, __LINE__);
    }

    INFO("                ||state0|| = " << Utils::norm(storage_.state0));
    INFO("                 ||state|| = " << Utils::norm(stateView_));
    INFO("              ||stateDot|| = " << Utils::norm(stateDot_));
    INFO("                  parDot   = " << parDot_);
    INFO("                      ds0  = " << storage_.ds0);
    INFO("                      ds   = " << ds_);
    INFO("Continuation: create tangent " << mode << " done");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
normalize()
{
    if (normalizeStrategy_ == 'O')
    {
        // -----------------------------------------------------------------
        // We normalize the tangent of the branch (stateDot, parDot)
        // 1) First we scale the tangent of the state
        // 2) Then we compute the norm of the combined vector
        // 3) Finally we normalize (stateDot, parDot)
        // -----------------------------------------------------------------

        // 1) scale stateDot
        zeta_ = tanScaling_ / Utils::norm(stateDot_);
        stateDot_->Scale(zeta_);

        // 2) compute norm of combined vector
        double nrm      = Utils::norm(stateDot_);
        double normComb = sqrt(nrm * nrm + 1);

        // 3) scale the components
        stateDot_->Scale(1.0 / normComb);
        parDot_ = 1.0 / normComb;
    }
    else if (normalizeStrategy_ == 'N')
    {
        // At this point we assume statedot contains a derivative with respect to
        // the parameter

        // compute ||dx/dpar||
        double nrm = Utils::norm(stateDot_);

        // compute sqrt(1+zeta*||dx/dpar||^2)
        double normComb = sqrt(zeta_ * nrm * nrm + 1);

        // dpar/ds = 1 / sqrt(1 + zeta*||dx/dpar||^2)
        parDot_ = 1.0 / normComb;

        // dx/ds = dx/dpar * dpar/ds
        stateDot_->Scale(parDot_);
    }
    else
    {
        WARNING(" undefined normalization strategy!",__FILE__, __LINE__);
    }
    
}

//======================================================================
template<typename Model, typename ParameterList>
int Continuation<Model, ParameterList>::
eulerPredictor()
{
    INFO("Continuation: predictor");
    // At the end of this function the model will be
    // in a 'predicted' state.

    // Apply predictor to the state in the model
    // Compute: state = state0 + ds * statedot
    //  - Note that at this point state0 and state are equal.
    stateView_->Update(ds_, *stateDot_, 1.0);

    // Compute  par = par0 + ds * pardot
    // - Note that at this point par0 and par are equal.
    par_ = par_ + ds_ * parDot_;

    INFO("   |                   old par: " << storage_.par0);
    INFO("   |             predicted par: " << par_);
    INFO("   |            norm old state: " << Utils::norm(storage_.state0));
    INFO("   |      norm predicted state: " << Utils::norm(stateView_));

    // Make sure the model has the same par
    model_->setPar(par_);

    // Test whether our prediction is acceptable
    model_->computeRHS();
    double rhsNrm = Utils::norm(rhsView_);
    INFO("   |                  norm rhs: " << rhsNrm);
    if (rhsNrm > predictorBound_)
    {
        INFO("   |      too big!");
        return 1;
    }
    else
        return 0;
}

//======================================================================
template<typename Model, typename ParameterList>
int Continuation<Model, ParameterList>::
newtonCorrector()
{
    INFO("Continuation: Newton corrector...");
    char      mode;          // mode for ComputeDFDPar()
    VectorPtr stateDir;      // direction for the state
    VectorPtr y;             // solution of solve with dFdPar
    VectorPtr z;             // solution of solve with F
    double    parDir;        // direction for the parameter

    double res0   = 100.0;
    double res    = 100.0;
    double normDX = 100.0;

    newtonIter_ = 0;
    while ( newtonIter_ < maxNewtonIterations_ )
    {
        // save the old residual
        res0 = res;

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
        normRHS_    = Utils::norm(rhsCopy_);

        R->Scale(-1.0);

        // Obtain the lower part (rbp in bag.f) of the continuation RHS,
        // (state1 - state0)
        VectorPtr stateDiff = model_->getState('C');
        stateDiff->Update(-1.0, *storage_.state0, 1.0);

        // (par2   - par0)
        double parDiff = par_ - storage_.par0;

        // Create normalization constraint
        double rbp;
        if (normalizeStrategy_ == 'O')
        {
            // rbp = ds - d/ds state^T * (state1 - state0) * zeta
            //             - d/ds par * (par1   - par0)
            rbp = ds_ - Utils::dot(stateDot_, stateDiff) * zeta_
                - parDot_ * parDiff;
        }
        else if (normalizeStrategy_ == 'N')
        {
            // rbp = ds*ds - (state1 - state0)^T * (state1 - state0) * zeta
            //             - (par2   - par0)^2
            rbp = (ds_ * ds_) - Utils::dot(stateDiff, stateDiff) * zeta_
                - (parDiff  * parDiff);
        }
        else
        {
            WARNING(" undefined normalization strategy!",__FILE__, __LINE__);
        }

        // At this point the model contains the predicted state and
        // parameter. The Jacobian will be computed based on the
        // predicted data.
        model_->computeJacobian();

        // Now we will perform 2 solves to solve the bordered system:
        // In both cases we obtain copies of the solution. Both copies
        // wil have their use either here or in the computation of the
        // next tangent.
        if (!newtChordHybr_)
        {
            model_->solve(dFdPar_);
            y = model_->getSolution('C');
        }

        model_->solve(R);
        z = model_->getSolution('C');

        // Determine the directions.....................................
        // First for the parameter:
        if (normalizeStrategy_ == 'O')
        {
            if (newtChordHybr_)
                parDir = (rbp - zeta_ * Utils::dot(stateDot_, z))
                    / (parDot_ + zeta_ * Utils::dot(stateDot_, stateDot_));
            else
                parDir = (rbp - zeta_ * Utils::dot(stateDot_, z))
                    / (parDot_ - zeta_ * Utils::dot(stateDot_, y));
            
        }
        else if (normalizeStrategy_ == 'N')
        {
            if (newtChordHybr_)
                parDir = (rbp - 2 * zeta_ * Utils::dot(stateDiff, z))
                    / (2 * parDiff + 2 * (zeta_ / parDiff) * Utils::dot(stateDiff, stateDiff));
            else
                parDir = (rbp - 2 * zeta_ * Utils::dot(stateDiff, z))
                    / (2 * parDiff - 2 * zeta_ * Utils::dot(stateDiff, y));
        }
        else
        {
            WARNING(" undefined normalization strategy!",__FILE__, __LINE__);
        }

        // Then for the state:
        //  we perform an update on z
        if (newtChordHybr_)
            z->Update(1.0 * parDir, *stateDot_, 1.0);
        else
            z->Update(-1.0 * parDir, *y, 1.0);

        //  let that be the new direction
        stateDir = z;
        
        // Update the state and the parameter in the model
        stateView_->Update(1.0, *stateDir, 1.0);
        par_ = par_ + parDir;  // update our parameter
        model_->setPar(par_);  // set it in the model
        
        ++newtonIter_;
        ++sumNewtonIter_;

        // test new RHS, if it satisfies the tolerance we exit the loop
        model_->computeRHS();
        normRHStest_ = Utils::norm(rhsView_);

        // if we see a drastic increase in |rhs| we quit the corrector
        // again, use predictor bound to select failures
        if (normRHStest_ > predictorBound_)
        {
            INFO(" norm too big! " << normRHStest_);
            return 1;
        }

        // if we don't see a decrease we run the backtracking routine
        if (backTracking_ and (normRHS_ < normRHStest_) )
        {
            int status = runBackTracking(stateDir, parDir);
            if (status)
                return 1; // exiting Newton when backtracking fails
        }

        // if we see a drastic increase in |dx| we also quit the corrector
        if (Utils::norm(stateDir) > 1e3 * Utils::norm(storage_.state0) &&
            Utils::norm(storage_.state0) > 0)
        {
            WARNING("  |dx| = " << Utils::norm(stateDir)
                    << " >> old |x| = " << Utils::norm(storage_.state0),
                    __FILE__, __LINE__);
            return 1;
        }

        // residual tests and printing --------------------------------------------
        
        // obtain the residual
        if (residualTest_ == 'R')
        {
            res = normRHStest_;
        }
        else if (residualTest_ == 'D')
        {
            normDX = Utils::normInf(stateDir);                              
            res = std::max(std::abs(parDir), normDX);
        }
        else
        {
            WARNING(" undefined residual test!",__FILE__, __LINE__);
            res = 999;
        }

        INFO("----------------------------------------------------------");
        INFO("           Newton corrector    iter: " << newtonIter_ );
        INFO("                       old ||R||_2 : " << normRHS_);
        
        if (residualTest_ == 'R')
        {
            INFO("                           ||R||_2 : " << res << " <? " << practicalTol_);
        }
        else if (residualTest_ == 'D')
        {
            INFO("                           ||R||_2 : " << normRHStest_);                   
            INFO("                     ||dx,dl||_inf : " << res << " <? " << practicalTol_);
        }
        
        INFO("                     old res / res : " << res0 / res);
        INFO("                          ||dx||_2 : " << Utils::norm(stateDir));
        INFO("                           ||x||_2 : " << Utils::norm(stateView_));
        INFO("                       old ||x||_2 : " << Utils::norm(storage_.state0));
        INFO("                               dl  : " << parDir);
        INFO("                                l  : " << par_);
        INFO("                            old l  : " << storage_.par0);
        INFO("----------------------------------------------------------");

        if (printImportantVectors_)
            Utils::save(stateDir, "dx");

        if (res < practicalTol_ && newtonIter_ >= minNewtonIterations_)
            break;
    }

    // Keep the last obtained y in stateDot, such that it can be used by
    // CreateTangent()
    if (!newtChordHybr_) stateDot_ = y;

    // Do some administrative work
    TRACK_ITERATIONS("Continuation: Newton iterations...", newtonIter_);

    if (res > practicalTol_)
    {
        WARNING("Continuation: Newton failed after " << newtonIter_ << " steps",
                __FILE__, __LINE__);
        if (rejectFailedNewton_)
            return 1;
        else
        {
            WARNING("Continuation: Proceeding with unconverged result",__FILE__, __LINE__);
        }
    }
    else
        INFO("Continuation: corrector converged in " << newtonIter_ << " steps");

    INFO("Continuation: Newton corrector... done");
    return 0;
}

//======================================================================
template<typename Model, typename ParameterList>
int Continuation<Model, ParameterList>::
runBackTracking(VectorPtr stateDir, double parDir)
{
    // Initialize reduction with -1/2
    double reduction = -1.0 / 2;
    double increase  = backTrackIncrease_;

    for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_)
    {
        // we allow for a some increase
        if (normRHStest_ < normRHS_ * increase)
            break;

        // Apply reduction to the model state and its parameter
        stateView_->Update(reduction, *stateDir, 1.0); // update the state
        par_ = par_ + reduction * parDir;              // update our parameter
        model_->setPar(par_);                          // set it in the model

        // Compute new RHS and obtain its norm
        model_->computeRHS();
        normRHStest_ = Utils::norm(rhsView_);

        INFO("    backtracking step " << backTrack_ <<
             ", norm: " << normRHStest_ << " < " << normRHS_ * increase << " ? ");

        reduction /= 2.0;   // Update reduction
    }

    TRACK_ITERATIONS("Continuation: backtracking steps...", backTrack_);

    if ( (normRHStest_ > normRHS_ * increase) &&
         (numBackTrackingSteps_ > 0) )
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
    size_t idx = 0;
    bool busy = false;
    for (auto &dest: destinations_)
    {
        idx++; // destination counter

        par_ = model_->getPar(); // just to be on the safe side
        double f0 = storage_.par0 - dest;
        double f1 =          par_ - dest;

        if (signMonitor_ == 0) // initialization
            signMonitor_ = SGN(f1);

        if (signMonitor_ != SGN(f1) && !secant_) // secant convergence
        {
            INFO("Continuation::detect():  sign switch in f1 detected");
            secant_     = true;
            activeDest_ = dest;
        }

        if (std::abs(f1) < 10 * destinationTolerance_ && !secant_)
        {
            INFO("Continuation::detect():  f1 near destination");
            secant_     = true;
            activeDest_ = dest;
        }

        busy = (secant_ && (std::abs(activeDest_ - dest) < 1e-7));

        if (busy)
        {
            INFO("  secant: f1 = " << f1 << " f0 = " << f0);
            INFO("      old ds = " << ds_);
            ds_ = -f1 * ds_ / (f1 - f0);   // secant method
            INFO("      new ds = " << ds_);
            createTangent('S');
        }

        signMonitor_ = SGN(f1);

        if (busy && std::abs(f1) < destinationTolerance_)
        {
            if (idx == destinations_.size())
                destinationReached_ = true;
            else
            {
                // Get the algorithm ready to proceed with the continuation
                secant_      = false;
                signMonitor_ = 0;
                ds_          = dsInit_;
                destinations_.erase(destinations_.begin() + idx - 1);
            }
        }

        if (busy)
        {
            INFO("  |par - dest| = " << std::abs(f1));
            break;
        }
    }
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
userDetect()
{
    if (userDetect_)
    {
        bool reached = model_->monitor();
        if (reached)
        {
            INFO("Continuation::userDetect() stopping criterion met");
            destinationReached_ = true;
        }
    }
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
adjustStep()
{
    // adjustStep is disabled during a secant process
    // or when requested with fixStepSize_
    if (secant_ || fixStepSize_)  
    {
        fixStepSize_ = false;
        return;   // Do nothing
    }

    if (newtonIter_ < minAllowedIterations_)
    {
        INFO("\n      adjusting ds... ");
        if (std::abs(ds_) == std::abs(dsMax_))
        {
            INFO("      ceiling found: " << ds_);
        }
        else
        {
            INFO("      old ds: " << ds_);
            ds_ = SGN(ds_) * std::min(std::abs(ds_) * scale1_, std::abs(dsMax_));
            INFO("      new ds: " << ds_);
        }
    }
    else if (newtonIter_ > maxAllowedIterations_)
    {
        INFO("      adjusting ds... ");
        if (std::abs(ds_) == std::abs(dsMin_))
        {
            INFO("      floor found: " << ds_);
        }
        else
        {
            INFO("      old ds: " << ds_);
            ds_ = SGN(ds_) * std::max(std::abs(ds_) / scale2_, std::abs(dsMin_));
            INFO("      new ds: " << ds_);
        }
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
    if (printImportantVectors_)
    {
        Utils::save(rhsView_, "failed_rhs"); // Print failed residual
        Utils::save(model_->getSolution('V'), "failed_sol"); // Print failed solution
        model_->dumpBlocks();
    }
    
    INFO("***********************************************");
    INFO("Continuation: reset... ");
    INFO("|                step: " << --step_); // decrease number of steps
    INFO("|           step size: " << ds_);
    INFO("|          norm state: " << Utils::norm(stateView_));
    INFO("|     norm d/ds state: " << Utils::norm(stateDot_));
    INFO("|          norm   rhs: " << Utils::norm(rhsView_));
    INFO("|  norm stored state0: " << Utils::norm(storage_.state0));
    INFO("|                 par: " << par_);
    INFO("| -------------------------------------------  ");

    restore(); // restore old model state

    // create new step size
    int s = SGN(ds_);
    ds_ = s * MAX(std::abs(ds_) / scale2_, std::abs(dsMin_));

    resetCounter_++;

    // fix step size for next step
    fixStepSize_ = true; 
    
    if (( std::abs(ds_) <= std::abs(dsMin_) ) && (( resetCounter_ >= 100 ) || giveUpAtdsMin_))
    {
        abortFlag_ = true;
        WARNING("Reached dsMin, continuation failed", __FILE__, __LINE__);
    }

    INFO("|      new norm state: " << Utils::norm(stateView_));
    INFO("| new norm d/ds state: " << Utils::norm(stateDot_));
    INFO("|             new par: " << par_);
    INFO("|       new step size: " << ds_);
    INFO("***********************************************");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
store()
{
    // We keep two previous states, state0 and state00
    storage_.state00   = storage_.state0;
    storage_.state0    = model_->getState('C');

    // We make a copy of the previous d/ds state
    storage_.stateDot0 = model_->getState('C');
    storage_.stateDot0->Update(1.0, *stateDot_, 0.0);

    // We keep two previous parameters and steps
    storage_.par00     = storage_.par0;
    storage_.par0      = model_->getPar();
    storage_.ds00      = storage_.ds0;
    storage_.ds0       = ds_;
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
restore()
{
    // Replace state in model with old state
    stateView_->Update(1.0, *storage_.state0, 0.0);

    model_->setPar(storage_.par0);

    par_      = storage_.par0;
    ds_       = storage_.ds0;
    stateDot_ = storage_.stateDot0;

    storage_.state0  = storage_.state00;
    storage_.state00 = model_->getState('C'); // why?
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
    
    VectorPtr view1 = model_->getState('V');
    VectorPtr copy1 = model_->getState('C');
    VectorPtr view2 = model_->getState('V');
    VectorPtr copy2 = model_->getState('C');
    
    INFO("view1      norm: " << Utils::norm(view1));
    INFO("view2      norm: " << Utils::norm(view2));
    INFO("copy1      norm: " << Utils::norm(copy1));
    INFO("copy2      norm: " << Utils::norm(copy2));

    INFO("Testing view");
    INFO("Scaling view1 with 2 ---------------");

    view1->Scale(2.0);
    
    INFO("view1      norm: " << Utils::norm(view1));
    INFO("view2      norm: " << Utils::norm(view2));

    assert(Utils::norm(view1) == Utils::norm(view2));

    INFO("Restore view: scale with 1/2");
    
    view1->Scale(0.5);
    
    INFO("Testing copy");
    INFO("Scaling copy1 with 3 ---------------");
    
    copy1->Scale(3.0);
    
    INFO("copy1      norm: " << Utils::norm(copy1));
    INFO("copy2      norm: " << Utils::norm(copy2));

    assert(Utils::norm(copy1) != Utils::norm(copy2));
    
    INFO("******************************************************");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
info()
{
    // get time of day
    std::time_t result = std::time(0); // hope this works on intel

    INFO(std::setprecision(8));
    INFO("-----------------------------------------");
    INFO("Continuation summary:  " << std::asctime(std::localtime(&result)));
    INFO(" destination reached:  " << destinationReached_);
    INFO("             aborted:  " << abortFlag_);
    INFO("           step size:  " << ds_);
    INFO("     parameter value:  " << par_);
    INFO("      starting value:  " << startingPar_);
    INFO("   destination value:  " << destinations_.back());
    INFO("           d/ds(par):  " << parDot_);

    if (maxSteps_ >= 0)
    {
        INFO("                step:  " << step_ << "/" << maxSteps_);
    }
    else
    {
        INFO("                step:  " << step_);
    }


    INFO("               ||x||:  " << Utils::norm(stateView_));
    INFO("         ||d/ds(x)||:  " << Utils::norm(stateDot_));
    INFO("        newton iters:  " << newtonIter_);
    INFO("    avg newton iters:  " << (double) sumNewtonIter_ / (step_ + resetCounter_));
    INFO("              resets:  " << resetCounter_);
    INFO("-----------------------------------------\n");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
modelInfo()
{
    INFO(std::setprecision(10));
    INFO("-----------------------------------------");
    INFO("Model view: ");
    INFO("     norm state :  " << Utils::norm(stateView_));
    INFO("     norm rhs   :  " << Utils::norm(rhsView_));
    INFO("     norm sol   :  " << Utils::norm(solView_));
    INFO("    parameter   :  " << model_->getPar());
    INFO("-----------------------------------------");
}

//======================================================================
template<typename Model, typename ParameterList>
void Continuation<Model, ParameterList>::
writeData(bool describe)
{
    // Write continuation data
    std::ostringstream cdatastring;

    if (describe) // write description of entries
    {
        cdatastring << std::setw(_FIELDWIDTH_)
                    << "par" 
                    << std::setw(_FIELDWIDTH_) 
                    << "ds"
                    << std::setw(_FIELDWIDTH_) 
                    << "||x||"
                    << std::setw(_FIELDWIDTH_) 
                    << "||F||"
                    << std::setw(_FIELDWIDTH_/2) 
                    << "NR"
                    << model_->writeData(describe);

        WRITECDATA(cdatastring.str());
    }

    cdatastring.str("");
    cdatastring.clear();
    cdatastring.precision(_PRECISION_);

    cdatastring << std::scientific
                << std::setw(_FIELDWIDTH_) <<  par_
                << std::setw(_FIELDWIDTH_) <<  ds_
                << std::setw(_FIELDWIDTH_) <<  Utils::norm(stateView_)
                << std::setw(_FIELDWIDTH_) <<  Utils::norm(rhsView_)
                << std::setw(_FIELDWIDTH_/2) <<  newtonIter_

                << model_->writeData();    

    WRITECDATA(cdatastring.str());
}

#endif
