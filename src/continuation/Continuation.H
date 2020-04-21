//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef CONTINUATION_H
#define CONTINUATION_H

//======================================================================
#include "ContinuationDecl.H"
#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include <math.h> // pow(), sqrt()
#include <ctime>
#include <iomanip>

//======================================================================
//Constructor
template<typename Model>
Continuation<Model>::
Continuation(Model model)
    : Continuation<Model>(model, Teuchos::rcp(new Teuchos::ParameterList()))
{}

template<typename Model>
Continuation<Model>::
Continuation(Model model, Teuchos::RCP<Teuchos::ParameterList> pars)
    : Continuation<Model>(model, *pars)
{}

template<typename Model>
Continuation<Model>::
Continuation(Model model, Teuchos::ParameterList& pars)
    : model_(model)
    , paramList_("Continuation Configuration")
{

#ifdef HAVE_JDQZPP
    // Create JDQZ generalized eigenvalue solver:
    // 1) Create a vector with complex arithmetic based on an Epetra_Vector
    Vector t = *model->getSolution('C');
    t.PutScalar(0.0);
    ComplexVector<Vector> z(t);

    // 2) Create JDQZ (mass) matrix and preconditioning interface
    JDQZInterface<Model, ComplexVector<Vector> > interface(model, z);

    // 3) Create JDQZ solver
    jdqz_ = std::make_shared<JDQZsolver>(interface, z);

    // Create a JDQZ sublist so the parameter list validation does not fail.
    // This seems like a Trilinos bug.
    pars.sublist("JDQZ");
#endif

    pars.validateParametersAndSetDefaults(getDefaultInitParameters());
    paramList_.setParameters(pars);

#ifdef HAVE_JDQZPP
    jdqz_->setParameters(paramList_.sublist("JDQZ"));
#endif

    parName_               = paramList_.get<std::string>("continuation parameter");
    dsInit_                = paramList_.get<double>("initial step size");
    dsMin_                 = paramList_.get<double>("minimum step size");
    dsMax_                 = paramList_.get<double>("maximum step size");
    scale1_                = paramList_.get<double>("increase step size");
    scale2_                = paramList_.get<double>("decrease step size");
    epsilon_               = paramList_.get<double>("epsilon increment");
    backTracking_          = paramList_.get<bool>("enable backtracking");
    numBackTrackingSteps_  = paramList_.get<int>("backtracking steps");
    backTrackIncrease_     = paramList_.get<double>("backtracking increase");
    maxSteps_              = paramList_.get<int>("maximum number of steps");
    maxNewtonIterations_   = paramList_.get<int>("maximum Newton iterations");
    minNewtonIterations_   = paramList_.get<int>("minimum Newton iterations");
    optNewtonIterations_   = paramList_.get<double>("optimal Newton iterations");
    newtonTolerance_       = paramList_.get<double>("Newton tolerance");
    destinationTolerance_  = paramList_.get<double>("destination tolerance");
    userDetect_            = paramList_.get<bool>("enable custom monitor");
    detectMode_            = paramList_.get<char>("detection of special points");
    tanScaling_            = paramList_.get<double>("state tangent scaling");
    normalizeStrategy_     = paramList_.get<char>("normalize strategy");
    eigenvalueAnalysis_    = paramList_.get<char>("eigenvalue analysis");
    rejectFailedNewton_    = paramList_.get<bool>("reject failed iteration");
    giveUpAtdsMin_         = paramList_.get<bool>("give up at minimum step size");
    newtChordHybr_         = paramList_.get<bool>("enable Newton Chord hybrid solve");
    tangentType_           = paramList_.get<char>("tangent type");
    residualTest_          = paramList_.get<char>("corrector residual test");
    initialTangent_        = paramList_.get<char>("initial tangent type");
    printImportantVectors_ = paramList_.get<bool>("print important vectors");
    postProcess_           = paramList_.get<std::string>("post processing");
    predictorBound_        = paramList_.get<double>("predictor bound");

    // Set the step size
    ds_      = dsInit_;
    dsStart_ = dsInit_;

    // Get the parameter destinations
    // We allow for convergence with destination tolerance at multiple destinations
    double dest;
    std::stringstream destID;
    for (int i = 0; i != maxNumDest_; ++i)
    {
        destID << "destination " << i;
        dest = paramList_.get<double>(destID.str().c_str());
        if (dest == -999.)
            break;

        destinations_.push_back(dest);
        destID.str("");
        destID.clear();
    }
    if (destinations_.empty())
        ERROR("No destinations given in continuation xml!",
              __FILE__, __LINE__);

    destinationsBackup_ = destinations_;

    INFO("Continuation: number of destinations = " << destinations_.size());

    pars = paramList_;
    INFO(std::endl << "Continuation parameter list: "
         << std::endl << pars << std::endl);

    if (postProcess_ != "at every point")
    {
        WARNING("Alternative post processing: " << postProcess_ << '\n' <<
                " This may affect synchronization", __FILE__, __LINE__);
    }
}

//======================================================================
//Destructor
template<typename Model>
Continuation<Model>::
~Continuation()
{
    INFO("Continuation destructor");
}

//======================================================================
template<typename Model>
void Continuation<Model>::
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
    par_         = model_->getPar(parName_);
    startingPar_ = model_->getPar(parName_);

    // Save the residual
    if (printImportantVectors_)
        Utils::save(rhsView_, "residual");

    // print some info on the model
    modelInfo();

    // initialize Storage struct
    storage_.ds0     = ds_;
    storage_.ds00    = ds_;
    storage_.par0    = par_;
    storage_.par00   = par_;
    storage_.parDot0 = 0.0;

    storage_.state0 = model_->getState('C');

    // initializations for detect()
    destinations_ = destinationsBackup_;
    signMonitor_  = std::vector<int>(destinations_.size(), 0);
    secant_       = false;

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
    reachedLastDest_    = false;
    abortFlag_          = false;
    fixStepSize_        = false;

    if (userDetect_) { INFO("Continuation: custom monitor enabled");}
    else { INFO("Continuation: custom monitor disabled");}
}

//======================================================================
template<typename Model>
int Continuation<Model>::
run()
{
    INFO("Continuation::run initialize...");
    initialize();    // core initializations
    INFO("Continuation::run initialize... done");

    TIMER_START("Continuation: run");

    createInitialTangent(); // Create the first tangent

    int status = 0;

    while (!reachedLastDest_ && (step_ != maxSteps_) && !abortFlag_)
    {
        ++step_;          // Increment step counter
        info();           // Print continuation summary
        store();          // Store the current state and parameter
        status = step();  // Perform a continuation step

        if (status) // If step failed we reset and continue
        {
            reset();      // Reset the model
            continue;     // Skip the rest
        }

        detect();     // Detect special points
        userDetect(); // Use additional targets provided by model class
        adjustStep(); // step size adjustment
    }
    TIMER_STOP("Continuation: run");

    if (abortFlag_)
    {
        WARNING("Continuation aborted!",__FILE__, __LINE__);
        return 1;
    }

    finalize();

    INFO("---------Finished continuation run--------------");
    return 0;
}

//======================================================================
template<typename Model>
int Continuation<Model>::
step()
{
    TIMER_START("Continuation: step");
    // Let the model do some administrative work at the beginning of a step
    model_->preProcess();

    int status = 0;
    status = eulerPredictor();  // Apply Euler predictor

    // If necessary reset the step, otherwise perform a normal
    // calculation of the tangent and step adjustment
    if (status)   // Failure
    {
        TIMER_STOP("Continuation: step");
        return 1; // Consider this step failed
    }

    TIMER_START("Continuation: Newton");
    status = newtonCorrector(); // Apply Newton corrector
    TIMER_STOP("Continuation: Newton");

    if (status)   // Failure
    {
        TIMER_STOP("Continuation: step");
        return 1; // Consider this step failed
    }

    if (printImportantVectors_)
    {
        TIMER_START("Continuation: step -> IO");
        // Save the residual and last solution from linear solve to binaries
        Utils::save(rhsView_, "residual");
        Utils::save(model_->getSolution('V'), "solution");
        TIMER_STOP("Continuation: step -> IO");
    }

    // Put the parameter and norm of the state in the history
    parHist_.push_back(par_);
    stateNormHist_.push_back(Utils::norm(stateView_));

    // Inspect the history for weird behaviour
    analyzeHist();

    // Create new tangents based on result from newtonCorrector
    createTangent(tangentType_);

    // stability scheme = P => eigenvalues computed at every converged step
    if (eigenvalueAnalysis_ == 'P')
        eigenSolver();

    // Let the model do some administrative work at the end of a succesful step
    if (postProcess_ == "at every point")
    {
        TIMER_START("Continuation: step -> postprocess");
        model_->postProcess();
        TIMER_STOP("Continuation: step -> postprocess");
    }

    // Write continuation data to file
    bool describe = (step_ == 1) ? true : false;
    writeData(describe);

    TIMER_STOP("Continuation: step");
    return 0; // Exiting normally
}


//======================================================================
template<typename Model>
void Continuation<Model>::
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
    if (printImportantVectors_)
    {
        Utils::save(rhsView_, "F0");
    }

    computeDFDPar('F');

    if (printImportantVectors_)
    {
        Utils::save(dFdPar_, "dFdpar");
        Utils::save(rhsView_, "F1");
    }

    // 2a) Solve J*dxdpar = -dFdPar_
    if ( (initialTangent_ == 'E') || (initialTangent_ == 'S'))
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
        Utils::save(stateDot_, fname.str());
    }

    // We scale the tangent of the state and parameter such they are
    // normalized.
    normalize();

    // Test the matrix and rhs
    VectorPtr tmp = model_->getSolution('C');
    VectorPtr one = model_->getSolution('C');
    tmp->Scale(0.0);
    one->PutScalar(1.0);

    model_->applyMatrix(*rhsView_, *tmp);

    INFO(std::setprecision(14));
    INFO("                   ||rhs|| = " << Utils::norm(rhsView_));
    INFO("                   ||J*F|| = " << Utils::norm(tmp));

    model_->applyMatrix(*one, *tmp);
    INFO("                   ||J*1|| = " << Utils::norm(tmp));

    INFO("                 ||state|| = " << Utils::norm(stateView_));
    INFO("              ||stateDot|| = " << Utils::norm(stateDot_));
    INFO("                  parDot   = " << parDot_);
    INFO("                      ds   = " << ds_);
    INFO(std::setprecision(8));

    if (printImportantVectors_)
    {
        Utils::save(tmp, "actionJ");
        Utils::save(model_->getSolution('V'), "stateDot");
    }
    INFO("Continuation: create initial tangent... done");
}

//======================================================================
template<typename Model>
void Continuation<Model>::
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
    model_->setPar(parName_, par_ + epsilon_);  // increment parameter --> par + eps
    model_->computeRHS();             // compute new RHS     --> F(par+eps)
    INFO("       |                            l = " << par_);
    INFO("       |                          eps = " << epsilon_);
    INFO("       |                  F(x, l+eps) = " << Utils::norm(rhsView_));
    model_->setPar(parName_, par_);             // restore parameter   --> par

    // Get a copy of the new RHS F(par+eps), put it in dFdPar_
    dFdPar_ = model_->getRHS('C');

    // Update dFdPar_ such that it contains (F(par+eps) - F(par))/eps
    dFdPar_->Update(-1.0 / epsilon_, *rhsCopy_, 1.0 / epsilon_);
    INFO("       |  (F(x,l+eps) - F(x,l)) / eps = " << Utils::norm(dFdPar_));
}

//======================================================================
template<typename Model>
void Continuation<Model>::
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

        // Get the previous state, par and ds from our storage
        // These are used to compute a secant tangent
        double    par0   = storage_.par0;
        double    ds0    = storage_.ds0;

        // Compute stateDot = (state1 - state0)/ds0
        stateDot_->Update(-1.0 / ds0, *storage_.state0, 1.0 / ds0);

        // Compute parDot = (par1 - par0)/ds
        par_    = model_->getPar(parName_);      // update our par_
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

    INFO("Continuation::createTangent(), ||xDot|| = " << Utils::norm(stateDot_)
         << " parDot = " << parDot_);
}

//======================================================================
template<typename Model>
void Continuation<Model>::
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
template<typename Model>
int Continuation<Model>::
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
    model_->setPar(parName_, par_);

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
template<typename Model>
int Continuation<Model>::
newtonCorrector()
{
    INFO("Continuation: Newton corrector...");
    char      mode;          // mode for ComputeDFDPar()
    VectorPtr stateDir;      // direction for the state
    VectorPtr y;             // solution of solve with dFdPar
    VectorPtr z;             // solution of solve with F
    double    parDir = 0;    // direction for the parameter

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
        double rbp = 0;
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
        model_->setPar(parName_, par_);  // set it in the model

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
            INFO("                           ||R||_2 : " << res << " <? " << newtonTolerance_);
        }
        else if (residualTest_ == 'D')
        {
            INFO("                           ||R||_2 : " << normRHStest_);
            INFO("                     ||dx,dl||_inf : " << res << " <? " << newtonTolerance_);
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

        if (res < newtonTolerance_ && newtonIter_ >= minNewtonIterations_)
            break;
    }

    // Keep the last obtained y in stateDot, such that it can be used by
    // CreateTangent()
    if (!newtChordHybr_) stateDot_ = y;

    // Do some administrative work
    TRACK_ITERATIONS("Continuation: Newton iterations...", newtonIter_);

    if (res > newtonTolerance_)
    {
        INFO("Continuation: Newton failed after " << newtonIter_ << " steps");
        if (rejectFailedNewton_)
            return 1;
        else
        {
            INFO("Continuation: Proceeding with unconverged result");
        }
    }
    else
        INFO("Continuation: corrector converged in " << newtonIter_ << " steps");

    INFO("Continuation: Newton corrector... done");
    return 0;
}

//======================================================================
template<typename Model>
int Continuation<Model>::
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
        model_->setPar(parName_, par_);                          // set it in the model

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
template<typename Model>
void Continuation<Model>::detect()
{
    double dest = destinations_[0];

    par_ = model_->getPar(parName_); // just to be on the safe side

    double f0, f1;
    if (detectMode_ == 'D')
    {
        f0 = storage_.par0 - dest;
        f1 =          par_ - dest;
    }
    else if (detectMode_ == 'P')
    {
        f0 = storage_.parDot0;
        f1 = parDot_;
    }
    else
        ERROR("Invalid detectMode " << detectMode_,
              __FILE__, __LINE__);

    if (f1 == f0)
    {
        WARNING("This should not happen: f1 = " << f1 << " == "
                << f0 << " = f0" , __FILE__, __LINE__);
    }

    if (signMonitor_[0] == 0) // only at initialization
        signMonitor_[0] = SGN(f1);

    if (signMonitor_[0] != SGN(f1) && !secant_) // enable secant iteration
    {
        INFO("Continuation::detect(): sign switch in f1 detected");
        INFO("        activated dest: " << dest);
        secant_     =  true;   // set flag to start secant process
        dsStart_    =  ds_;    // save step size at beginning of secant
    }
    else
        signMonitor_[0] = SGN(f1);

    if (secant_)
    {
        INFO("    secant: f1 = " << f1 << " f0 = " << f0);
        INFO("        old ds = " << ds_);
        ds_ = -f1 * ds_ / (f1 - f0);   // secant method
        INFO("        new ds = " << ds_ << std::endl);
        createTangent('S');
    }

    // destination reached
    if (secant_ && std::abs(f1) < destinationTolerance_)
    {
        INFO("Continuation::detect(): destination " << dest << " reached. \n");

        // eigenvalue analysis = E => eigenvalues computed at every converged destination
        if (eigenvalueAnalysis_ == 'E')
            eigenSolver();

        // Get the algorithm ready to proceed with the continuation
        secant_ = false;     // disable secant method
        ds_     = dsStart_;  // reset step size to before secant

        // do not adjust step size for the next step
        fixStepSize_ = true;

        // remove reached destination and corresponding sign monitor
        destinations_.erase(destinations_.begin());
        signMonitor_.erase(signMonitor_.begin());

        if (destinations_.empty())
            reachedLastDest_  = true;
        else // initialize next sign monitor
            signMonitor_[0] = SGN(par_ - destinations_[0]);
    }
}

//======================================================================
template<typename Model>
void Continuation<Model>::
userDetect()
{
    if (userDetect_)
    {
        bool reached = model_->monitor();
        if (reached)
        {
            INFO("Continuation::userDetect() stopping criterion met");
            reachedLastDest_  = true;
        }
    }
}

//======================================================================
template<typename Model>
void Continuation<Model>::
adjustStep()
{
    // adjustStep is disabled during a secant process
    // or when requested with fixStepSize_
    if (secant_ || fixStepSize_)
    {
        fixStepSize_ = false;
        return;   // Do nothing
    }
    INFO("Continuation::adjustStep(), old ds: " << ds_);

    // step size control, see [Seydel p 188.]
    double factor = optNewtonIterations_ / (double) newtonIter_;

    // set some bounds for this factor
    factor = (factor < 0.5) ? 0.5 : factor;
    factor = (factor > 2.0) ? 2.0 : factor;

    INFO("                           control: " << factor);
    // adjust step size
    ds_ *= factor;

    // set some bounds for the step size
    ds_ = (std::abs(ds_) > std::abs(dsMax_)) ? SGN(ds_) * std::abs(dsMax_) : ds_;
    ds_ = (std::abs(ds_) < std::abs(dsMin_)) ? SGN(ds_) * std::abs(dsMin_) : ds_;

    INFO("                            new ds: " << ds_);

}

//======================================================================
template<typename Model>
void Continuation<Model>::
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
template<typename Model>
void Continuation<Model>::
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

    if (( std::abs(ds_) <= std::abs(dsMin_) ) &&
        (( resetCounter_ >= 100 ) || giveUpAtdsMin_))
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
template<typename Model>
void Continuation<Model>::
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
    storage_.par0      = model_->getPar(parName_);
    storage_.ds00      = storage_.ds0;
    storage_.ds0       = ds_;

    // We keep the previous parameter tangent
    storage_.parDot0   = parDot_;
}

//======================================================================
template<typename Model>
void Continuation<Model>::
restore()
{
    // Replace state in model with old state
    stateView_->Update(1.0, *storage_.state0, 0.0);

    model_->setPar(parName_, storage_.par0);

    par_      = storage_.par0;
    ds_       = storage_.ds0;
    stateDot_ = storage_.stateDot0;

    storage_.state0  = storage_.state00;
    storage_.state00 = model_->getState('C'); // why?
    storage_.par0    = storage_.par00;
    storage_.ds0     = storage_.ds00;
}

//=====================================================================
template<typename Model>
void Continuation<Model>::
eigenSolver()
{
    if (eigenvalueAnalysis_ != 'N')
    {
#ifdef HAVE_JDQZPP
        jdqz_->solve();

        // save eigenvectors
        std::stringstream ss;
        ss << "ev_step_" << step_;

        Utils::saveEigenvectors(jdqz_, ss.str());
#else
        WARNING("JDQZPP has not been installed!", __FILE__, __LINE__);
#endif
    }
    else
    {
        WARNING("Faulty stability scheme!", __FILE__, __LINE__);
    }
}

//=====================================================================
template<typename Model>
const Teuchos::ParameterList&
Continuation<Model>::
getParameters()
{ return paramList_; }

template<typename Model>
void Continuation<Model>::
setParameters(Teuchos::ParameterList& newParams)
{
    newParams.validateParameters(getDefaultParameters());
}

//======================================================================
template<typename Model>
void Continuation<Model>::
test()
{
    model_->computeRHS();
    testCopyView();
}

//======================================================================
template<typename Model>
void Continuation<Model>::
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

    copy1->Random();

    INFO("copy1      norm: " << Utils::norm(copy1));
    INFO("copy2      norm: " << Utils::norm(copy2));

    assert(Utils::norm(copy1) != Utils::norm(copy2));

    INFO("******************************************************");
}

//======================================================================
template<typename Model>
void Continuation<Model>::
info()
{
    // get time of day
    std::time_t result = std::time(0); // hope this works on intel

    INFO(std::setprecision(8));
    INFO("-----------------------------------------");
    INFO("Continuation summary:  " << std::asctime(std::localtime(&result)));
    INFO(" destination reached:  " << reachedLastDest_ );
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

    if ((step_ + resetCounter_) > 0)
    {
        INFO("    avg newton iters:  " << (double) sumNewtonIter_ / (step_ + resetCounter_));
    }

    INFO("              resets:  " << resetCounter_);
    INFO("-----------------------------------------\n");
}

//======================================================================
template<typename Model>
void Continuation<Model>::
finalize()
{
    info();
    modelInfo();

    if (printImportantVectors_)
        model_->dumpBlocks();

    if (postProcess_ == "at final point")
    {
        TIMER_START("Continuation: step -> postprocess");
        model_->postProcess();
        TIMER_STOP("Continuation: step -> postprocess");
    }
}

//======================================================================
template<typename Model>
void Continuation<Model>::
modelInfo()
{
    INFO(std::setprecision(10));
    INFO("-----------------------------------------");
    INFO("Model view: ");
    INFO("     norm state :  " << Utils::norm(stateView_));
    INFO("     norm rhs   :  " << Utils::norm(rhsView_));
    INFO("     norm sol   :  " << Utils::norm(solView_));
    INFO("    parameter   :  " << model_->getPar(parName_));
    INFO("-----------------------------------------");
}

//======================================================================
template<typename Model>
void Continuation<Model>::
writeData(bool describe)
{
    // Write continuation data
    std::ostringstream cdatastring;

    if (describe) // write description of entries
    {
        cdatastring << "#" << std::setw(_FIELDWIDTH_-1)
                    << "par"
                    << std::setw(_FIELDWIDTH_ * 3/4)
                    << "ds"
                    << std::setw(_FIELDWIDTH_ * 3/4)
                    << "||x||"
                    << std::setw(_FIELDWIDTH_ * 3/4)
                    << "||F||"
                    << std::setw(_FIELDWIDTH_/3)
                    << "NR"
                    << model_->writeData(describe);

        WRITECDATA(cdatastring.str());
    }

    cdatastring.str("");
    cdatastring.clear();

    cdatastring << std::scientific
                << std::setw(_FIELDWIDTH_) << std::setprecision(_PRECISION_)
                <<  par_
                << std::setw(_FIELDWIDTH_ * 3/4) << std::setprecision(_PRECISION_ / 2)
                <<  ds_
                << std::setw(_FIELDWIDTH_ * 3/4) << std::setprecision(_PRECISION_ / 2)
                << Utils::norm(stateView_)
                << std::setw(_FIELDWIDTH_ * 3/4) << std::setprecision(_PRECISION_ / 2)
                <<  Utils::norm(rhsView_)
                << std::setw(_FIELDWIDTH_ / 3) << std::setprecision(_PRECISION_ / 2)
                <<  newtonIter_

                << model_->writeData();

    WRITECDATA(cdatastring.str());
}

template<typename Model>
Teuchos::ParameterList
Continuation<Model>::
getDefaultInitParameters()
{
    Teuchos::ParameterList result = getDefaultParameters();
    result.setName("Default Init Continuation List");

#ifdef HAVE_JDQZPP
    JDQZsolver::getDefaultParameters(result.sublist("JDQZ"));
#endif

    result.get("continuation parameter", "Combined Forcing");
    result.get("initial step size", 1.0e-2);
    result.get("minimum step size", 1.0e-8);
    result.get("maximum step size", 1.0e3);
    result.get("increase step size", 1.25);
    result.get("decrease step size", 2.0);
    result.get("epsilon increment", 1.0e-5);
    result.get("enable backtracking", false);
    result.get("backtracking steps", 0);
    result.get("backtracking increase", 0.0);
    result.get("maximum number of steps", -1);
    result.get("maximum Newton iterations", 7);
    result.get("minimum Newton iterations", 1);
    result.get("optimal Newton iterations", 3.5);
    result.get("Newton tolerance", 1.0e-4);
    result.get("destination tolerance", 1.0e-7);
    result.get("enable custom monitor", false);
    result.get("detection of special points", 'D');
    result.get("state tangent scaling", 1.0e0);
    result.get("normalize strategy", 'N');
    result.get("eigenvalue analysis", 'N');
    result.get("reject failed iteration", true);
    result.get("give up at minimum step size", true);
    result.get("enable Newton Chord hybrid solve", false);
    result.get("tangent type", 'S');
    result.get("corrector residual test", 'D');
    result.get("initial tangent type", 'E');
    result.get("print important vectors", false);

    Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int>> post_processing_validator(
        new Teuchos::StringToIntegralParameterEntryValidator<int>(
            Teuchos::tuple<std::string> ("at every point", "at final point"), "post processing"));

    result.set("post processing", "at every point",
               "When to call post processing as implemented by the Model",
               post_processing_validator);

    result.get("predictor bound", 1e3);

    std::stringstream destID;
    for (int i = 0; i != maxNumDest_; ++i)
    {
        destID << "destination " << i;
        result.set(destID.str().c_str(), -999.0);
        destID.str("");
        destID.clear();
    }

    return result;
}

template<typename Model>
Teuchos::ParameterList
Continuation<Model>::
getDefaultParameters()
{
    Teuchos::ParameterList result("Default Continuation List");
    return result;
}
#endif
