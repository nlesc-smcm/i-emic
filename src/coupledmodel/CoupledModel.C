#include "CoupledModel.H"
#include "Ocean.H"
#include "AtmospherePar.H"
#include "Combined_MultiVec.H"

#include <vector>
#include <memory>
#include <functional>

#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

//==================================================================
// constructor
CoupledModel::CoupledModel(std::shared_ptr<Ocean> ocean,
                           std::shared_ptr<AtmospherePar> atmos,
                           Teuchos::RCP<Teuchos::ParameterList> params)
    :
    ocean_(ocean),
    atmos_(atmos),

    stateView_(std::make_shared<Combined_MultiVec>
               (ocean->getState('V'), atmos->getState('V'))),
    solView_  (std::make_shared<Combined_MultiVec>
               (ocean->getSolution('V'), atmos->getSolution('V'))),
    rhsView_  (std::make_shared<Combined_MultiVec>
               (ocean->getRHS('V'), atmos->getRHS('V'))),

    parName_          (params->get("Continuation parameter",
                                   "Combined Forcing")),
    solvingScheme_    (params->get("Solving scheme", 'G')),

    iterGS_           (params->get("Max GS iterations", 10)),
    toleranceGS_      (params->get("GS tolerance", 1e-1)),
    syncCtr_          (0),
    solverInitialized_(false)
{
    // Let the sub-models know our continuation parameter
    ocean_->setParName(parName_);
    atmos_->setParName(parName_);

    // Communicate surface landmask
    LandMask mask = ocean_->getLandMask();
    atmos_->setLandMask(mask);

    C12_ = CouplingBlock<std::shared_ptr<Ocean>,
                         std::shared_ptr<AtmospherePar> >(ocean_, atmos_);

    C21_ = CouplingBlock<std::shared_ptr<AtmospherePar>,
                         std::shared_ptr<Ocean> >(atmos_, ocean_);

    // Output parameters
    INFO(*params);

    // Synchronize state
    synchronize();
}

//------------------------------------------------------------------
// Here we throw Epetra_Vectors around. Inside the models we should check
// that the maps are correct, which in fact checks that the
// domain decompositions are compatible.
void CoupledModel::synchronize()
{

    TIMER_START("CoupledModel: synchronize...");

    syncCtr_++; // Keep track of synchronizations

    // Set atmosphere data in the ocean
    // ocean_->synchronize(atmos_);

    // Set ocean data in atmosphere
    // atmos_->synchronize(ocean_);

    TIMER_STOP("CoupledModel: synchronize...");
}

//------------------------------------------------------------------
void CoupledModel::computeJacobian()
{
    TIMER_START("CoupledModel: compute Jacobian");

    // Synchronize the states
    if (solvingScheme_ != 'D') { synchronize(); }

    ocean_->computeJacobian();  // Ocean
    atmos_->computeJacobian();  // Atmosphere

    TIMER_STOP("CoupledModel: compute Jacobian");
}

//------------------------------------------------------------------
void CoupledModel::computeRHS()
{
    TIMER_START("CoupledModel compute RHS");

    // Synchronize the states in the fully coupled case
    if (solvingScheme_ != 'D') { synchronize(); }

    ocean_->computeRHS();   // Ocean
    atmos_->computeRHS();   // Atmosphere

#ifdef DEBUGGING_NEW
    INFO("CoupledModel::computeRHS ocean ||rhs|| = " << Utils::norm(ocean_->getRHS('V')));
    INFO("CoupledModel::computeRHS atmos ||rhs|| = " << Utils::norm(atmos_->getRHS('V')));

    Utils::print(ocean_->getRHS('V'), "oceanRHS");
    Utils::print(atmos_->getRHS('V'), "atmosRHS");
#endif

    TIMER_STOP("CoupledModel compute RHS");
}

//====================================================================
void CoupledModel::initializeFGMRES()
{
    INFO("CoupledModel: initialize FGMRES...");

    Teuchos::RCP<Teuchos::ParameterList> solverParams =
        rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("solver_params.xml", solverParams.ptr());

    // Construct matrix operator
    Teuchos::RCP<BelosOp<CoupledModel> > coupledMatrix =
        Teuchos::rcp(new BelosOp<CoupledModel>(*this, false) );

    // Construct preconditioning operator
    Teuchos::RCP<BelosOp<CoupledModel> > coupledPrec =
        Teuchos::rcp(new BelosOp<CoupledModel>(*this, true) );

    // Construct non-owning rcps to vectors
    Teuchos::RCP<Combined_MultiVec> solV =
        Teuchos::rcp(&(*solView_), false);

    Teuchos::RCP<Combined_MultiVec> rhsV =
        Teuchos::rcp(&(*rhsView_), false);

    // Construct linear problem
    problem_ =
        Teuchos::rcp(new Belos::LinearProblem
                     <double, Combined_MultiVec,
                     BelosOp<CoupledModel> >
                     (coupledMatrix, solV, rhsV));

    // Set preconditioning
    problem_->setRightPrec(coupledPrec);

    int gmresIters  = solverParams->get("FGMRES iterations", 400);
    double gmresTol = solverParams->get("FGMRES tolerance", 1e-2);
    int maxrestarts = solverParams->get("FGMRES restarts", 2);
    int output      = solverParams->get("FGMRES output", 20);

    int NumGlobalElements = stateView_->GlobalLength();
    int blocksize         = 1; // number of vectors in rhs
    int maxiters          = NumGlobalElements / blocksize - 1;

    // Create Belos parameterlist
    Teuchos::RCP<Teuchos::ParameterList> belosParamList =
        rcp(new Teuchos::ParameterList());

    belosParamList->set("Block Size", blocksize);
    belosParamList->set("Flexible Gmres", true);
    belosParamList->set("Adaptive Block Size", true);
    belosParamList->set("Num Blocks", gmresIters);
    belosParamList->set("Maximum Restarts", maxrestarts);
    belosParamList->set("Orthogonalization","DGKS");
    belosParamList->set("Output Frequency", output);
    belosParamList->set("Verbosity", Belos::TimingDetails +
                        Belos::Errors +
                        Belos::Warnings +
                        Belos::StatusTestDetails );
    belosParamList->set("Maximum Iterations", maxiters);
    belosParamList->set("Convergence Tolerance", gmresTol);

    // Belos block FGMRES setup
    belosSolver_ =
        Teuchos::rcp(new Belos::BlockGmresSolMgr
                     <double, Combined_MultiVec, BelosOp<CoupledModel> >
                     (problem_, belosParamList));

    solverInitialized_ = true;

    INFO("CoupledModel: initialize FGMRES done");
}

//------------------------------------------------------------------
void CoupledModel::solve(std::shared_ptr<Combined_MultiVec> rhs)
{
    // Start solve
    TIMER_START("CoupledModel: solve...");

    if (solvingScheme_ == 'D') // fully decoupled solve
    {
        ocean_->solve(rhs->First());
        atmos_->solve(rhs->Second());
    }
    else if (solvingScheme_ == 'B') // backward block GS solve
        blockGSSolve(rhs);
    else if (solvingScheme_ == 'F') // FGMRES (Belos) on complete matrix
        FGMRESSolve(rhs);
    else
        WARNING("(CoupledModel::Solve()) Invalid mode!",
                __FILE__, __LINE__);

    TIMER_STOP("CoupledModel: solve...");
}

//------------------------------------------------------------------
void CoupledModel::FGMRESSolve(std::shared_ptr<Combined_MultiVec> rhs)
{
    if (!solverInitialized_)
        initializeFGMRES();

    ocean_->buildPreconditioner();

    Teuchos::RCP<Combined_MultiVec> solV =
        Teuchos::rcp(&(*solView_), false);

    Teuchos::RCP<Combined_MultiVec> rhsV =
        Teuchos::rcp(&(*rhs), false);

    solV->PutScalar(0.0);

    bool set = problem_->setProblem(solV, rhsV);

    TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
                               "*** Belos::LinearProblem failed to setup");
    try
    {
        belosSolver_->solve();      // Solve
    }
    catch (std::exception const &e)
    {
        INFO("Ocean: exception caught: " << e.what());
    }

    int iters  = belosSolver_->getNumIters();
    bool loa   = belosSolver_->isLOADetected();
    if (loa)
        INFO(" CoupledModel: FGMRES loss of accuracy detected");

    double tol = belosSolver_->achievedTol();

    INFO(" CoupledModel: FGMRES, iters = " << iters << ", ||r|| = " << tol);

#ifdef DEBUGGING_NEW
    // compute residual
    std::shared_ptr<Combined_MultiVec> r = getSolution('C');
    applyMatrix(*solView_, *r);
    r->Update(1.0, *rhs, -1.0);
    double normb = Utils::norm(rhs);
    INFO(" CoupledModel: FGMRES ||x1|| = "
         << Utils::norm(solView_->First()));
    INFO(" CoupledModel: FGMRES ||x2|| = "
         << Utils::norm(solView_->Second()));
    INFO(" CoupledModel: FGMRES ||r1|| = "
         << Utils::norm(r->First()) / normb);
    INFO(" CoupledModel: FGMRES ||r2|| = "
         << Utils::norm(r->Second()) / normb);
#endif
}

//------------------------------------------------------------------
void CoupledModel::blockGSSolve(std::shared_ptr<Combined_MultiVec> rhs)
{
    // ***************************************************************
    // Notation: J = [A,B;C,D], x = [x1;x2], b = [b1;b2]
    //           M = [A, 0; 0, D], E = [0, 0; -C, 0], F = [0, -B; 0, 0]
    //
    // Symmetric block GS: (M-F)*x^{k+1/2} = E*x^{k} + b
    //                     (M-E)*x^{k+1)   = F*x^{k+1/2) + b
    //
    // This leads to iteratively solving   D*x2 = -C*x1 + b2
    //                                     A*x1 = -B*x2 + b1
    //
    // After the iteration we do a final solve with  D*x2 = -C*x1 + b2
    //  (because it's cheap)
    // ***************************************************************

    ERROR("blockGSSolve not implemented, look in git history", __FILE__, __LINE__);

}

//------------------------------------------------------------------
//      out = [J1 C12; C21 J2] * [v1; v2]
void CoupledModel::applyMatrix(Combined_MultiVec const &v,
                               Combined_MultiVec &out, char mode)
{
    TIMER_START("CoupledModel: apply matrix...");

    // Initialize output
    out.PutScalar(0.0);

    // Apply the diagonal blocks
    // *out.First()  = *v.First();
    ocean_->applyMatrix(*v.First(),  *out.First());

    //*out.Second() = *v.Second();
    atmos_->applyMatrix(*v.Second(), *out.Second());


    if (mode == 'C')
    {
        // Obtain temporary vector
        Combined_MultiVec z(v);
        z.PutScalar(0.0);

        // Apply coupling blocks
        // C12_.applyMatrix(*v.Second(), *z.First());
        // C21_.applyMatrix(*v.First(),  *z.Second());

        out.Update(1.0, z, 1.0);
    }
    TIMER_STOP("CoupledModel: apply matrix...");

    #ifdef DEBUGGING_NEW
    #endif
}

//------------------------------------------------------------------
void CoupledModel::applyPrecon(Combined_MultiVec const &v,
                               Combined_MultiVec &out, char mode)
{
    TIMER_START("CoupledModel: apply preconditioner2...");

    out.PutScalar(0.0);     // Initialize output

    ocean_->applyPrecon(*v.First(),  *out.First() );

    atmos_->applyPrecon(*v.Second(), *out.Second());

    // *out.First()  = *v.First();
    // *out.Second() = *v.Second();

    TIMER_STOP("CoupledModel: apply preconditioner2...");
}

//------------------------------------------------------------------
double CoupledModel::computeResidual(std::shared_ptr<Combined_MultiVec> rhs)
{

    Combined_MultiVec b = *getSolution('C');
    b.PutScalar(0.0);

    applyMatrix(*solView_, b);

    double rhsNorm = Utils::norm(rhs);

    b.Update(1, *rhs, -1); //  b-Jx
    b.Scale(1.0 / rhsNorm);
    double relResidual = Utils::norm(&b); // ||b-Jx||/||b||

    return relResidual;
}

//------------------------------------------------------------------
std::shared_ptr<Combined_MultiVec> CoupledModel::getSolution(char mode)
{
    // obtain solution based on mode
    if (mode == 'V') // View
        return solView_;
    else if (mode == 'C') // Copy
    {
        return std::make_shared<Combined_MultiVec>(
            ocean_->getSolution('C'),
            atmos_->getSolution('C'));
    }
    else
    {
        WARNING("Invalid mode", __FILE__, __LINE__);
        return std::shared_ptr<Combined_MultiVec>();
    }
}

//------------------------------------------------------------------
std::shared_ptr<Combined_MultiVec> CoupledModel::getState(char mode)
{
    if (mode == 'V') // View
        return stateView_;
    else if (mode == 'C') // Copy
    {
        return std::make_shared<Combined_MultiVec>(
            ocean_->getState('C'),
            atmos_->getState('C'));
    }
    else
    {
        WARNING("Invalid mode", __FILE__, __LINE__);
        return std::shared_ptr<Combined_MultiVec>();
    }
}

//------------------------------------------------------------------
std::shared_ptr<Combined_MultiVec> CoupledModel::getRHS(char mode)
{
    if (mode == 'V') // View
        return rhsView_;
    else if (mode == 'C') // Copy
    {
        return std::make_shared<Combined_MultiVec>(
            ocean_->getRHS('C'),
            atmos_->getRHS('C'));
    }
    else
    {
        WARNING("Invalid mode", __FILE__, __LINE__);
        return std::shared_ptr<Combined_MultiVec>();
    }
}

//------------------------------------------------------------------
double CoupledModel::getPar()
{
    double par_ocean = ocean_->getPar(parName_);
    double par_atmos = atmos_->getPar(parName_);

    // In the case that the internal parameters are not the same,
    // we return the maximum. This happens when we perform continuations
    // in parameters that do not exist in all models.
    return std::max(par_ocean, par_atmos);
}

//------------------------------------------------------------------
void CoupledModel::setPar(double value)
{
    ocean_->setPar(parName_, value);
    atmos_->setPar(parName_, value);
}

//------------------------------------------------------------------
void CoupledModel::preProcess()
{
    ocean_->preProcess();
    atmos_->preProcess();
}

//------------------------------------------------------------------
void CoupledModel::postProcess()
{
    // If the solver is completely decoupled, this is the right
    // moment to synchronize
    if (solvingScheme_ == 'D')
        synchronize();

    // Let the models do their own post-processing
    ocean_->postProcess();
    atmos_->postProcess();
}
