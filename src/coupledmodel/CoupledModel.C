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
    
    solvingScheme_    (params->get("Solving scheme", 'C')),
    
    precScheme_       (params->get("Preconditioning", 'F')),

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
    INFO("CoupledModel parameters:");
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

    // Set ocean data in atmosphere
    atmos_->synchronize(ocean_);

    // Set atmosphere data in the ocean
    ocean_->synchronize(atmos_);

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

    if (solvingScheme_ == 'C')
    {
        C12_.computeBlock();   // Recompute Ocean <- Atmos dependence
        C21_.computeBlock();   // Recompute Atmos <- Ocean dependence
    }

    TIMER_STOP("CoupledModel: compute Jacobian");
}

//------------------------------------------------------------------
void CoupledModel::computeRHS()
{
    TIMER_START("CoupledModel compute RHS");

    // Synchronize the states in the fully coupled case
    if (solvingScheme_ != 'D') { synchronize(); }

    ocean_->computeRHS();   // Ocean
    atmos_->computeRHS();   // Atmosphere/

#ifdef DEBUGGING_NEW
    // INFO("CoupledModel::computeRHS ocean ||rhs|| = " << Utils::norm(ocean_->getRHS('V')));
    // INFO("CoupledModel::computeRHS atmos ||rhs|| = " << Utils::norm(atmos_->getRHS('V')));

    // Utils::print(ocean_->getRHS('V'), "oceanRHS");
    // Utils::print(atmos_->getRHS('V'), "atmosRHS");
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

    int gmresIters  = solverParams->get("FGMRES iterations", 200);
    double gmresTol = solverParams->get("FGMRES tolerance", 1e-2);
    int maxrestarts = solverParams->get("FGMRES restarts", 0);
    int output      = solverParams->get("FGMRES output", 1000);

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

    // FGMRES with the coupled system.
    // The type of coupling is determined in applyMatrix() and applyPrecon().
    FGMRESSolve(rhs);

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
        INFO("CoupledModel: exception caught: " << e.what());
    }

    int iters  = belosSolver_->getNumIters();
    bool loa   = belosSolver_->isLOADetected();
    if (loa)
        INFO(" CoupledModel: FGMRES loss of accuracy detected");

    double tol = belosSolver_->achievedTol();

    INFO(" CoupledModel: FGMRES, iters = " << iters << ", ||r|| = " << tol);

}

//------------------------------------------------------------------
//      out = [J1 C12; C21 J2] * [v1; v2]
void CoupledModel::applyMatrix(Combined_MultiVec const &v, Combined_MultiVec &out)
{
    TIMER_START("CoupledModel: apply matrix...");

    // Initialize output
    out.PutScalar(0.0);

    // Apply the diagonal blocks
    ocean_->applyMatrix(*v.First(),  *out.First());
    atmos_->applyMatrix(*v.Second(), *out.Second());

    if (solvingScheme_ == 'C')
    {
        // Obtain temporary vector
        Combined_MultiVec z(v);
        z.PutScalar(0.0);

        // Apply off-diagonal coupling blocks
        C12_.applyMatrix(*v.Second(), *z.First());
        C21_.applyMatrix(*v.First(),  *z.Second());

        out.Update(1.0, z, 1.0);
    }
    TIMER_STOP("CoupledModel: apply matrix...");
}

//------------------------------------------------------------------
void CoupledModel::applyMassMat(Combined_MultiVec const &v, Combined_MultiVec &out)
{
    TIMER_START("CoupledModel: apply mass matrix...");

    // Initialize output
    out.PutScalar(0.0);

    // Apply mass matrix
    ocean_->applyMassMat(*v.First(), *out.First());
    atmos_->applyMassMat(*v.Second(), *out.Second());
    
    TIMER_STOP("CoupledModel: apply mass matrix...");    
}

//------------------------------------------------------------------
void CoupledModel::applyPrecon(Combined_MultiVec const &x, Combined_MultiVec &z)
{
    TIMER_START("CoupledModel: apply preconditioner...");

    z.PutScalar(0.0);  // Initialize output

    if (precScheme_ == 'D' || solvingScheme_ != 'C')
    {
        atmos_->applyPrecon(*x.Second(), *z.Second());
        ocean_->applyPrecon(*x.First() , *z.First() );
    }
    else if ( (precScheme_ == 'B' || precScheme_ == 'C') && solvingScheme_ == 'C')
    {
        Combined_MultiVec tmp(x);
        tmp.PutScalar(0.0);

        atmos_->applyPrecon(*x.Second(), *z.Second()); //  z2   = inv(M2)*x2
        C12_.applyMatrix(*z.Second(), *tmp.First());   //  tmp1 = C12*x2
        tmp.First()->Update(1.0, *x.First(), -1.0);    //  tmp1 = x1 - C12*x2
        ocean_->applyPrecon(*tmp.First(), *z.First()); //  z1   = inv(M1)*tmp1
        if (precScheme_ == 'C')
        {
            C21_.applyMatrix(*z.First(), *tmp.Second());     // tmp2 = C21*x1
            tmp.Second()->Update(1.0, *x.Second(), 1.0);
            atmos_->applyPrecon(*tmp.Second(), *z.Second());
        }
            
    }
    else if ( (precScheme_ == 'F' || precScheme_ == 'G') && solvingScheme_ == 'C')
    {
        Combined_MultiVec tmp(x);
        tmp.PutScalar(0.0);
        ocean_->applyPrecon(*x.First(), *z.First());     // z1   = inv(M1) * x1
        C21_.applyMatrix(*z.First(), *tmp.Second());     // tmp2 = C21*z1
        tmp.Second()->Update(1.0, *x.Second(), -1.0);    // tmp2 = x2 - C21*z1
        atmos_->applyPrecon(*tmp.Second(), *z.Second()); // z2   = inv(M2)*tmp2 

        if (precScheme_ == 'G')
        {
            C12_.applyMatrix(*z.Second(), *tmp.First());     // tmp1 = C12*z2
            tmp.First()->Update(1.0, *x.First(), 1.0);       // tmp1 = x1 + C12*z2
            ocean_->applyPrecon(*tmp.First(), *z.First());   //   z1 = inv(M1) * tmp1
            C21_.applyMatrix(*z.First(), *tmp.Second());     // tmp2 = C21*z1
            tmp.Second()->Update(1.0, *x.Second(), -1.0);    // tmp2 = x2 - C21*z1
            atmos_->applyPrecon(*tmp.Second(), *z.Second()); // z2   = inv(M2)*tmp2 
        }        
    }
    else
    {
        WARNING("Invalid prec/solve scheme: " << precScheme_ << " " << solvingScheme_,
                __FILE__, __LINE__);
    }

    TIMER_STOP("CoupledModel: apply preconditioner...");
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

//------------------------------------------------------------------
void CoupledModel::dumpBlocks()
{
    DUMPMATLAB("C11", *(ocean_->getJacobian()));
    DUMPMATLAB("C22", *(atmos_->getJacobian()));
    DUMPMATLAB("C12", *(C12_.getBlock()));
    DUMPMATLAB("C21", *(C21_.getBlock()));
}
