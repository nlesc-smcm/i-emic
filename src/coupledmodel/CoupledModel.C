
#include "CoupledModel.H"

#include <functional>

#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
// TODO: - factorize constructors for ocean+atmos, ocean+atmos+seaice
//         - create setup routine that builds the identifiers
//       - generalize this code

//==================================================================
// constructor
CoupledModel::CoupledModel(std::shared_ptr<Ocean> ocean,
                           std::shared_ptr<AtmospherePar> atmos,
                           std::shared_ptr<SeaIce> seaice,
                           Teuchos::RCP<Teuchos::ParameterList> params)
    :
    OCEAN(-1),
    ATMOS(-1),
    SEAICE(-1),
    ocean_(ocean),
    atmos_(atmos),
    seaice_(seaice),
    parName_          (params->get("Continuation parameter",
                                   "Combined Forcing")),
    
    solvingScheme_    (params->get("Solving scheme", 'C')),
    
    precScheme_       (params->get("Preconditioning", 'F')),

    useOcean_         (params->get("Use ocean",      true)),
    useAtmos_         (params->get("Use atmosphere", true)),
    useSeaIce_        (params->get("Use sea ice",    false)),
    
    syncCtr_          (0),
    solverInitialized_(false)
{
    // Check xml sanity
    if (!useOcean_ && !useAtmos_ && !useSeaIce_)
    {
        ERROR("At least one model should be active",
              __FILE__, __LINE__);
    }
    
    // set models and identifiers
    int ident = 0;

    if (useOcean_)
    {
        models_.push_back(ocean);
        OCEAN = ident++;
    }

    if (useAtmos_)
    {
        models_.push_back(atmos);
        ATMOS = ident++;
    }

    if (useSeaIce_)
    {
        models_.push_back(seaice);
        SEAICE = ident++;
    }

    // default construction
    stateView_ = std::make_shared<Combined_MultiVec>();
    solView_   = std::make_shared<Combined_MultiVec>();
    rhsView_   = std::make_shared<Combined_MultiVec>();

    for (auto &model: models_)
    {
        // create our collection of vector views
        stateView_->AppendVector(model->getState('V'));
        solView_->AppendVector(model->getSolution('V'));
        rhsView_->AppendVector(model->getRHS('V'));
        
        // notify the models of the current continuation parameter
        model->setParName(parName_);
    }

    // The landmask interface is still in the fortran code, so Ocean
    // is responsible. In the case we don't have an ocean there is
    // also no landmask. Communicate surface landmask:
    if (useOcean_)
    {
        LandMask mask = models_[OCEAN]->getLandMask();
        // Start at first model beyond Ocean
        for (size_t i = 1; i <  models_.size(); ++i)
            models_[i]->setLandMask(mask);
    }    


    // Initialize CouplingBlock matrix
    using Block = CouplingBlock<std::shared_ptr<Model>,
                                std::shared_ptr<Model> >;

    C_ = std::vector<std::vector<Block> >(models_.size(),
                                          std::vector<Block>(models_.size()));
    
    for (size_t i = 0; i != models_.size(); ++i)
        for (size_t j = 0; j != models_.size(); ++j)
        {
            if (i != j) // only off-diagonal blocks
            {
                C_[i][j] = Block(models_[i], models_[j]);
                INFO("Created CouplingBlock: " << C_[i][j].name());
            }
        }
   
    // Output parameters
    INFO("\nCoupledModel parameters:");
    INFO(*params);
    INFO("\n");

    // -->We could ask some model specific information here, but that
    // -->would need dynamic_pointer_casts...
    if (useOcean_)
    {
        auto oceanPtr = std::dynamic_pointer_cast<Ocean>(models_[OCEAN]);
        if (oceanPtr)
        {
            INFO("Ocean couplings: coupled_T = " << oceanPtr->getCoupledT() );
            INFO("                 coupled_S = " << oceanPtr->getCoupledS() );
            INFO("--------------------------------------\n");
        }
        else
        {
            ERROR("CoupledModel downcasting failed", __FILE__, __LINE__);
        }
    }

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

    for (size_t i = 0; i != models_.size(); ++i)
        for (size_t j = 0; j != models_.size(); ++j)
        {
            if (models_[i] != models_[j])
                models_[i]->synchronize(models_[j]);
        }
    
    TIMER_STOP("CoupledModel: synchronize...");
}

//------------------------------------------------------------------
void CoupledModel::computeJacobian()
{
    TIMER_START("CoupledModel: compute Jacobian");

    // Synchronize the states
    if (solvingScheme_ != 'D') { synchronize(); }

    for (size_t i = 0; i != models_.size(); ++i)
    {
        models_[i]->computeJacobian();  // Ocean
        if (solvingScheme_ == 'C')
        {
            for (size_t j = 0; j != models_.size(); ++j)
            {
                if (i != j)
                    C_[i][j].computeBlock();
            }
        }
    }

    TIMER_STOP("CoupledModel: compute Jacobian");
}

//------------------------------------------------------------------
void CoupledModel::computeRHS()
{
    TIMER_START("CoupledModel compute RHS");

    // Synchronize the states in the fully coupled case
    if (solvingScheme_ != 'D') { synchronize(); }

    for (auto &model: models_)
        model->computeRHS();

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
    bool testExpl   = solverParams->get("FGMRES explicit residual test",
                                        false);

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
    belosParamList->set("Explicit Residual Test", testExpl);
    belosParamList->set("Implicit Residual Scaling",
                        "Norm of Preconditioned Initial Residual");

    // Belos block FGMRES setup
    belosSolver_ =
        Teuchos::rcp(new Belos::BlockGmresSolMgr
                     <double, Combined_MultiVec, BelosOp<CoupledModel> >
                     (problem_, belosParamList) );

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

    for (auto &model: models_)
        model->buildPreconditioner();

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
    for (size_t i = 0; i != models_.size(); ++i)
        models_[i]->applyMatrix(*v(i), *out(i));

    if (solvingScheme_ == 'C')
    {
        // Obtain temporary vector
        Combined_MultiVec z(v);
        
        // Apply off-diagonal coupling blocks
        for (size_t j = 0; j != models_.size(); ++j)
        {
            // reset z
            z.PutScalar(0.0);

            // fill z
            for (size_t i = 0; i != models_.size(); ++i)
            {
                if (i != j)
                    C_[i][j].applyMatrix(*v(j), *z(i));
            }
            
            out.Update(1.0, z, 1.0);
        }            
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
    for (size_t i = 0; i != models_.size(); ++i)
        models_[i]->applyMassMat(*v(i), *out(i));
    
    TIMER_STOP("CoupledModel: apply mass matrix...");    
}

//------------------------------------------------------------------
void CoupledModel::applyPrecon(Combined_MultiVec const &x, Combined_MultiVec &z)
{
    TIMER_START("CoupledModel: apply preconditioner...");

    z.PutScalar(0.0);  // Initialize output

    if (precScheme_ == 'D' || solvingScheme_ != 'C')
    {
        atmos_->applyPrecon(*x(ATMOS), *z(ATMOS));
        ocean_->applyPrecon(*x(OCEAN), *z(OCEAN) );
    }
    else if ( (precScheme_ == 'B' || precScheme_ == 'C') && solvingScheme_ == 'C')
    {
        Combined_MultiVec tmp(x);
        tmp.PutScalar(0.0);

        atmos_->applyPrecon(*x(ATMOS), *z(ATMOS)); //  z2   = inv(M2)*x2
        C_[OCEAN][ATMOS].applyMatrix(*z(ATMOS), *tmp(OCEAN));   //  tmp1 = C12*x2
        tmp(OCEAN)->Update(1.0, *x(OCEAN), -1.0);    //  tmp1 = x1 - C12*x2
        ocean_->applyPrecon(*tmp(OCEAN), *z(OCEAN)); //  z1   = inv(M1)*tmp1

        if (precScheme_ == 'C')
        {
            C_[ATMOS][OCEAN].applyMatrix(*z(OCEAN), *tmp(ATMOS));     // tmp2 = C21*x1
            tmp(ATMOS)->Update(1.0, *x(ATMOS), 1.0);
            atmos_->applyPrecon(*tmp(ATMOS), *z(ATMOS));
        }
            
    }
    else if ( (precScheme_ == 'F' || precScheme_ == 'G') && solvingScheme_ == 'C')
    {
        /*
        //!--------------------------------------------------
        //!Forward Block Gauss-Seidel (0-based indexing):
        //!M_k z_k^{j+1} = x_k + \sum_{i=k+1}^{n-1} C_{ki} z_i^{j}
        //!                    - \sum_{ i=0 }^{k-1} C_{ki} z_i^{j+1}
        //!              =: b_k
        //!--------------------------------------------------
        */
         
        Combined_MultiVec tmp(x);          // create temporary array
        Combined_MultiVec b(x);            // initialize b_k with x_k

        double sign = 0.0;

        //--> this should be a parameter in xml and we should get rid
        //--> of 'G' and 'C'
        int maxPrecIters = (precScheme_ == 'G') ? 2 : 1;
            
        // for (int iters = 0; iters != maxPrecIters; ++iters)
        //     for (size_t k = 0; k != models_.size(); ++k)    // iterate over rows
        //     {
        //         for (size_t i = 0; i < models_.size(); ++i)  // create \sum_{i=0}^{k-1} C_{ki} z_i
        //         {
        //             if (i < k)
        //                 sign = -1.0;
        //             else if ( (i > k) && (iters > 0 ) ) // assuming z is initialised 0
        //                 sign =  1.0;
        //             else
        //                 continue; // skip i == k
                
        //             tmp(k)->PutScalar(0.0);                 // just to be sure
        //             C_[k][i].applyMatrix(*z(i), *tmp(k));   // MV
        //             b(k)->Update(1.0, *tmp(k), sign);       // add to b_k
        //         }
        //         models_[k]->applyPrecon(*b(k), *z(k));             // solve M_k  
        //     }
        

        tmp.PutScalar(0.0);
        ocean_->applyPrecon(*x(OCEAN), *z(OCEAN));     // z1   = inv(M1) * x1
        C_[ATMOS][OCEAN].applyMatrix(*z(OCEAN), *tmp(ATMOS));     // tmp2 = C21*z1
        tmp(ATMOS)->Update(1.0, *x(ATMOS), -1.0);    // tmp2 = x2 - C21*z1
        atmos_->applyPrecon(*tmp(ATMOS), *z(ATMOS)); // z2   = inv(M2)*tmp2 

        // if (precScheme_ == 'G')
        // {
        //     C_[OCEAN][ATMOS].applyMatrix(*z(ATMOS), *tmp(OCEAN));     // tmp1 = C12*z2
        //     tmp(OCEAN)->Update(1.0, *x(OCEAN), 1.0);       // tmp1 = x1 + C12*z2
        //     ocean_->applyPrecon(*tmp(OCEAN), *z(OCEAN));   //   z1 = inv(M1) * tmp1
        //     C_[ATMOS][OCEAN].applyMatrix(*z(OCEAN), *tmp(ATMOS));     // tmp2 = C21*z1
        //     tmp(ATMOS)->Update(1.0, *x(ATMOS), -1.0);    // tmp2 = x2 - C21*z1
        //     atmos_->applyPrecon(*tmp(ATMOS), *z(ATMOS)); // z2   = inv(M2)*tmp2 
        // }        
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

    // Parameter values are equal to the continuation parameter or 0.
    
    double parvalue = 0.0;
    parvalue = (std::abs(par_ocean) > 0.0) ? par_ocean : parvalue;
    parvalue = (std::abs(par_atmos) > 0.0) ? par_atmos : parvalue;

    return parvalue;
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
    DUMPMATLAB("C12", *(C_[OCEAN][ATMOS].getBlock()));
    DUMPMATLAB("C21", *(C_[ATMOS][OCEAN].getBlock()));
}
