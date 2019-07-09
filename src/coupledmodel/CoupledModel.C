#include "CoupledModel.H"
#include "Ocean.H"
#include "Atmosphere.H"
#include "SeaIce.H"

#include <functional>

#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

//==================================================================
// constructor
CoupledModel::CoupledModel(std::shared_ptr<Model> ocean,
                           std::shared_ptr<Model> atmos,
                           std::shared_ptr<Model> seaice,
                           Teuchos::RCP<Teuchos::ParameterList> params)
    :
    OCEAN              (-1),
    ATMOS              (-1),
    SEAICE             (-1),
    syncCtr_           (0),
    solverInitialized_ (false)
{
    // set xml parameters
    setParameters(params);

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

    // common setup
    setup();
}

//------------------------------------------------------------------
// constructor
CoupledModel::CoupledModel(std::shared_ptr<Model> ocean,
                           std::shared_ptr<Model> atmos,
                           Teuchos::RCP<Teuchos::ParameterList> params)
    :
    OCEAN              (-1),
    ATMOS              (-1),
    SEAICE             (-1),
    syncCtr_           (0),
    solverInitialized_ (false)
{
    // set xml parameters
    setParameters(params);

    // in this setup seaice is not included
    useSeaIce_ = false;

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

    // common setup
    setup();
}

//------------------------------------------------------------------
void CoupledModel::setParameters(Teuchos::RCP<Teuchos::ParameterList> params)
{
    solvingScheme_ = params->get("Solving scheme",'C');
    precScheme_    = params->get("Preconditioning",'F');
    useOcean_      = params->get("Use ocean",true);
    useAtmos_      = params->get("Use atmosphere",true);
    useSeaIce_     = params->get("Use sea ice",false);
}

//------------------------------------------------------------------
void CoupledModel::setup()
{
    // Sanity check
    if ( (OCEAN < 0) &&
         (ATMOS < 0) &&
         (SEAICE < 0) )
    {
        ERROR("At least one model should be active",
              __FILE__, __LINE__);
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
    }

    // Create the GID2Coord mapping where we use the model ordering
    // that is in models_.
    createGID2CoordMap();

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

    // Synchronize state
    synchronize();
}

//------------------------------------------------------------------
void CoupledModel::createGID2CoordMap()
{
    int N, M, L, dof, aux, modelIdent = 0;

    gid2coord_.clear();

    for (auto &model: models_)
    {
        N = model->getDomain()->GlobalN();
        M = model->getDomain()->GlobalM();
        L = model->getDomain()->GlobalL();

        dof = model->getDomain()->Dof();

        // Auxiliary unknowns do not have a grid coordinate and are
        // appended at the end of an ordinary map.
        aux = model->getDomain()->Aux();

        for (int k = 0; k != L; ++k)
            for (int j = 0; j != M; ++j)
                for (int i = 0; i != N; ++i)
                    for (int xx = 0; xx < dof; ++xx)
                    {
                        gid2coord_.push_back({modelIdent, i, j, k, xx});
                    }

        if (aux > 0)
        {
            for (int a = 0; a < aux; ++a)
                gid2coord_.push_back({modelIdent, 0, 0, 0, dof + a});
        }

        modelIdent++;
    }

    assert((int) gid2coord_.size() == stateView_->GlobalLength());
}

//------------------------------------------------------------------
void CoupledModel::gid2coord(int const &gid, int &modelIdent,
                             int &i, int &j, int &k, int &xx)
{
    modelIdent = gid2coord_[gid][0];
    i          = gid2coord_[gid][1];
    j          = gid2coord_[gid][2];
    k          = gid2coord_[gid][3];
    xx         = gid2coord_[gid][4];
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
            if ( models_[i] != models_[j] )
                 models_[i]->synchronize<>(models_[j]);
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
    belosParamList->set("Verbosity",
                        Belos::Errors + Belos::Warnings);
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

    // initialize effort counter
    effortCtr_ = 0;
    effort_    = 0.0;

    INFO("CoupledModel: initialize FGMRES done");
}

//------------------------------------------------------------------
void CoupledModel::solve(std::shared_ptr<const Combined_MultiVec> rhs)
{
    // Start solve
    TIMER_START("CoupledModel: solve...");

    // FGMRES with the coupled system.
    // The type of coupling is determined in applyMatrix() and applyPrecon().
    FGMRESSolve(rhs);

    TIMER_STOP("CoupledModel: solve...");
}

//------------------------------------------------------------------
void CoupledModel::FGMRESSolve(std::shared_ptr<const Combined_MultiVec> rhs)
{
    INFO("CoupledModel: FGMRES solve");

    if (!solverInitialized_)
        initializeFGMRES();

    for (auto &model: models_)
        model->buildPreconditioner();

    Teuchos::RCP<Combined_MultiVec> solV =
        Teuchos::rcp(&(*solView_), false);

    Teuchos::RCP<const Combined_MultiVec> rhsV =
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

    // project checkerboard modes from solution
    // if (useOcean_)
    // {
    //     Teuchos::RCP<Epetra_Vector> solView = models_[OCEAN]->getSolution('V');
    //     models_[OCEAN]->pressureProjection(solView);
    // }

    int iters  = belosSolver_->getNumIters();
    bool loa   = belosSolver_->isLOADetected();
    if (loa)
        INFO(" CoupledModel: FGMRES loss of accuracy detected");

    double tol = belosSolver_->achievedTol();

    double normb = Utils::norm(rhs);
    double nrm = explicitResNorm(rhs);
    INFO("           ||b||         = " << normb);
    INFO("           ||x||         = " << Utils::norm(solView_));
    INFO("        ||b-Ax|| / ||b|| = " << nrm / normb);

    if ((tol > 0) && (normb > 0) && ( (nrm / normb / tol) > 10))
    {
        WARNING("Actual residual norm ten times larger: "
              << (nrm / normb) << " > " << tol
              , __FILE__, __LINE__);
    }

    // keep track of effort
    if (effortCtr_ == 0)
        effort_ = 0;

    effortCtr_++;
    effort_ = (effort_ * (effortCtr_ - 1) + iters ) / effortCtr_;

    INFO("CoupledModel: FGMRES, iters = " << iters << ", ||r|| = " << tol);
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

            // fill components of z
            for (size_t i = 0; i != models_.size(); ++i)
            {
                if (i != j)
                    C_[i][j].applyMatrix(*v(j), *z(i));
            }

            // update with complete z
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
//--> some code repetition can be resolved here
void CoupledModel::applyPrecon(Combined_MultiVec const &x, Combined_MultiVec &z)
{
    TIMER_START("CoupledModel: apply preconditioner...");

    z.PutScalar(0.0);  // Initialize output

    if (precScheme_ == 'D' || solvingScheme_ != 'C')
    {
        for (size_t i = 0; i != models_.size(); ++i)
            models_[i]->applyPrecon(*x(i), *z(i));
    }
    else if ( (precScheme_ == 'B' || precScheme_ == 'C') && solvingScheme_ == 'C')
    {
        /*
        //!--------------------------------------------------
        //!Backward Block Gauss-Seidel (0-based indexing):
        //!M_k z_k^{j+1} = x_k + \sum_{ i=0 }^{k-1} C_{ki} z_i^{j}
        //!                    - \sum_{i=k+1}^{n-1} C_{ki} z_i^{j+1}
        //!              =: b_k
        //!--------------------------------------------------
        */

        Combined_MultiVec tmp(x);  // create temporary array
        tmp.PutScalar(0.0);
        Combined_MultiVec b(x);    // create b

        //--> this should be a parameter in xml and we should get rid
        //--> of 'G' and 'C'
        int maxPrecIters = (precScheme_ == 'C') ? 2 : 1;

        double sign;
        for (int iters = 0; iters != maxPrecIters; ++iters)
            for (int k = (int) models_.size()-1; k >= 0; --k)
            {
                *b(k) = *x(k); // reinitialize b(k) with x(k)
                for (int i = 0; i < (int) models_.size(); ++i) // create summation
                {
                    if ( (i < k) && (iters > 0 ) ) // assuming z is initialised 0
                        sign = 1.0;
                    else if (i > k)
                        sign = -1.0;
                    else
                        continue;

                    tmp(k)->PutScalar(0.0);                // zero result
                    C_[k][i].applyMatrix(*z(i), *tmp(k));  // MV
                    b(k)->Update(sign, *tmp(k), 1.0);      // add MV to b
                }
                if ( (precScheme_ == 'C') &&
                     (iters == maxPrecIters-1) &&
                     (k == OCEAN) ) break; // skip final solve (ocean)

                models_[k]->applyPrecon(*b(k), *z(k));     // solve M
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
        tmp.PutScalar(0.0);
        Combined_MultiVec b(x);            // initialize b_k with x_k

        double sign = 0.0;

        //--> this should be a parameter in xml and we should get rid
        //--> of 'G' and 'C'
        int maxPrecIters = (precScheme_ == 'G') ? 2 : 1;


        for (int iters = 0; iters != maxPrecIters; ++iters)
            for (size_t k = 0; k != models_.size(); ++k)    // iterate over rows
            {
                *b(k) = *x(k);  // reinitialize b(k) with x(k)
                for (size_t i = 0; i < models_.size(); ++i)  // create summation
                {
                    if (i < k)
                        sign = -1.0;
                    else if ( (i > k) && (iters > 0 ) ) // assuming z is initialised 0
                        sign =  1.0;
                    else
                        continue;

                    tmp(k)->PutScalar(0.0);                 // just to be sure
                    C_[k][i].applyMatrix(*z(i), *tmp(k));   // MV
                    b(k)->Update(sign, *tmp(k), 1.0);       // add MV to b_k
                }
                models_[k]->applyPrecon(*b(k), *z(k));      // solve M_k
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
double CoupledModel::explicitResNorm(std::shared_ptr<const Combined_MultiVec> rhs)
{

    Combined_MultiVec b = *getSolution('C');
    b.PutScalar(0.0);

    applyMatrix(*solView_, b);          // A*x
    b.Update(1, *rhs, -1);              // b-A*x
    double resnorm = Utils::norm(&b);   // ||b-A*x||

    // Utils::save(b, "lsresidual");
    return resnorm;
}

//------------------------------------------------------------------
std::shared_ptr<Combined_MultiVec> CoupledModel::getSolution(char mode)
{
    // obtain solution based on mode
    if (mode == 'V') // View
        return solView_;
    else if (mode == 'C') // Copy
    {
        std::shared_ptr<Combined_MultiVec> out =
            std::make_shared<Combined_MultiVec>();

        for (auto &model: models_)
            out->AppendVector(model->getSolution('C'));

        return out;
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
        std::shared_ptr<Combined_MultiVec> out =
            std::make_shared<Combined_MultiVec>();

        for (auto &model: models_)
            out->AppendVector(model->getState('C'));

        return out;
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
        std::shared_ptr<Combined_MultiVec> out =
            std::make_shared<Combined_MultiVec>();

        for (auto &model: models_)
            out->AppendVector(model->getRHS('C'));

        return out;
    }
    else
    {
        WARNING("Invalid mode", __FILE__, __LINE__);
        return std::shared_ptr<Combined_MultiVec>();
    }
}

//------------------------------------------------------------------
double CoupledModel::getPar(std::string const &parName)
{
    // Parameter values are equal to the continuation parameter or 0.
    double par, out = 0.0;
    for (auto &model: models_)
    {
        par = model->getPar(parName);
        out = (std::abs(par) > 0.0) ? par : out;
    }

    return out;
}

//------------------------------------------------------------------
void CoupledModel::setPar(std::string const &parName, double value)
{
    for (auto &model: models_)
        model->setPar(parName, value);
}

//------------------------------------------------------------------
void CoupledModel::initializeState()
{
    for (auto &model: models_)
    {
        synchronize();
        model->initializeState();
    }
}

//------------------------------------------------------------------
void CoupledModel::preProcess()
{
    for (auto &model: models_)
        model->preProcess();
}

//------------------------------------------------------------------
void CoupledModel::postProcess()
{
    // If the solver is completely decoupled, this is the right
    // moment to synchronize
    if (solvingScheme_ == 'D')
        synchronize();

    // Let the models do their own post-processing
    for (auto &model: models_)
        model->postProcess();
}

//------------------------------------------------------------------
void CoupledModel::saveStateToFile(std::string const &filename)
{
    for (auto &model: models_)
    {
        std::stringstream outFile;
        outFile << model->name() << "_" << filename;
        model->saveStateToFile(outFile.str());
    }
}

//------------------------------------------------------------------
//! Gather important continuation data to use in summary file
std::string const CoupledModel::writeData(bool describe)
        {
            std::ostringstream datastring;
            if (describe)
            {
                datastring << std::setw(_FIELDWIDTH_/ 3)
                           << "MV";
                for (auto &model: models_)
                    datastring << model->writeData(describe) << " ";
            }
            else
            {
                datastring.precision(_PRECISION_);
                datastring << std::setw(_FIELDWIDTH_/3)
                           << std::round(effort_);
                effortCtr_ = 0;

                for (auto &model: models_)
                    datastring << model->writeData(describe) << " ";
            }

            return datastring.str();
        }

//------------------------------------------------------------------
void CoupledModel::dumpBlocks()
{
    std::stringstream ss;
    for (size_t i = 0; i != models_.size(); ++i)
    {
        ss.str("");
        ss << "J_" << models_[i]->name();
        DUMPMATLAB(ss.str().c_str(), *(models_[i]->getJacobian()));
        for (size_t j = 0; j != models_.size(); ++j)
        {
            ss.str("");
            if (i != j)
            {
                ss << "C_" << C_[i][j].name();
                DUMPMATLAB(ss.str().c_str(), *(C_[i][j].getBlock()));
            }
        }
    }
}
