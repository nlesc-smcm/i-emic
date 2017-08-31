#include "TestDefinitions.H"
#include <limits>
//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    std::shared_ptr<Ocean>         ocean;
    std::shared_ptr<AtmospherePar> atmos;
    std::shared_ptr<CoupledModel>  coupledModel;
    RCP<Teuchos::ParameterList> oceanParams;
    RCP<Teuchos::ParameterList> atmosphereParams;
    RCP<Teuchos::ParameterList> coupledmodelParams;
    RCP<Teuchos::ParameterList> continuationParams;
}

//------------------------------------------------------------------
TEST(ParameterLists, Initialization)
{
    bool failed = false;
    try
    {
        // Create parameter object for Ocean
        oceanParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
        oceanParams->setName("Ocean parameters");

        // Create parameter object for Atmosphere
        atmosphereParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());
        atmosphereParams->setName("Atmosphere parameters");

        // Create parameter object for CoupledModel
        coupledmodelParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("coupledmodel_params.xml", coupledmodelParams.ptr());
        coupledmodelParams->setName("CoupledModel parameters");

        // Create parameter object for Continuation
        continuationParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params.xml", continuationParams.ptr());
        continuationParams->setName("Continuation parameters");

        INFO('\n' << "Overwriting:");
        // The Continuation and CoupledModel parameterlists overwrite settings
        Utils::overwriteParameters(oceanParams,        coupledmodelParams);
        Utils::overwriteParameters(atmosphereParams,   coupledmodelParams);

        Utils::overwriteParameters(oceanParams,        continuationParams);
        Utils::overwriteParameters(atmosphereParams,   continuationParams);
        Utils::overwriteParameters(coupledmodelParams, continuationParams);
        INFO('\n');
    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed,false);
}

//------------------------------------------------------------------
TEST(Ocean, Initialization)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        ocean = std::make_shared<Ocean>(comm, oceanParams);
    }
    catch (...)
    {
        failed = true;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, Initialization)
{
    bool failed = false;
    try
    {
        // Create atmosphere
        atmos = std::make_shared<AtmospherePar>(comm, atmosphereParams);
    }
    catch (...)
    {
        failed = true;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Initialization)
{
    bool failed = false;
    try
    {
        // Create coupledmodel
        coupledModel = std::make_shared<CoupledModel>(ocean,atmos,coupledmodelParams);
    }
    catch (...)
    {
        failed = true;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, inspectState)
{
    bool failed = false;
    try
    {
        std::shared_ptr<Combined_MultiVec> state = coupledModel->getState('V');
        state->Random();

        int firstL  = state->First()->GlobalLength();
        int secondL = state->Second()->GlobalLength();
        int stateL  = state->GlobalLength();

        INFO(" global 1: " << firstL << " 2: " << secondL
             << " 1+2: " << stateL);

        EXPECT_EQ(firstL + secondL, stateL);

        firstL  = state->First()->MyLength();
        secondL = state->Second()->MyLength();
        stateL  = state->MyLength();

        INFO( " local 1: " << firstL << " 2: " << secondL
             << " 1+2: " << stateL );
        EXPECT_EQ(firstL + secondL, stateL);

        double firstNrm = Utils::norm(state->First());
        double secndNrm = Utils::norm(state->Second());
        double stateNrm = Utils::norm(state);

        EXPECT_NEAR(stateNrm, sqrt(pow(firstNrm,2) + pow(secndNrm,2)), 1e-7);

        INFO( " norm 1: " << firstNrm << " 2: " << secndNrm
              << " 1+2: " << stateNrm << std::endl );
    }
    catch (...)
    {
        failed = true;
    }

    EXPECT_EQ(failed, false);
}


//------------------------------------------------------------------
TEST(CoupledModel, computeJacobian)
{
    bool failed = false;
    try
    {
        coupledModel->computeJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> atmosJac = atmos->getJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> oceanJac = ocean->getJacobian();

        Utils::print(atmosJac, "atmosJac");
        Utils::print(oceanJac, "oceanJac");

        Utils::print(&atmosJac->ColMap(), "atmosJacColMap");
        Utils::print(&atmosJac->DomainMap(), "atmosJacDomainMap");

    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
// We need this information from THCM
extern "C" _SUBROUTINE_(getooa)(double*, double*);

TEST(CoupledModel, applyMatrix)
{
    bool failed = false;

    try
    {
        std::shared_ptr<Combined_MultiVec> x = coupledModel->getState('C');
        std::shared_ptr<Combined_MultiVec> y = coupledModel->getState('C');

        double normIn = Utils::norm(y);

        coupledModel->applyMatrix(*x, *y);

        double normOut = Utils::norm(y);

        EXPECT_NE(normIn, normOut);

        // Test the coupling blocks in the matrix separately

        CouplingBlock<std::shared_ptr<Ocean>,
                      std::shared_ptr<AtmospherePar> > C12(ocean, atmos);

        CouplingBlock<std::shared_ptr<AtmospherePar>,
                      std::shared_ptr<Ocean> > C21(atmos, ocean);

        Teuchos::RCP<Epetra_MultiVector> oceanVec = x->First();
        Teuchos::RCP<Epetra_MultiVector> atmosVec = x->Second();

        double value[3] = {1.234, 2.12, -23.5};

        // Test atmos -> ocean coupling
        for (int v = 0; v != 3; ++v)
        {
            atmosVec->PutScalar(value[v]);
            oceanVec->PutScalar(0.0);

            C12.applyMatrix(*atmosVec, *oceanVec);

            // Get ocean parameters
            double Ooa, Os;
            FNAME(getooa)(&Ooa, &Os );

            double max = 0;
            double el;
            for (int i = 0; i != oceanVec->MyLength(); ++i)
            {
                el = std::abs( (*oceanVec)[0][i] );
                max = (el > max) ? el : max;
            }

            EXPECT_NEAR( std::abs(-Ooa * value[v]), max, 1e-7);
        }

        // Test ocean -> atmos coupling
        for (int v = 0; v != 3; ++v)
        {
            oceanVec->PutScalar(0.0);
            for (int i = 4; i < oceanVec->MyLength(); i+=6)
            {
                (*oceanVec)[0][i] = value[v];
            }

            atmosVec->PutScalar(0.0);

            C21.applyMatrix(*oceanVec, *atmosVec);

            double max = 0;
            double el;
            for (int i = 0; i != atmosVec->MyLength(); ++i)
            {
                el = std::abs( (*atmosVec)[0][i] );
                max = (el > max) ? el : max;
            }
            EXPECT_NEAR( std::abs(value[v]), max, 1e-7);

            oceanVec->PutScalar(0.0);
            for (int i = 3; i < oceanVec->MyLength(); i+=6)
            {
                (*oceanVec)[0][i] = value[v];
            }

            atmosVec->PutScalar(0.0);
            C21.applyMatrix(*oceanVec, *atmosVec);

            max = 0;
            for (int i = 0; i != atmosVec->MyLength(); ++i)
            {
                el = std::abs( (*atmosVec)[0][i] );
                max = (el > max) ? el : max;
            }
            EXPECT_NE( std::abs(value[v]), max);

        }

    }
    catch (...)
    {
        failed = true;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, View)
{
    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    std::shared_ptr<Combined_MultiVec> rhsV =
        coupledModel->getRHS('V');

    stateV->PutScalar(0.0);
    coupledModel->computeRHS();

    double norm1 = Utils::norm(rhsV);

    stateV->PutScalar(1.0);
    coupledModel->computeRHS();

    double norm2 = Utils::norm(rhsV);

    EXPECT_NE(norm1, norm2);
}

//------------------------------------------------------------------
TEST(CoupledModel, Synchronization)
{
    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    stateV->First()->PutScalar(1.234);
    stateV->Second()->PutScalar(2.345);

    coupledModel->computeRHS();

    Teuchos::RCP<Epetra_Vector> oceanAtmosT = ocean->getLocalAtmosT();

    double maxValue, minValue, meanValue;
    oceanAtmosT->MaxValue(&maxValue);
    EXPECT_EQ(maxValue, 2.345);
    oceanAtmosT->MinValue(&minValue);
    EXPECT_EQ(minValue, 2.345);

    if (oceanAtmosT->Map().UniqueGIDs())
    {
        oceanAtmosT->MeanValue(&meanValue);
        EXPECT_NEAR(meanValue, 2.345, 1e-7);
    }
    else
        INFO(" oceanAtmosT (overl.) GID's not unique, which is fine in parallel.");

    Teuchos::RCP<Epetra_Vector> atmosOceanT = atmos->getLocalOceanT();

    Utils::print(atmosOceanT, "atmosOceanT" + std::to_string(comm->MyPID()));

    atmosOceanT->MaxValue(&maxValue);
    EXPECT_EQ(maxValue, 1.234);
    atmosOceanT->MinValue(&minValue);
    EXPECT_EQ(minValue, 1.234);

    if (atmosOceanT->Map().UniqueGIDs())
    {
        atmosOceanT->MeanValue(&meanValue);
        EXPECT_NEAR(meanValue, 1.234, 1e-7);
    }
    else
        INFO(" atmosOceanT (overl.) GID's not unique, which is fine in parallel.");
}

//------------------------------------------------------------------
TEST(CoupledModel, Newton)
{
    // One step in a 'natural continuation'
    // initialize state in model
    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    stateV->PutScalar(0.0);

    std::shared_ptr<Combined_MultiVec> solV =
        coupledModel->getSolution('V');

    solV->PutScalar(0.0);

    // set parameter
    coupledModel->setPar(0.005);

    // try to converge
    int maxit = 10;
    std::shared_ptr<Combined_MultiVec> b;
    int niter = 0;
    for (; niter != maxit; ++niter)
    {
        coupledModel->computeRHS();

        coupledModel->computeJacobian();

        b = coupledModel->getRHS('C');

        INFO(" ocean F  = " << Utils::norm(coupledModel->getRHS('V')->First()) );
        INFO(" atmos F  = " << Utils::norm(coupledModel->getRHS('V')->Second()) );

        CHECK_ZERO(b->Scale(-1.0));

        double normb = Utils::norm(b);

        coupledModel->solve(b);

        std::shared_ptr<Combined_MultiVec> x = coupledModel->getSolution('C');
        std::shared_ptr<Combined_MultiVec> y = coupledModel->getSolution('C');

        INFO(" ocean x  = " << Utils::norm(stateV->First()) );
        INFO(" atmos x  = " << Utils::norm(stateV->Second()) );
        INFO(" ocean dx = " << Utils::norm(x->First()) );
        INFO(" atmos dx = " << Utils::norm(x->Second()) );

        stateV->Update(1.0, *x, 1.0); // x = x + dx;

        coupledModel->applyMatrix(*x, *y);

        y->Update(1.0, *b, -1.0);
        y->Scale(1./normb);

        Utils::print(y, "residual");

        INFO(" ocean ||r|| / ||b||  = " << Utils::norm(y->First()));
        INFO(" atmos ||r|| / ||b||  = " << Utils::norm(y->Second()));
        INFO(" total ||r|| / ||b||  = " << Utils::norm(y));

        if (Utils::norm(coupledModel->getRHS('V')) < 0.01)
            break;
    }

    EXPECT_LT(Utils::norm(coupledModel->getRHS('V')), 0.01);
    EXPECT_LT(niter, 10);
    INFO("CoupledModel, Newton converged in " << niter << " iterations");
}

//------------------------------------------------------------------
// 1st integral condition test for atmosphere
TEST(CoupledModel, AtmosphereIntegralCondition1)
{
    Teuchos::RCP<Epetra_Vector> intCoeff = atmos->getIntCoeff();
    Teuchos::RCP<Epetra_Vector> atmosX   = atmos->getState('C');

    double result = Utils::dot(intCoeff, atmosX);

    INFO("  atmosphere state norm: " << Utils::norm(atmosX));
    INFO("  atmosphere integral condition on q: " << result);
    
    EXPECT_NEAR(result, 0.0, 1e-7);
}

//------------------------------------------------------------------
TEST(CoupledModel, AtmosphereEPfields)
{
    bool failed = false;
    try
    {
        Teuchos::RCP<Epetra_Vector> E = atmos->getE();
        Teuchos::RCP<Epetra_Vector> P = atmos->getP();
        std::ofstream fileE, fileP;
        fileE.open("fileE.txt");
        fileP.open("fileP.txt");
        E->Print(fileE);
        P->Print(fileP);
        fileE.close(); fileP.close();

        double tol = 1e-12;
    
        int nnzE = Utils::nnz(E, tol);
        int nnzP = Utils::nnz(P, tol);
        EXPECT_EQ(nnzE, nnzP);

        int nelE = E->Map().NumMyElements();
        int nelP = P->Map().NumMyElements();

        EXPECT_EQ(nelE, nelP);

        // expecting data distributions to be equal
        for (int i = 0; i != nelE; ++i)
        {
            if ( std::abs((*E)[i]) > tol )
            {
                EXPECT_GT( std::abs((*P)[i]), tol );
            }
        }
    }
    catch (...)
    {
        failed = true;
    }

    EXPECT_EQ(failed, false);    
}


//------------------------------------------------------------------
// full continuation
TEST(CoupledModel, Continuation)
{
    bool failed = false;
    try
    {
        // One step in an arclength continuation
        // initialize state in model
        std::shared_ptr<Combined_MultiVec> stateV =
            coupledModel->getState('V');
        
        stateV->PutScalar(0.0);
        
        std::shared_ptr<Combined_MultiVec> solV =
            coupledModel->getSolution('V');
        
        solV->PutScalar(0.0);
        
        // set initial parameter
        coupledModel->setPar(0.0);
        
        // Create continuation
        Continuation<std::shared_ptr<CoupledModel>,
                     Teuchos::RCP<Teuchos::ParameterList> >
            continuation(coupledModel, continuationParams);

        // Run continuation
        continuation.run();
    }
    catch (...)
    {
        failed = true;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
// 2nd integral condition test
TEST(CoupledModel, AtmosphereIntegralCondition2)
{
    Teuchos::RCP<Epetra_Vector> intCoeff = atmos->getIntCoeff();
    Teuchos::RCP<Epetra_Vector> atmosX   = atmos->getState('C');

    double result = Utils::dot(intCoeff, atmosX);

    INFO("  atmosphere state norm: " << Utils::norm(atmosX));
    INFO("  atmosphere integral condition on q: " << result);

    EXPECT_NEAR(result, 0.0, 1e-7);
}

//------------------------------------------------------------------
// using the solution from the previous continuation we check the
// integrity of the Jacobian matrix
TEST(CoupledModel, SmallPerturbation)
{
    std::shared_ptr<Combined_MultiVec> x  = coupledModel->getState('V');
    std::shared_ptr<Combined_MultiVec> xp = coupledModel->getState('C');
    xp->Scale(0.01);                // perturbation
    double nrmxp = Utils::norm(xp); 
    x->Update(1.0, *xp, 1.0);       // perturb state
    coupledModel->computeRHS();

    // temporary
    Combined_MultiVec tmp(*x);
    tmp.PutScalar(0.0);

    double nrmp = Utils::norm(coupledModel->getRHS('V')); // perturbed norm
    INFO("Perturbed norm: " << nrmp / nrmxp);

    x->Update(-1.0, *xp, 1.0);   // unperturb state
    coupledModel->computeRHS();
    coupledModel->computeJacobian();

    coupledModel->applyMatrix(*xp, tmp);
    tmp.Update(1.0,*coupledModel->getRHS('V'),1.0);

    double nrm = Utils::norm(&tmp); // perturbed norm
    INFO("Linearized norm: " << nrm / nrmxp);

    EXPECT_NEAR(nrm / nrmxp ,nrmp / nrmxp,1e-03);
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    initializeEnvironment(argc, argv);
    if (outFile == Teuchos::null)
        throw std::runtime_error("ERROR: Specify output streams");

    ::testing::InitGoogleTest(&argc, argv);

    // -------------------------------------------------------
    // TESTINGv
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    // Get rid of possibly parallel objects for a clean ending.
    ocean        = std::shared_ptr<Ocean>();
    atmos        = std::shared_ptr<AtmospherePar>();
    coupledModel = std::shared_ptr<CoupledModel>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;


    MPI_Finalize();
    return out;
}
