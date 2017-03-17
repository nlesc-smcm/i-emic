#include "TestDefinitions.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    std::shared_ptr<Ocean>         ocean;
    std::shared_ptr<AtmospherePar> atmos;
    std::shared_ptr<CoupledModel>  coupledModel;
}

//------------------------------------------------------------------
TEST(Ocean, Initialization)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        RCP<Teuchos::ParameterList> oceanParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());

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
        RCP<Teuchos::ParameterList> atmosphereParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());

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
        RCP<Teuchos::ParameterList> coupledModelParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("coupledmodel_params.xml", coupledModelParams.ptr());
        coupledModel = std::make_shared<CoupledModel>(ocean,atmos,coupledModelParams);
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
    coupledModel->setPar(0.01);

    // try to converge
    int maxit = 3;
    std::shared_ptr<Combined_MultiVec> b;
    for (int i = 0; i != maxit; ++i)
    {
        coupledModel->computeRHS();

        coupledModel->computeJacobian();

        b = coupledModel->getRHS('C');

        CHECK_ZERO(b->Scale(-1.0));

        double normb = Utils::norm(b);
        double modelNorm = Utils::norm(coupledModel->getRHS('V'));

        coupledModel->solve(b);

        modelNorm = Utils::norm(coupledModel->getRHS('V'));

        INFO(" norm b              " << normb);
        INFO(" norm b in model     " << modelNorm);

        INFO(" ocean F  = " << Utils::norm(coupledModel->getRHS('V')->First()) );
        INFO(" atmos F  = " << Utils::norm(coupledModel->getRHS('V')->Second()) );

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
        INFO("               ||b||  = " << normb);
    }

    EXPECT_LT(Utils::norm(coupledModel->getRHS('V')), 0.01);
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
    // TESTING
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
