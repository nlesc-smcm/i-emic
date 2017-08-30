#include "TestDefinitions.H"

// Testing the serial and parallel atmosphere

//------------------------------------------------------------------
namespace
{
    std::shared_ptr<Atmosphere>    atmos;
    std::shared_ptr<AtmospherePar> atmosPar;
}

//------------------------------------------------------------------
TEST(Atmosphere, Initialization)
{
    bool failed = false;
    // Create atmosphere parameters
    RCP<Teuchos::ParameterList> atmosphereParams =
        rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());

    try
    {
        atmos = std::make_shared<Atmosphere>(atmosphereParams);
    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);

    failed = false;
    try
    {
        atmosPar = std::make_shared<AtmospherePar>(comm, atmosphereParams);
    }
    catch (std::exception const &e)
    {
        INFO("TEST(Atmosphere, Initialization) exception: " << e.what());
        failed = true;
    }
    catch (int error)
    {
        INFO("TEST(Atmosphere, Initialization) exception: error code = " << error);
        failed = true;
    }
    catch (...)
    {
        INFO("TEST(Atmosphere, Initialization) some exception thrown...");
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, RHS)
{
    bool failed = false;
    try
    {
        // Check dimensions
        int dim1 = atmos->dim();
        int dim2 = atmosPar->dim();
        EXPECT_EQ(dim1, dim2);

        // Create random state from atmosPar standardMap state
        Teuchos::RCP<Epetra_Vector> state_par = atmosPar->getState('C');

        state_par->SetSeed(rand() % 100);
        INFO("TEST(Atmosphere, RHS): seed = " << state_par->Seed());

        // Randomize and scale standardMap state
        state_par->Random();

        // Gather state_par into std::vector
        std::shared_ptr<std::vector<double> > state =
            getGatheredVector(state_par);

        // Put the full state in the serial atmosphere
        atmos->setState(state);

        // Put the standardMap parallel state into the parallel atmosphere
        atmosPar->setState(state_par);

        // Compute RHS in both models
        atmos->computeRHS();
        atmosPar->computeRHS();

        // Obtain RHS from both models
        Teuchos::RCP<Epetra_Vector> rhsPar = atmosPar->getRHS('V');
        std::shared_ptr<std::vector<double> > rhsSer = atmos->getRHS('V');

        // Check norms parallel and serial rhs
        EXPECT_NEAR(norm(rhsPar), norm(rhsSer), 1e-7);
        INFO("TEST(Atmosphere, RHS): serial norm = " << norm(rhsSer));
        INFO("TEST(Atmosphere, RHS): parall norm = " << norm(rhsPar));
    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, State)
{
    bool failed = false;
    try
    {
        // Obtain state from both models
        Teuchos::RCP<Epetra_Vector> statePar = atmosPar->getState('V');
        std::shared_ptr<std::vector<double> > stateSer = atmos->getState('V');

        // Check norms parallel and serial rhs
        EXPECT_NEAR(norm(statePar), norm(stateSer), 1e-7);
    }
    catch(...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, SurfaceTemperature)
{
    bool failed = false;
    try
    {
        // Obtain map for surface valeus
        Teuchos::RCP<Epetra_BlockMap> surfmap = atmosPar->getStandardSurfaceMap();

        // Create random vector from surface map
        Teuchos::RCP<Epetra_Vector> surfvals = Teuchos::rcp(new Epetra_Vector(*surfmap));
        surfvals->SetSeed(rand() % 100);
        INFO("TEST(Atmosphere, SurfaceTemperature): seed = " << surfvals->Seed());

        surfvals->Random();
        surfvals->Scale(1e1);

        // Create gather map
        Teuchos::RCP<Epetra_BlockMap> gmap =
            Utils::AllGather(*surfmap);

        // Create import strategy
        Teuchos::RCP<Epetra_Import> gimp =
            Teuchos::rcp(new Epetra_Import(*gmap, *surfmap));

        // Create gather vector
        Teuchos::RCP<Epetra_Vector> gathered = Teuchos::rcp(new Epetra_Vector(*gmap));

        // Gather <surfvals> into <gathered>
        gathered->Import(*surfvals, *gimp, Insert);

        // Create vector for serial atmosphere
        int size = surfmap->NumGlobalElements();
        std::shared_ptr<std::vector<double> > fullsurfvals =
            std::make_shared<std::vector<double> >(size);

        // Extract gathered
        gathered->ExtractCopy(&(*fullsurfvals)[0], size);

        // Put the randomized ocean surface values in both models
        atmos->setOceanTemperature(*fullsurfvals);
        atmosPar->setOceanTemperature(surfvals);

        // Compute RHS in both models
        atmos->computeRHS();
        atmosPar->computeRHS();

        // Obtain RHS from both models
        Teuchos::RCP<Epetra_Vector> rhsPar = atmosPar->getRHS('V');
        std::shared_ptr<std::vector<double> > rhsSer = atmos->getRHS('V');

        // Check norms parallel and serial rhs
        EXPECT_NEAR(norm(rhsPar), norm(rhsSer), 1e-6);
        INFO("TEST(Atmosphere, SurfaceTemperature): serial norm = " << norm(rhsSer));
        INFO("TEST(Atmosphere, SurfaceTemperature): parall norm = " << norm(rhsPar));

    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, Jacobian)
{
    bool failed = false;
    try
    {
        atmos->computeJacobian();
        atmosPar->computeJacobian();

        // Obtain vector and randomize
        Teuchos::RCP<Epetra_Vector> x = atmosPar->getState('C');
        x->Random();
        x->Scale(1e1);

        // Gather randomized vector
        std::shared_ptr<std::vector<double> > gx = getGatheredVector(x);

        // Check the gather routine
        EXPECT_NEAR(norm(gx), norm(x), 1e-7);

        std::shared_ptr<std::vector<double> > outSer =
            std::make_shared<std::vector<double> >(gx->size(), 0.0);

        // apply serial Jacobian
        atmos->applyMatrix(*gx, *outSer);

        Teuchos::RCP<Epetra_Vector> outPar = atmosPar->getState('C');
        outPar->PutScalar(0.0);

        // apply parallel Jacobian
        atmosPar->applyMatrix(*x, *outPar);

        // check results
        EXPECT_NEAR(norm(outSer), norm(outPar), 1e-7);

        INFO("TEST(Atmosphere, apply Jacobian): ||serial out|| = " << norm(outSer));
        INFO("TEST(Atmosphere, apply Jacobian): ||parall out|| = " << norm(outPar));
    }
    catch (std::exception const &e)
    {
        INFO("TEST(Atmosphere, Initialization) exception: " << e.what());
        failed = true;
    }
    catch (int error)
    {
        INFO("TEST(Atmosphere, Initialization) exception: error code = " << error);
        failed = true;
    }
    catch (...)
    {
        INFO("TEST(Atmosphere, Initialization) some exception thrown...");
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, EPfields)
{
    atmos->computeRHS();
    atmosPar->computeRHS();
    
    // get serial integral coefficients and evaporation field
    std::shared_ptr<std::vector<double> > serPco = atmos->getPrecipIntCo();
    std::shared_ptr<std::vector<double> > serE = atmos->getE();

    // get parallel integral coefficients and evaporation field
    Teuchos::RCP<Epetra_Vector> parPco = atmosPar->getPrecipIntCo();
    Teuchos::RCP<Epetra_Vector> parE = atmosPar->getE();

    // compute area
    double serArea, parArea;
    serArea = Utils::sum(*serPco);
    parPco->Norm1(&parArea);
    EXPECT_NEAR(serArea, parArea,1e-7);

    // compute dot products / integrals
    double serInt, parInt;
    serInt = Utils::dot(*serPco, *serE) / serArea;
    parInt = Utils::dot(parPco, parE) / parArea;
    EXPECT_NEAR(serInt, parInt, 1e-7);
    EXPECT_NE(serInt, 0.0);
    EXPECT_NE(parInt, 0.0);
    
}

//------------------------------------------------------------------
TEST(Atmosphere, Newton)
{
    Teuchos::RCP<Epetra_Vector> state = atmosPar->getState('V');
    state->PutScalar(0.0);
    atmosPar->setPar(0.5);

    Teuchos::RCP<Epetra_Vector> b = atmosPar->getRHS('V');
    Teuchos::RCP<Epetra_Vector> x = atmosPar->getSolution('V');

    int maxit = 10;
    int niter = 0;
    for (; niter != maxit; ++niter)
    {
        INFO(" old ||b|| = " << Utils::norm(b) );
        atmosPar->computeRHS();
        INFO(" new ||b|| = " << Utils::norm(b) );
        atmosPar->computeJacobian();

        Teuchos::RCP<Epetra_Vector> r = atmosPar->getSolution('C');

        b->Scale(-1.0);
        r->PutScalar(0.0);

        INFO(" old ||x|| = " << Utils::norm(x) );
        atmosPar->solve(b);
        INFO(" new ||x|| = " << Utils::norm(x) );

        atmosPar->applyMatrix(*x, *r);
        r->Update(1.0, *b, -1.0);

        EXPECT_NEAR(Utils::norm(r), 0, 1e-1);

        INFO("     ||r|| = " << Utils::norm(r) << std::endl);

        state->Update(1.0, *x, 1.0);
        b->Scale(-1.0);
        if (Utils::norm(b) < 1e-7)
            break;
    }

    INFO("Newton converged in " << niter << " iterations.");

    EXPECT_NEAR(Utils::norm(b), 0, 1e-7);
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
    atmos     = std::shared_ptr<Atmosphere>();
    atmosPar  = std::shared_ptr<AtmospherePar>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
