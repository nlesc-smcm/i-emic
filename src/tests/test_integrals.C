#include "TestDefinitions.H"

#include "NumericalJacobian.H"
#include "Ocean.H"
#include "Atmosphere.H"
#include "SeaIce.H"
#include "CoupledModel.H"
#include "Continuation.H"

namespace
{
    Teuchos::RCP<Epetra_Comm>             comm;
    std::shared_ptr<Ocean>               ocean;
    std::shared_ptr<Atmosphere>          atmos;
    std::shared_ptr<SeaIce>             seaice;
    std::shared_ptr<CoupledModel> coupledModel;
    std::vector<Teuchos::RCP<Teuchos::ParameterList> > params;
    enum Ident { OCEAN, ATMOS, SEAICE, COUPLED, CONT};
}

//------------------------------------------------------------------
TEST(ParameterLists, Initialization)
{
    bool failed = false;
    try
    {
        std::vector<std::string> files = {"ocean_params.xml",
                                          "atmosphere_params.xml",
                                          "seaice_params.xml",
                                          "coupledmodel_params.xml",
                                          "continuation_params.xml"};

        std::vector<std::string> names = {"Ocean parameters",
                                          "Atmosphere parameters",
                                          "Sea ice parameters",
                                          "CoupledModel parameters",
                                          "Continuation parameters"};

        for (int i = 0; i != (int) files.size(); ++i)
            params.push_back(Utils::obtainParams(files[i], names[i]));

        params[OCEAN]->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");

        INFO('\n' << "Overwriting:");

        // Allow dominant parameterlists. Note that this trick uses a
        // 'flattened' hierarchy. The Continuation and CoupledModel
        // parameterlists are allowed to overwrite settings.
        Utils::overwriteParameters(params[OCEAN],  params[COUPLED]);
        Utils::overwriteParameters(params[ATMOS],  params[COUPLED]);
        Utils::overwriteParameters(params[SEAICE], params[COUPLED]);

        Utils::overwriteParameters(params[OCEAN],  params[CONT]);
        Utils::overwriteParameters(params[ATMOS],  params[CONT]);
        Utils::overwriteParameters(params[SEAICE], params[CONT]);

        Utils::overwriteParameters(params[COUPLED], params[CONT]);
        INFO('\n');
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed,false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Initialization)
{
    bool failed = false;
    try
    {
        ocean  = std::make_shared<Ocean>(comm, params[OCEAN]);
        atmos  = std::make_shared<Atmosphere>(comm, params[ATMOS]);
        seaice = std::make_shared<SeaIce>(comm, params[SEAICE]);
        
        coupledModel =
            std::make_shared<CoupledModel>(ocean, atmos, seaice, params[COUPLED]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, RHS)
{
    double nrm;
    atmos->computeRHS();
    nrm = Utils::norm(atmos->getRHS());
    std::cout << "||F||atmos   = " << nrm << std::endl;
    EXPECT_LT(nrm, 1e-7);
    seaice->computeRHS();
    nrm = Utils::norm(seaice->getRHS());
    std::cout << "||F||seaice  = " << nrm << std::endl;
    EXPECT_LT(nrm, 1e-7);
    ocean->computeRHS();
    nrm = Utils::norm(ocean->getRHS());
    std::cout << "||F||ocean   = " << nrm << std::endl;
    EXPECT_LT(nrm, 1e-7);
    coupledModel->computeRHS();
    nrm = Utils::norm(coupledModel->getRHS());
    std::cout << "||F||coupled = " << nrm << std::endl;
    EXPECT_LT(nrm, 1e-7);
}

//------------------------------------------------------------------
TEST(CoupledModel, Continuation)
{
    bool failed = false;
    try
    {
        // Create continuation
        Continuation<std::shared_ptr<CoupledModel>>
            continuation(coupledModel, params[CONT]);

        // Run continuation        
        int status = continuation.run();
        EXPECT_EQ(status, 0);
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Integrate_E_min_P)
{
    Teuchos::RCP<Epetra_Vector> dA = atmos->getPIntCoeff('C');
    Teuchos::RCP<Epetra_Vector> E  = atmos->interfaceE();
    Teuchos::RCP<Epetra_Vector> P  = atmos->interfaceP();

    E->Update(-1.0, *P, 1.0);
    double I = Utils::dot(E, dA);
    EXPECT_LT(std::abs(I), 1e-7);

    Teuchos::RCP<Epetra_Vector> Msi = seaice->interfaceM();

    // restrict dA to sea ice
    Teuchos::RCP<Epetra_Vector> dAr = Teuchos::rcp(new Epetra_Vector(*dA));
    dAr->Multiply(1.0, *Msi, *dA, 0.0);

    // total sublimation
    E = atmos->interfaceE();
    double S = Utils::dot(E, dAr);
    std::cout << " total sublimation = " << S << std::endl;

    // FIXME: Should this always be greater than 0?
}

//------------------------------------------------------------------
TEST(CoupledModel, SeaIceCorrection)
{
    double oceanSCorr = ocean->getSCorr();
    Teuchos::RCP<Epetra_Vector> siCorr = seaice->interfaceG();

    std::cout << " Integral correction calculated by Ocean:  "
              << std::setprecision(12) << oceanSCorr << std::endl;
    std::cout << " Integral correction calculated by SeaIce: "
              << std::setprecision(12) << (*siCorr)[0] << std::endl;
    
    EXPECT_NEAR(oceanSCorr, (*siCorr)[0], 1e-10);
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    comm = initializeEnvironment(argc, argv);
    if (outFile == Teuchos::null)
        throw std::runtime_error("ERROR: Specify output streams");

    ::testing::InitGoogleTest(&argc, argv);

    // -------------------------------------------------------
    // TESTING
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    // Get rid of possibly parallel objects for a clean ending.
    ocean        = std::shared_ptr<Ocean>();
    atmos        = std::shared_ptr<Atmosphere>();
    seaice       = std::shared_ptr<SeaIce>();
    coupledModel = std::shared_ptr<CoupledModel>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();

    // print profile
    if (comm->MyPID() == 0)
        printProfile();

    return out;
}
