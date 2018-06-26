#include "TestDefinitions.H"
#include "NumericalJacobian.H"

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
        std::vector<string> files = {"ocean_params.xml",
                                     "atmosphere_params.xml",
                                     "seaice_params.xml",
                                     "coupledmodel_params.xml",
                                     "continuation_params.xml"};

        std::vector<string> names = {"Ocean parameters",
                                     "Atmosphere parameters",
                                     "Sea ice parameters",
                                     "CoupledModel parameters",
                                     "Continuation parameters"};

        for (int i = 0; i != (int) files.size(); ++i)
            params.push_back(obtainParams(files[i], names[i]));

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
TEST(CoupledModel, Continuation)
{
    bool failed = false;
    try
    {
        // Create continuation
        Continuation<std::shared_ptr<CoupledModel>,
                     Teuchos::RCP<Teuchos::ParameterList> >
            continuation(coupledModel, params[CONT]);

        // Run continuation        
        continuation.run();
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Ocean, Integrate_E_min_P)
{
    Teuchos::RCP<Epetra_Vector> dA = atmos->getPIntCoeff('C');
    Teuchos::RCP<Epetra_Vector> E  = ocean->interfaceE();
    Teuchos::RCP<Epetra_Vector> P  = atmos->interfaceP();
    
    E->Update(-1.0, *P, 1.0);
    double I = Utils::dot(E, dA);
    EXPECT_LT(std::abs(I), 1e-7);

    Teuchos::RCP<Epetra_Vector> Msi = seaice->interfaceM();
    std::cout << *Msi << std::endl;
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
    return out;
}
