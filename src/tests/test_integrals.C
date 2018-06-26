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

        // Allow dominant parameterlists. Not that this trick uses a
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

