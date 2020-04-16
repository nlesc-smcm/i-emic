#include "TestDefinitions.H"

#include "Ocean.H"
#include "TransientFactory.H"

//------------------------------------------------------------------
namespace
{
    Teuchos::RCP<Epetra_Comm>               comm;
    std::shared_ptr<Ocean>         ocean;
    std::vector<Teuchos::RCP<Teuchos::ParameterList> > params;

    enum Ident { OCEAN, TIME};
}

//------------------------------------------------------------------
TEST(ParameterLists, Initialization)
{
    bool failed = false;
    try
    {
        std::vector<std::string> files = {"test_oceantransient.xml",
                                          "timestepper_params.xml"};

        std::vector<std::string> names = {"Ocean parameters",
                                          "Time stepper parameters"};

        for (int i = 0; i != (int) files.size(); ++i)
            params.push_back(Utils::obtainParams(files[i], names[i]));

        params[OCEAN]->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed,false);
}

//------------------------------------------------------------------
TEST(ThetaStepper, Run)
{
    bool failed = false;
    try
    {
        // Create Theta<Ocean> object
        Teuchos::RCP<Ocean> ocean =
            Teuchos::rcp(new Ocean(comm, params[OCEAN]));

        // Create ThetaStepper
        auto stepper = TransientFactory(ocean, params[TIME]);

        // Run ThetaStepper
        int status = stepper->run();
        int sumK   = stepper->total_newton_steps();
        EXPECT_EQ(status, 0);

        double normState = Utils::norm(ocean->getState('V'));

        // Checking some values with original behavior
        EXPECT_NEAR(normState, 37.03750142, 1e-4);
        EXPECT_EQ(sumK, 30);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
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
    ocean = std::shared_ptr<Ocean>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
