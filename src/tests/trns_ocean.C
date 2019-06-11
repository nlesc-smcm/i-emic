#include "TestDefinitions.H"
#include "ThetaStepper.H"
#include "Theta.H"

//------------------------------------------------------------------
namespace 
{
    RCP<Epetra_Comm>               comm;
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
        Teuchos::RCP<Theta<Ocean> > oceanTheta =
            Teuchos::rcp(new Theta<Ocean>(comm, params[OCEAN]));
        
        // Create ThetaStepper
        ThetaStepper<Teuchos::RCP<Theta<Ocean> >,
                     Teuchos::RCP<Teuchos::ParameterList> >
            stepper(oceanTheta, params[TIME]);

        // Run ThetaStepper
        int status = stepper.run();
        int sumK   = stepper.getSumK();
        EXPECT_EQ(status, 0);
        
        double normState = Utils::norm(oceanTheta->getState('V'));

        // Checking some values with original behavior
        EXPECT_NEAR(normState, 22.3251803, 1e-4);
        EXPECT_EQ(sumK, 32);

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
