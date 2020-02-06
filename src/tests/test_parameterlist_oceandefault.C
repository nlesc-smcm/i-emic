#include "TestDefinitions.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    RCP<Ocean> ocean;
    RCP<Epetra_Comm>  comm;
}

//------------------------------------------------------------------
TEST(ParameterListOcean, DefaultInitialization)
{
    bool failed = false;
    try
    {
        // Copy of the default parameters
        const Teuchos::ParameterList defaultParams = Ocean::getDefaultInitParameters();

        // Empty parameter configuration
        Teuchos::ParameterList oceanParams("Test List");

        ::testing::internal::CaptureStdout();
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
        ::testing::internal::GetCapturedStdout();

        // Copy of the configuration reported by Ocean
        const Teuchos::ParameterList& currentParams = ocean->getParameters();

        // Check that every parameter is used
        //
        // This MUST be before checkParameterListAgainstDefaultAndOverrides
        // because it marks everything as used...
        EXPECT_TRUE(checkParameters(oceanParams, checkUsedParameterEntry));

        // Check that every parameter in oceanParams matches defaultParams, and
        // check that they're all reported as defaulted
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(oceanParams, defaultParams));
        EXPECT_TRUE(checkParameters(oceanParams, checkDefaultParameterEntry));

        // Check that every parameter is used
        //
        // This MUST be before checkParameterListAgainstDefaultAndOverrides
        // because it marks everything as used...
        EXPECT_TRUE(checkParameters(currentParams, checkUsedParameterEntry));

        // Check that every parameter reported by Ocean matches defaultParams,
        // and check that they're all reported as defaulted
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(currentParams, defaultParams));
        EXPECT_TRUE(checkParameters(currentParams, checkDefaultParameterEntry));
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
    ocean = Teuchos::null;

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    if (comm->MyPID() == 0)
        printProfile();

    MPI_Finalize();
    return out;
}
