#include "TestDefinitions.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    RCP<Ocean> ocean;
    RCP<Epetra_Comm>  comm;
}

//------------------------------------------------------------------
TEST(ParameterListOcean, Initialization)
{
    bool failed = false;
    try
    {
        // Copy of the input parameters to validate against
        const Teuchos::ParameterList startParams(*Utils::obtainParams("ocean_params.xml", "Ocean parameters"));

        // Copy of the default parameters
        const Teuchos::ParameterList defaultParams(Ocean::getDefaultInitParameters());
        // The input parameters to validate with
        Teuchos::ParameterList oceanParams(*Utils::obtainParams("ocean_params.xml", "Ocean parameters"));

        ::testing::internal::CaptureStdout();
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
        ::testing::internal::GetCapturedStdout();

        // Parameters currently reported by Ocean
        const Teuchos::ParameterList& currentParams = ocean->getParameters();

        // Check that every parameter is used
        //
        // This MUST be before checkParameterListAgainstDefaultAndOverrides
        // because it marks everything as used...
        EXPECT_TRUE(checkParameters(oceanParams, checkUsedParameterEntry));

        // Check that every entry in oceanParams corresponds to the value in
        // startParams, missing entries are compared against defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(oceanParams, defaultParams, startParams));

        // Check that every parameter is used
        //
        // This MUST be before checkParameterListAgainstDefaultAndOverrides
        // because it marks everything as used...
        EXPECT_TRUE(checkParameters(currentParams, checkUsedParameterEntry));

        // Check that every entry reported by Ocean corresponds to the value in
        // startParams, missing entries are compared against defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(currentParams, defaultParams, startParams));
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
