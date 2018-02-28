#include "TestDefinitions.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    RCP<Epetra_Comm> comm;
    RCP<Ocean> ocean;
}

//------------------------------------------------------------------
TEST(Ocean, Continuation1)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        RCP<Teuchos::ParameterList> oceanParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("ocean_params_intt_init.xml", oceanParams.ptr());
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));

        // Create continuation params
        RCP<Teuchos::ParameterList> continuationParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params_intt_init.xml",
                                    continuationParams.ptr());

        // Create contination
        Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
            continuation(ocean, continuationParams);

        // Run continuation
        continuation.run();

        // Copy ocean_output.h5 to ocean_input.h5
        std::ifstream src("ocean_output.h5", std::ios::binary);
        std::ofstream dst("ocean_input.h5", std::ios::binary);
        dst << src.rdbuf();
        
        // Copy info
        src = std::ifstream("info_0.txt", std::ios::binary);
        dst = std::ofstream("info_0.old", std::ios::binary);
        dst << src.rdbuf();

        // Copy cdata
        src = std::ifstream("cdata.txt", std::ios::binary);
        dst = std::ofstream("cdata.old", std::ios::binary);
        dst << src.rdbuf();

        // Delete cdata        
        src.close();
        src.open("cdata.txt", std::ofstream::out | std::ofstream::trunc);
        src.close();
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Ocean, Continuation2)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        RCP<Teuchos::ParameterList> oceanParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("ocean_params_intt_cont.xml", oceanParams.ptr());
        ocean = Teuchos::null;
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
        
        // Create continuation params
        RCP<Teuchos::ParameterList> continuationParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params_intt_cont.xml",
                                    continuationParams.ptr());

        // Create contination
        Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
            continuation(ocean, continuationParams);

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
    ocean = Teuchos::null;
    
    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
