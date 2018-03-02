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
        std::ifstream src1("ocean_output.h5", std::ios::binary);
        std::ofstream dst1("ocean_input.h5", std::ios::binary);
        dst1 << src1.rdbuf();
        
        // Copy info
        std::ifstream src2("info_0.txt", std::ios::binary);
        std::ifstream dst2("info_0.old", std::ios::binary);
        dst2 << src2.rdbuf();

        // Copy cdata
        std::ifstream src3("cdata.txt", std::ios::binary);
        std::ifstream dst3("cdata.old", std::ios::binary);
        dst3 << src3.rdbuf();

        // Delete cdata        
        src3.close();
        src3.open("cdata.txt", std::ofstream::out | std::ofstream::trunc);
        src3.close();
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
