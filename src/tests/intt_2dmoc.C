#include "TestDefinitions.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    RCP<Epetra_Comm> comm;
    RCP<Ocean> ocean;
}

//------------------------------------------------------------------
TEST(Ocean, Continuation)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        RCP<Teuchos::ParameterList> oceanParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));

        // Create continuation params
        RCP<Teuchos::ParameterList> continuationParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params.xml",
                                    continuationParams.ptr());

        {
            // Create contination
            Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
                continuation(ocean, continuationParams);

            // Run continuation in the forcing
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "CMPR");
        continuationParams->set("destination 0", -0.2);
        continuationParams->set("initial step size", -0.5);

        {
            // Create contination
            Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
                continuation(ocean, continuationParams);

            // Run continuation in the asymmetry parameter
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "Salinity Forcing");
        continuationParams->set("destination 0", 0.02);
        continuationParams->set("initial step size", 0.5);

        {
            // Create contination
            Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
                continuation(ocean, continuationParams);

            // Run continuation in the salinity forcing
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "CMPR");
        continuationParams->set("destination 0", 0.0);
        continuationParams->set("initial step size", 0.5);

        {
            // Create contination
            Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
                continuation(ocean, continuationParams);

            // Run continuation in the asymmetry parameter to get back to the symmetric state
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "Salinity Forcing");
        continuationParams->set("destination 0", 0.04);
        continuationParams->set("initial step size", 0.5);

        {
            // Create contination
            Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
                continuation(ocean, continuationParams);

            // Run continuation in the salinity forcing to state 1
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        double psiMin1, psiMax1;
        ocean->getPsiM(psiMin1, psiMax1);
        EXPECT_NEAR(psiMax1, 14.7, 1e-1);
        EXPECT_NEAR(psiMin1, 0, 1e-4);

        continuationParams->set("continuation parameter", "Salinity Forcing");
        continuationParams->set("destination 0", 0.03);
        continuationParams->set("initial step size", -0.5);

        {
            // Create contination
            Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
                continuation(ocean, continuationParams);

            // Run continuation in the salinity forcing
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "Salinity Forcing");
        continuationParams->set("destination 0", 0.04);
        continuationParams->set("initial step size", -0.5);

        {
            // Create contination
            Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
                continuation(ocean, continuationParams);

            // Run continuation in the salinity forcing to state 2
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        double psiMin2, psiMax2;
        ocean->getPsiM(psiMin2, psiMax2);
        EXPECT_NEAR(psiMin1, -psiMax2, 1e-4);
        EXPECT_NEAR(psiMax1, -psiMin2, 1e-4);
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

    if (comm->MyPID() == 0)
        printProfile();

    MPI_Finalize();
    return out;
}
