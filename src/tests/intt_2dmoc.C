#include "TestDefinitions.H"

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Continuation.H"
#include "Ocean.H"
#include "Utils.H"

#ifdef HAVE_RAILS
#include "LyapunovModel.H"
#endif

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Epetra_Comm> comm;
    Teuchos::RCP<Ocean> ocean;
}

//------------------------------------------------------------------
TEST(Ocean, Continuation)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        Teuchos::RCP<Teuchos::ParameterList> oceanParams =
            Utils::obtainParams("ocean_params.xml", "Ocean parameters");
        oceanParams->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));

        // Create continuation params
        Teuchos::RCP<Teuchos::ParameterList> continuationParams =
            Teuchos::rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params.xml",
                                    continuationParams.ptr());

        {
            // Create contination
            Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

            // Run continuation in the forcing
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "CMPR");
        continuationParams->set("destination 0", -0.2);
        continuationParams->set("initial step size", -0.5);

        {
            // Create contination
            Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

            // Run continuation in the asymmetry parameter
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "Salinity Forcing");
        continuationParams->set("destination 0", 0.02);
        continuationParams->set("initial step size", 0.5);

        {
            // Create contination
            Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

            // Run continuation in the salinity forcing
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "CMPR");
        continuationParams->set("destination 0", 0.0);
        continuationParams->set("initial step size", 0.5);

        {
            // Create contination
            Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

            // Run continuation in the asymmetry parameter to get back to the symmetric state
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "Salinity Forcing");
        continuationParams->set("destination 0", 0.04);
        continuationParams->set("initial step size", 0.5);

        {
            // Create contination
            Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

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
            Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

            // Run continuation in the salinity forcing
            int status = continuation.run();
            EXPECT_EQ(status, 0);
        }

        continuationParams->set("continuation parameter", "Salinity Forcing");
        continuationParams->set("destination 0", 0.04);
        continuationParams->set("initial step size", -0.5);

        {
            // Create contination
            Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

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
TEST(Ocean, Lyapunov)
{
    bool failed = false;
    try
    {
#ifdef HAVE_RAILS
        std::vector<double> ev, ev2;

        // Create the Ocean object wrapped in a Lyapunov model
        Teuchos::RCP<LyapunovModel<Ocean> > lyap = Teuchos::rcp(
            new LyapunovModel<Ocean>(*ocean));

        lyap->computeCovarianceMatrix();
        ev = lyap->getEigenvalues();
        EXPECT_NEAR(ev[0], 10.08, 1e-2);
        EXPECT_NEAR(ev[1], 1.76, 1e-2);

        // Create continuation params
        Teuchos::RCP<Teuchos::ParameterList> continuationParams =
            Teuchos::rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params.xml",
                                    continuationParams.ptr());

        continuationParams->set("continuation parameter", "Salinity Forcing");
        continuationParams->set("destination 0", 0.041);
        continuationParams->set("initial step size", 0.1);

        // Create continuation
        Continuation<Teuchos::RCP<LyapunovModel<Ocean>>> continuation(lyap, continuationParams);

        // Run continuation
        int status = continuation.run();
        assert(status == 0);

        // Check the new eigenvalues
        ev = lyap->getEigenvalues();
        EXPECT_NEAR(ev[0], 21.81, 1e-2);
        EXPECT_NEAR(ev[1], 4.09, 1e-2);

        // Check that this was already computed during the continuation
        lyap->computeCovarianceMatrix();
        ev2 = lyap->getEigenvalues();
        EXPECT_NEAR(ev[0], ev2[0], 1e-2);
        EXPECT_NEAR(ev[1], ev2[1], 1e-2);
#else
#ifdef GTEST_SKIP
        GTEST_SKIP();
#else
        WARNING("RAILS not found", __FILE__, __LINE__);
#endif
#endif
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
