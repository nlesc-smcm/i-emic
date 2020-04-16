#include "TestDefinitions.H"

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Ocean.H"
#include "Continuation.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Epetra_Comm> comm;
    Teuchos::RCP<Ocean> ocean;
}

//------------------------------------------------------------------
TEST(Ocean, Continuation1)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        Teuchos::RCP<Teuchos::ParameterList> oceanParams =
            Utils::obtainParams("ocean_params_intt_init.xml", "Ocean parameters");
        oceanParams->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));

        // Create continuation params
        Teuchos::RCP<Teuchos::ParameterList> continuationParams =
            Teuchos::rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params_intt_init.xml",
                                    continuationParams.ptr());

        // Create contination
        Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

        // Run continuation
        int status = continuation.run();
        EXPECT_EQ(status, 0);

        // After the run get the Jacobian matrix
        Teuchos::RCP<Epetra_CrsMatrix> jac = ocean->getJacobian();
        DUMPMATLAB("ocean_jac_init", *jac);

        // Also get the mass matrix
        Teuchos::RCP<Epetra_Vector> diagB = ocean->getMassMat();
        EXPECT_NE(Utils::norm(diagB), 0.0);
        DUMP_VECTOR("ocean_B_init", *diagB);


        // Copy ocean_output.h5 to ocean_input.h5
        std::ifstream src1("ocean_output.h5", std::ios::binary);
        std::ofstream dst1("ocean_input.h5", std::ios::binary);
        dst1 << src1.rdbuf();

        // Copy info
        std::ifstream src2("info_0.txt", std::ios::binary);
        std::ofstream dst2("info_0.old", std::ios::binary);
        dst2 << src2.rdbuf();

        // Copy cdata
        std::ifstream src3("cdata.txt", std::ios::binary);
        std::ofstream dst3("cdata.old", std::ios::binary);
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
        Teuchos::RCP<Teuchos::ParameterList> oceanParams =
            Utils::obtainParams("ocean_params_intt_cont.xml", "Ocean parameters");
        oceanParams->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");
        ocean = Teuchos::null;
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));

        // Create continuation params
        Teuchos::RCP<Teuchos::ParameterList> continuationParams =
            Teuchos::rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params_intt_cont.xml",
                                    continuationParams.ptr());

        // Create contination
        Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

        // Run continuation
        int status = continuation.run();
        EXPECT_EQ(status, 0);

        // After the run get the Jacobian matrix
        Teuchos::RCP<Epetra_CrsMatrix> jac = ocean->getJacobian();
        DUMPMATLAB("ocean_jac_cont", *jac);

        // Also get the mass matrix
        Teuchos::RCP<Epetra_Vector> diagB = ocean->getMassMat();
        EXPECT_NE(Utils::norm(diagB), 0.0);
        DUMP_VECTOR("ocean_B_cont", *diagB);

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
