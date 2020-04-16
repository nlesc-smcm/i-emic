#include "TestDefinitions.H"

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Continuation.H"
#include "Ocean.H"

#include "TRIOS_Domain.H"

#include "Epetra_Import.h"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Teuchos::ParameterList> oceanParams;
    Teuchos::RCP<Ocean> ocean;
    Teuchos::RCP<Epetra_Comm> comm;
}

//------------------------------------------------------------------
TEST(Ocean, Initialization)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        oceanParams = Utils::obtainParams("reft_ocean_params.xml", "Ocean parameters");
        oceanParams->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Ocean, Continuation)
{
    // Create continuation params
    Teuchos::RCP<Teuchos::ParameterList> continuationParams =
        rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("reft_continuation_params.xml",
                                continuationParams.ptr());

    // Create contination
    Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

    // Run continuation
    int status = continuation.run();
    EXPECT_EQ(status, 0);
}

//------------------------------------------------------------------
TEST(Ocean, CheckReference)
{
    // create vectors and load reference solution
    Teuchos::RCP<Epetra_Vector> x = ocean->getState();
    Teuchos::RCP<Epetra_Vector> y = ocean->getState();
    y->PutScalar(0.0);
    Utils::load(y, "ocean_reference");

    // test values of fields that are not determined up to a constant
    // (u,v,T,S)
    std::vector<int> unknowns = {1, 2, 5, 6};
    Teuchos::RCP<Epetra_Map> rowMap = ocean->getDomain()->GetSolveMap();
    Teuchos::RCP<Epetra_Map> mapid;
    Teuchos::RCP<Epetra_Import> importid;
    for (auto &id: unknowns)
    {
        mapid = Utils::CreateSubMap(*rowMap, 6, id);
        importid =  Teuchos::rcp(new Epetra_Import(*rowMap, *mapid));

        Epetra_Vector xid(*mapid);
        CHECK_ZERO(xid.Export(*x, *importid, Zero));
        Epetra_Vector yid(*mapid);
        CHECK_ZERO(yid.Export(*y, *importid, Zero));

        double normxid, normyid;
        xid.Norm2(&normxid);
        yid.Norm2(&normyid);

        EXPECT_NEAR(normxid, normyid, 1e-3);
    }
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
