#include "TestDefinitions.H"

#include "Ocean.H"
#include "Atmosphere.H"
#include "CoupledModel.H"
#include "Combined_MultiVec.H"
#include "Continuation.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Epetra_Comm>      comm;
    std::shared_ptr<Ocean>         ocean;
    std::shared_ptr<Atmosphere>    atmos;
    std::shared_ptr<CoupledModel>  coupledModel;
    std::vector<Teuchos::RCP<Teuchos::ParameterList> > params;
    enum Ident { OCEAN, ATMOS, COUPLED, CONT};
    std::shared_ptr<Combined_MultiVec> state1, state2;
}

//------------------------------------------------------------------
TEST(ParameterLists, Initialization)
{
    bool failed = false;
    try
    {
        std::vector<std::string> files = {"test_oa_ocean.xml",
                                          "test_oa_atmos.xml",
                                          "test_oa_cpldm.xml",
                                          "test_oa_contn.xml"};

        std::vector<std::string> names = {"Ocean parameters",
                                          "Atmosphere parameters",
                                          "CoupledModel parameters",
                                          "Continuation parameters"};

        for (int i = 0; i != (int) files.size(); ++i)
            params.push_back(Utils::obtainParams(files[i], names[i]));

        params[OCEAN]->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");

        INFO('\n' << "Overwriting:");

        // Allow dominant parameterlists. Note that this trick uses a
        // 'flattened' hierarchy. The Continuation and CoupledModel
        // parameterlists are allowed to overwrite settings.
        Utils::overwriteParameters(params[OCEAN],  params[COUPLED]);
        Utils::overwriteParameters(params[ATMOS],  params[COUPLED]);

        Utils::overwriteParameters(params[OCEAN],  params[CONT]);
        Utils::overwriteParameters(params[ATMOS],  params[CONT]);

        Utils::overwriteParameters(params[COUPLED], params[CONT]);
        INFO('\n');
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed,false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Initialization_1)
{
    bool failed = false;
    try
    {
        // Specify landmask
        params[OCEAN]->sublist("THCM").set("Land Mask", "test6x12x4_1");
        // Create parallel Ocean
        ocean = std::make_shared<Ocean>(comm, params[OCEAN]);
        // Create atmosphere
        atmos = std::make_shared<Atmosphere>(comm, params[ATMOS]);
        // Create coupledmodel
        coupledModel = std::make_shared<CoupledModel>(ocean,
                                                      atmos,
                                                      params[COUPLED]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Continuation_1)
{
    bool failed = false;
    try
    {
        coupledModel->setPar("Combined Forcing", 0.0);
        coupledModel->initializeState();

        // Create continuation
        Continuation<std::shared_ptr<CoupledModel>>
            continuation(coupledModel, params[CONT]);

        // Run continuation
        int status = continuation.run();
        EXPECT_EQ(status, 0);

        state1 = coupledModel->getState('C');
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Initialization_2)
{
    bool failed = false;
    try
    {
        // Specify landmask
        params[OCEAN]->sublist("THCM").set("Land Mask", "test6x12x4_2");
        // Destroy everything
        ocean        = std::shared_ptr<Ocean>();
        atmos        = std::shared_ptr<Atmosphere>();
        coupledModel = std::shared_ptr<CoupledModel>();

        // Create everything again
        ocean = std::make_shared<Ocean>(comm, params[OCEAN]);
        atmos = std::make_shared<Atmosphere>(comm, params[ATMOS]);
        coupledModel = std::make_shared<CoupledModel>(ocean,
                                                      atmos,
                                                      params[COUPLED]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Continuation_2)
{
    bool failed = false;
    try
    {
        coupledModel->setPar("Combined Forcing", 0.0);
        coupledModel->initializeState();

        // Create continuation
        Continuation<std::shared_ptr<CoupledModel>>
            continuation(coupledModel, params[CONT]);

        // Run continuation
        int status = continuation.run();
        EXPECT_EQ(status, 0);

        state2 = coupledModel->getState('C');

    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Ocean, Periodicity)
{
    std::vector<int> unknowns = {1,2,5,6};
    double norm1, norm2;
    for (auto &i : unknowns)
    {
        norm1 = Utils::normOfField((*state1)(OCEAN), ocean->getDomain(), i);
        norm2 = Utils::normOfField((*state2)(OCEAN), ocean->getDomain(), i);
        EXPECT_NEAR(norm1, norm2, 1e-4);
    }
}

//------------------------------------------------------------------
TEST(Atmosphere, Periodicity)
{
    std::vector<int> unknowns = {1,2};
    double norm1, norm2;
    for (auto &i : unknowns)
    {
        norm1 = Utils::normOfField((*state1)(ATMOS), atmos->getDomain(), i);
        norm2 = Utils::normOfField((*state2)(ATMOS), atmos->getDomain(), i);
        EXPECT_NEAR(norm1, norm2, 1e-4);
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
    ocean        = std::shared_ptr<Ocean>();
    atmos        = std::shared_ptr<Atmosphere>();
    coupledModel = std::shared_ptr<CoupledModel>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
