#include "TestDefinitions.H"

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Topo.H"
#include "Ocean.H"
#include "Continuation.H"
#include "Utils.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Ocean>  ocean;
    Teuchos::RCP<Topo<Teuchos::RCP<Ocean>, Teuchos::RCP<Teuchos::ParameterList> > > topo;
    Teuchos::RCP<Epetra_Comm>  comm;
}

//------------------------------------------------------------------
TEST(Ocean, Initialization)
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
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Topo, Initialization)
{
    bool failed = false;
    try
    {
        // Create topography class
        Teuchos::RCP<Teuchos::ParameterList> topoParams =
            Utils::obtainParams("topo_params.xml", "Topo parameters");

        topoParams->sublist("Belos solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");
        topo = Teuchos::rcp(new Topo<Teuchos::RCP<Ocean>, Teuchos::RCP<Teuchos::ParameterList> >
                   (ocean, topoParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Topo, Arrays)
{
    std::vector<int> a = topo->getA();
    std::vector<int> b = topo->getB();
    int N = topo->nMasks();

    std::vector<int> aref = {0,1,2,3,4,5,6,7};
    std::vector<int> bref = {1,2,3,4,5,6,7,8};

    for (int i = 0; i != std::min(N, 8); ++i)
    {
        EXPECT_EQ(a[i], aref[i]);
        EXPECT_EQ(b[i], bref[i]);
    }
}

//------------------------------------------------------------------
TEST(Topo, Copy)
{
    double tol = 1e-12;
    double nrmC, nrmV;
    Ocean::VectorPtr solCopy = topo->getSolution('C');
    Ocean::VectorPtr solView = topo->getSolution('V');
    nrmC = Utils::norm(solCopy);
    nrmV = Utils::norm(solView);
    EXPECT_NEAR(nrmC, nrmV, tol);

    solCopy->PutScalar(3.14);
    nrmC = Utils::norm(solCopy);
    nrmV = Utils::norm(solView);

    EXPECT_NE(nrmC, nrmV);
}

//------------------------------------------------------------------
TEST(Topo, View)
{
    double tol = 1e-12;
    double nrmC, nrmV;
    Ocean::VectorPtr solCopy = topo->getSolution('C');
    Ocean::VectorPtr solView = topo->getSolution('V');
    nrmC = Utils::norm(solCopy);
    nrmV = Utils::norm(solView);
    EXPECT_NEAR(nrmC, nrmV, tol);

    solView->PutScalar(3.14);

    nrmC = Utils::norm(solCopy);
    nrmV = Utils::norm(solView);

    solView->PutScalar(0.0);

    EXPECT_NE(nrmC, nrmV);
}

//------------------------------------------------------------------
TEST(Topo, SpinupContinuation)
{
    bool failed = false;
    try
    {
        // not sure if needed
        ocean->setPar("Combined Forcing", 0.0);
        ocean->getState('V')->PutScalar(0.0);
        topo->getSolution('V')->PutScalar(0.0); // superfluous

        // Create parameter object for continuation
        Teuchos::RCP<Teuchos::ParameterList> continuationParams =
            Teuchos::rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("spinup_continuation_params.xml",
                                    continuationParams.ptr());

        // Create spinup continuation
        Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

        // Run continuation
        INFO("--**-- Topo test: running spinup...");
        int status = continuation.run();
        EXPECT_EQ(status, 0);
        INFO("--**-- Topo test: running spinup... done");

        Utils::MaskStruct mask = ocean->getLandMask();
        Utils::printSurfaceMask(mask.global_surface, "mask_before", mask.n);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Topo, RHS)
{
    bool failed = false;
    try
    {
        topo->setPar("Delta", 0.3);
        topo->computeRHS();
        topo->computeJacobian();
        topo->preProcess();
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Topo, TopoContinuation)
{
    bool failed = false;

    // ocean norm spinup topography
    int nMasks;
    std::vector<double> norms;
    try
    {
        // not sure if needed
        topo->setPar("Delta", 0.0);

        // Create parameter object for continuation
        Teuchos::RCP<Teuchos::ParameterList> continuationParams =
            Teuchos::rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("topo_continuation_params.xml",
                                    continuationParams.ptr());

        // Create topo continuation
        Continuation<Teuchos::RCP<Topo<Teuchos::RCP<Ocean>, Teuchos::RCP<Teuchos::ParameterList> > >>
            continuation(topo, continuationParams);

        // Run topo continuation
        nMasks    = topo->nMasks();
        int startMask = topo->startMaskIdx();

        INFO(" Running topo cont...");// Run continuation

        // We do a couple of homotopy continuations and expect to end
        // up with the original landmask. We do mask 0..4,
        // where 0,4 and 1,3 are equal (see topo_params.xml).
        norms = std::vector<double>(nMasks);
        norms[0] = Utils::norm(ocean->getState('V'));
        for (int maskIdx = startMask; maskIdx != nMasks-1; maskIdx++)
        {
            topo->setMaskIndex(maskIdx);
            topo->reinitialize();
            topo->predictor();

            int status = continuation.run();
            EXPECT_EQ(status, 0);

            norms[maskIdx+1] = Utils::norm(ocean->getState('V'));
        }
        INFO(" Running topo cont... done");
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);

    for (int k = 0; k < std::floor(nMasks / 2); ++k)
    {
        // test equal topos
        EXPECT_LT(std::abs(norms[k]-norms[nMasks-k-1]), 1e-7);

        // test different (subsequent) topos
        EXPECT_GT(std::abs(norms[k+1]-norms[k]), 1e-7);
        EXPECT_GT(std::abs(norms[nMasks-k-1]-norms[nMasks-k-2]), 1e-7);
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
    topo  = Teuchos::null;

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
