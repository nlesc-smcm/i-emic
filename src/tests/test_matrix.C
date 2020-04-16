#include "TestDefinitions.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Ocean.H"
#include "Atmosphere.H"
#include "SeaIce.H"
#include "CoupledModel.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Epetra_Comm>      comm;

    std::shared_ptr<Ocean>         ocean;
    std::shared_ptr<Atmosphere>    atmos;
    std::shared_ptr<SeaIce>        seaice;
    std::shared_ptr<CoupledModel>  coupledModel;
    std::vector<Teuchos::RCP<Teuchos::ParameterList> > params;
    enum Ident { OCEAN, ATMOS, SEAICE, COUPLED, CONT};
}

//------------------------------------------------------------------
TEST(ParameterLists, Initialization)
{
    bool failed = false;
    try
    {
        std::vector<std::string> files = {"ocean_params.xml",
                                          "atmosphere_params.xml",
                                          "seaice_params.xml",
                                          "coupledmodel_params.xml",
                                          "continuation_params.xml"};

        std::vector<std::string> names = {"Ocean parameters",
                                          "Atmosphere parameters",
                                          "Sea ice parameters",
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
        Utils::overwriteParameters(params[SEAICE], params[COUPLED]);

        Utils::overwriteParameters(params[OCEAN],  params[CONT]);
        Utils::overwriteParameters(params[ATMOS],  params[CONT]);
        Utils::overwriteParameters(params[SEAICE], params[CONT]);

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
TEST(Ocean, Initialization)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        ocean = std::make_shared<Ocean>(comm, params[OCEAN]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, Initialization)
{
    bool failed = false;
    try
    {
        // Create atmosphere
        atmos = std::make_shared<Atmosphere>(comm, params[ATMOS]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(SeaIce, Initialization)
{
    bool failed = false;
    try
    {
        // Create atmosphere
        seaice = std::make_shared<SeaIce>(comm, params[SEAICE]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Initialization)
{
    bool failed = false;
    try
    {
        // Create coupledmodel
        coupledModel = std::make_shared<CoupledModel>(ocean,
                                                      atmos,
                                                      seaice,
                                                      params[COUPLED]);

        coupledModel->initializeState();
        coupledModel->postProcess();
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
// Apply the vector and store. Load vectors from other runs with other
// #cores and compare their norm.
TEST(CoupledModel, Matrix)
{
    std::shared_ptr<Combined_MultiVec> s = coupledModel->getState('V');
    // s->SetSeed(22);
    // s->Random();
    // s->Scale(0.01);
    INFO("   ||s|| = " << Utils::norm(s));
    coupledModel->setPar("Combined Forcing", 0.01);
    coupledModel->computeJacobian();

    std::shared_ptr<Combined_MultiVec> x = coupledModel->getSolution('C');
    std::shared_ptr<Combined_MultiVec> y = coupledModel->getSolution('C');
    std::shared_ptr<Combined_MultiVec> o = coupledModel->getSolution('C');
    x->PutScalar(0.0);
    o->PutScalar(1.0);

    coupledModel->applyMatrix(*o, *x);

    double normx = Utils::norm(x);

    INFO(" # cores = " << comm->NumProc());
    INFO("   ||x|| = " << normx);

    std::stringstream fname_save;
    fname_save << "MV_procs" << comm->NumProc() << "_";
    Utils::save(x, fname_save.str());

    std::vector<int> procs = {1,2,4,8};
    for (auto &p : procs)
    {
        if (p >= comm->NumProc())
            break;

        std::stringstream fname_load, diff;
        fname_load << "MV_procs" << p << "_";

        Utils::load(y, fname_load.str());
        INFO("   ||y|| = " << Utils::norm(y));
        EXPECT_NEAR(Utils::norm(y), normx, 1e-8);

        y->Update(-1.0,*x,1.0);
        diff << "diff_" << comm->NumProc() << "-" << p;
        Utils::save(y, diff.str());
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
    seaice       = std::shared_ptr<SeaIce>();
    coupledModel = std::shared_ptr<CoupledModel>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
