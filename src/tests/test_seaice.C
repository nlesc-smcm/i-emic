#include "TestDefinitions.H"

// Testing standalone seaice

//------------------------------------------------------------------
namespace
{
    std::shared_ptr<SeaIce>     seaIce;
    RCP<Epetra_Comm>            comm;
    RCP<Teuchos::ParameterList> seaIceParams;
}

//------------------------------------------------------------------
TEST(SeaIce, Initialization)
{
    bool failed = false;

    // Create seaice parameters
    seaIceParams = rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("seaice_params.xml", seaIceParams.ptr());

    try
    {
        seaIce = std::make_shared<SeaIce>(comm, seaIceParams);
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
// compute rhs under idealized conditions
TEST(SeaIce, computeRHS)
{
    Teuchos::RCP<Epetra_Vector> state = seaIce->getState('V');
    Teuchos::RCP<Epetra_Vector> rhs   = seaIce->getRHS('V');
    seaIce->idealizedForcing();
    state->PutScalar(0.0);     
    seaIce->computeRHS();

    // we know a few norms for the idealized case
    EXPECT_NEAR(Utils::norm(rhs), 6.1874305776402e+03, 1e-7);

    state->PutScalar(1.0);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), 4.825751331172339e+04, 1e-7);

    state->PutScalar(1.234);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), 5.810287606919227e+04, 1e-7);
}

//------------------------------------------------------------------
TEST(SeaIce, computeJacobian)
{
    Teuchos::RCP<Epetra_Vector> state = seaIce->getState('V');
    state->PutScalar(0.0);
    seaIce->computeJacobian();
    
    Teuchos::RCP<Epetra_CrsMatrix> jac = seaIce->getJacobian();
   
    DUMPMATLAB("seaice_jac", *jac);
    double normOne, normInf, normFrob;
    normInf = jac->NormInf();
    EXPECT_NEAR(normInf, 2.712570151601840e+03, 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne, 2.722944891390629e+03, 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, 4.275055233767215e+04, 1e-7);

    state->PutScalar(1.0);
    seaIce->computeJacobian();
    jac = seaIce->getJacobian();

    normInf = jac->NormInf();
    EXPECT_NEAR(normInf,  2.712570151601840e+03, 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne,  2.672944891390629e+03, 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, 4.274306639884469e+04, 1e-7);

    state->PutScalar(1.234);
    seaIce->computeJacobian();
    jac = seaIce->getJacobian();

    normInf = jac->NormInf();
    EXPECT_NEAR(normInf, 2.712570151601840e+03, 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne, 2.672944891390629e+03, 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, 4.274306639884469e+04, 1e-7);    
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
    seaIce = std::shared_ptr<SeaIce>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
