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
    seaIce->setPar("Combined Forcing", 1.0);
    seaIce->computeRHS();

    // we know a few norms for the idealized case
    EXPECT_NEAR(Utils::norm(rhs), 6.029707151443e+03, 1e-7);

    state->PutScalar(1.0);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), 3.375072979252e+04, 1e-7);

    state->PutScalar(1.234);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), 4.141771047244e+04, 1e-7);
    seaIce->setPar("Combined Forcing", 0.1);
    state->PutScalar(0.0);        
    seaIce->computeRHS();

    // we know a few norms for the idealized case
    EXPECT_NEAR(Utils::norm(rhs), 6.029707151443e+02, 1e-7);

    state->PutScalar(1.0);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), 3.346451471720e+03, 1e-7);

    state->PutScalar(1.234);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), 4.130427375713e+03, 1e-7);
}


//------------------------------------------------------------------
TEST(SeaIce, computeJacobian)
{
    Teuchos::RCP<Epetra_Vector> state = seaIce->getState('V');
    state->PutScalar(0.0);
    seaIce->setPar("Combined Forcing", 1.0);        
    seaIce->computeJacobian();
    
    Teuchos::RCP<Epetra_CrsMatrix> jac = seaIce->getJacobian();
   
    DUMPMATLAB("seaice_jac", *jac);
    double normOne, normInf, normFrob;
    normInf = jac->NormInf();
    EXPECT_NEAR(normInf, 2.079199128921e+03, 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne, 2.127721683616e+03, 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, 3.322972397910e+04, 1e-7);

    state->PutScalar(1.0);
    seaIce->computeJacobian();
    jac = seaIce->getJacobian();

    normInf = jac->NormInf();
    EXPECT_NEAR(normInf,  2.079199128921e+03, 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne,  2.077721683616e+03, 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, 3.322009265079e+04, 1e-7);

    state->PutScalar(1.234);
    seaIce->computeJacobian();
    jac = seaIce->getJacobian();

    normInf = jac->NormInf();
    EXPECT_NEAR(normInf, 2.079199128921e+03, 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne, 2.077721683616e+03, 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, 3.322009265079e+04, 1e-7);

    state->PutScalar(0.0);
    seaIce->setPar("Combined Forcing", 0.1);        
    seaIce->computeJacobian();
    
    normInf = jac->NormInf();
    EXPECT_NEAR(normInf, 2.088199128921e+02, 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne, 2.144154293601e+02, 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, 3.379287916717e+03, 1e-7);
}


//------------------------------------------------------------------
TEST(SeaIce, Solve)
{
    Teuchos::RCP<Epetra_Vector> state = seaIce->getState('V');
    Teuchos::RCP<Epetra_Vector> b, x, r;
    double rNorm, bNorm, xNorm;
    
    for (int i = 0; i != 5; ++i)
    {
        state->Random();
        
        seaIce->computeJacobian();
        seaIce->computeRHS();

        b = seaIce->getRHS('C');

        seaIce->solve(b);
        
        x = seaIce->getSolution('C');
        r = seaIce->getSolution('C');

        r->PutScalar(0.0);

        seaIce->applyMatrix(*x, *r);
        r->Update(1.0, *b, -1.0);

        rNorm = Utils::norm(r);
        bNorm = Utils::norm(b);
        xNorm = Utils::norm(x);

        std::cout << " || b - Ax || = " << rNorm << std::endl;
        std::cout << "        ||b|| = " << bNorm << std::endl;
        std::cout << "        ||x|| = " << xNorm << std::endl;
    
        EXPECT_NEAR(rNorm, 0, 1e-4);
    }
}


//------------------------------------------------------------------
TEST(SeaIce, Newton)
{
    Teuchos::RCP<Epetra_Vector> state = seaIce->getState('V');
    Teuchos::RCP<Epetra_Vector> b, x, r;

    state->Random();
    state->Scale(100);
    seaIce->computeRHS();
    for (int i = 0; i != 5; ++i)
    {
        seaIce->computeJacobian();
        b = seaIce->getRHS('C');
        b->Scale(-1.0);
        seaIce->solve(b);
        x = seaIce->getSolution('V');
        state->Update(1.0, *x, 1.0);
        seaIce->computeRHS();
        std::cout << "||dx|| = " << Utils::norm(x) << std::endl;
    }
    EXPECT_LT(Utils::norm(x), 1e-8);
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
