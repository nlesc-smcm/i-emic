#include "TestDefinitions.H"

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "NumericalJacobian.H"
#include "SeaIce.H"

// Testing standalone seaice

//------------------------------------------------------------------
namespace
{
    std::shared_ptr<SeaIce>              seaIce;
    Teuchos::RCP<Epetra_Comm>            comm;
    Teuchos::RCP<Teuchos::ParameterList> seaIceParams;
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

/*
//------------------------------------------------------------------
// compute rhs under idealized conditions
TEST(SeaIce, computeRHS)
{
    Teuchos::RCP<Epetra_Vector> state = seaIce->getState('V');
    Teuchos::RCP<Epetra_Vector> rhs   = seaIce->getRHS('V');
    seaIce->idealizedForcing();
    state->PutScalar(0.0);
    
    std::vector<double> testResults = {
        1.783263089276e+02,
        6.165632554697e+02,
        7.289000573802e+02,
        1.714310465387e+02,
        6.012200920047e+02,
        7.110217880537e+02
    };
    
    seaIce->setPar("Combined Forcing", 1.0);
    seaIce->computeRHS();

    // we know a few norms for the idealized case
    EXPECT_NEAR(Utils::norm(rhs), testResults[0], 1e-7);

    state->PutScalar(1.0);     
    seaIce->computeRHS();
    
    EXPECT_NEAR(Utils::norm(rhs), testResults[1], 1e-7);
                                              
    state->PutScalar(1.234);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), testResults[2], 1e-7);

    seaIce->setPar("Combined Forcing", 0.1);
    state->PutScalar(0.0);
    seaIce->computeRHS();

    // we know a few norms for the idealized case
    EXPECT_NEAR(Utils::norm(rhs), testResults[3], 1e-7);
    
    state->PutScalar(1.0);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), testResults[4], 1e-7);

    state->PutScalar(1.234);     
    seaIce->computeRHS();

    EXPECT_NEAR(Utils::norm(rhs), testResults[5], 1e-7);
}

//------------------------------------------------------------------
TEST(SeaIce, computeJacobian)
{
    Teuchos::RCP<Epetra_Vector> state = seaIce->getState('V');
    state->PutScalar(0.0);
    seaIce->setPar("Combined Forcing", 1.0);        
    seaIce->computeJacobian();
    Teuchos::RCP<Epetra_CrsMatrix> jac = seaIce->getJacobian();

    std::vector<double> testResults = {
        5.099999983333e+01,
        9.650722444261e+01,
        6.264955752426e+02,
        4.947134256694e+01,
        4.650722460927e+01,
        4.821791232159e+02,
        4.947134256694e+01,
        4.650722460927e+01,
        4.821791232159e+02,
        5.099999983333e+01,
        9.650722444261e+01,
        6.259851360683e+02,
        4.947134256694e+01,
        4.650722460927e+01,
        4.815157231546e+02,
        4.947134256694e+01,
        4.650722460927e+01,
        4.815157231546e+02
    };
    
    double normOne, normInf, normFrob;
    normInf = jac->NormInf();
    EXPECT_NEAR(normInf, testResults[0], 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne, testResults[1], 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, testResults[2], 1e-7);

    state->PutScalar(1.0);
    seaIce->computeJacobian();
    jac = seaIce->getJacobian();

    normInf = jac->NormInf();
    EXPECT_NEAR(normInf, testResults[3], 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne, testResults[4], 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, testResults[5], 1e-7);

    state->PutScalar(1.234);
    seaIce->computeJacobian();
    jac = seaIce->getJacobian();

    normInf = jac->NormInf();
    EXPECT_NEAR(normInf, testResults[6], 1e-7);

    normOne = jac->NormOne();
    EXPECT_NEAR(normOne, testResults[7], 1e-7);

    normFrob = jac->NormFrobenius();
    EXPECT_NEAR(normFrob, testResults[8], 1e-7);

    // state->PutScalar(0.0);
    // seaIce->setPar("Combined Forcing", 0.1);        
    // seaIce->computeJacobian();
    
    // normInf = jac->NormInf();
    // EXPECT_NEAR(normInf, testResults[9], 1e-7);

    // normOne = jac->NormOne();
    // EXPECT_NEAR(normOne, testResults[10], 1e-7);

    // normFrob = jac->NormFrobenius();
    // EXPECT_NEAR(normFrob, testResults[11], 1e-7);

}
*/

//------------------------------------------------------------------
TEST(SeaIce, numericalJacobian)
{
    seaIce->setPar("Combined Forcing", 0.85);

    Teuchos::RCP<Epetra_Vector> x = seaIce->getState('V');

    x->SetSeed(1.0);
    x->Random();

    // int myEl = x->Map().NumMyElements();
    // x->PutScalar(0.0);
    // // gradually changing thickness deviation
    // for (int i = SEAICE_HH_-1; i < myEl; i += SEAICE_NUN_)
    //     (*x)[i] = -.5 + ( (double) i / (myEl-1) );

    seaIce->computeRHS();
    seaIce->computeJacobian();

    Teuchos::RCP<Epetra_CrsMatrix> seaIceJac  = seaIce->getJacobian();
    DUMPMATLAB("seaIceJac", *seaIceJac);

    // only do this test for small problems in serial
    int nmax = 2e3;

    if ( (comm->NumProc() == 1) &&
         (seaIce->getState('V')->GlobalLength() < nmax) )
    {
        bool failed = false;
        try
        {
            INFO("compute njC");
            NumericalJacobian<std::shared_ptr<SeaIce>,
                              Teuchos::RCP<Epetra_Vector> > numJac;

            numJac.setTolerance(1e-12);
            numJac.seth(1e-6);
            numJac.compute(seaIce, seaIce->getState('V'));
            numJac.print("seaIceNumJac");

            NumericalJacobian<std::shared_ptr<SeaIce>,
                              Teuchos::RCP<Epetra_Vector> >::CCS ccs;
            numJac.fillCCS(ccs);

            EXPECT_NE(ccs.beg.back(), 0);

            Teuchos::RCP<Epetra_Vector> x = seaIce->getState('C');

            testEntries(seaIce, ccs, x);

        }
        catch (...)
        {
            failed = true;
            throw;
        }


        EXPECT_EQ(failed, false);
    }
    if (comm->NumProc() != 1)
    {
        std::cout << ("****Numerical Jacobian test cannot run in parallel****\n") ;
        INFO("****Numerical Jacobian test cannot run in parallel****");
    }

    if (seaIce->getState('V')->GlobalLength() > nmax)
    {
        std::cout << ("****Numerical Jacobian test cannot run for this problem size****\n");
        INFO("****Numerical Jacobian test cannot run for this problem size****");
    }
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
    seaIce->computeRHS();
    for (int i = 0; i != 15; ++i)
    {
        seaIce->computeJacobian();
        b = seaIce->getRHS('C');
        b->Scale(-1.0);
        seaIce->solve(b);
        x = seaIce->getSolution('V');
        state->Update(1.0, *x, 1.0);
        seaIce->computeRHS();
        std::cout.precision(3);

        std::cout << "i = " << std::setw(2) << i
                  << ", ||dx||=" << std::setw(8) << std::scientific << Utils::norm(x)
                  << ", ||H|| =" << std::setw(8) << std::scientific << Utils::norm(seaIce->interface(x, SEAICE_HH_))
                  << ", ||Q|| =" << std::setw(8) << std::scientific << Utils::norm(seaIce->interface(x, SEAICE_QQ_))
                  << ", ||M|| =" << std::setw(8) << std::scientific << Utils::norm(seaIce->interface(x, SEAICE_MM_))
                  << ", ||T|| =" << std::setw(8) << std::scientific << Utils::norm(seaIce->interface(x, SEAICE_TT_))
                  << ", ||G|| =" << std::setw(8) << std::scientific << Utils::norm(seaIce->interface(x, SEAICE_GG_))
                  << std::endl;

        if (Utils::norm(x) < 1e-8)
            break;
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
