#include "TestDefinitions.H"

#include "NumericalJacobian.H"
#include "Ocean.H"
#include "Atmosphere.H"
#include "SeaIce.H"
#include "CoupledModel.H"
#include "Continuation.H"

#ifdef HAVE_JDQZPP
#include "ComplexVector.H"
#include "JDQZInterface.H"
#include "jdqz.hpp"
#endif

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Epetra_Comm>               comm;
    std::shared_ptr<Ocean>         ocean;
    std::shared_ptr<Atmosphere>    atmos;
    std::shared_ptr<SeaIce>        seaice;
    std::shared_ptr<CoupledModel>  coupledModel;
    std::vector<Teuchos::RCP<Teuchos::ParameterList> > params;
    enum Ident { OCEAN, ATMOS, SEAICE, COUPLED, CONT, EIGEN};
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
                                          "continuation_params.xml",
                                          "jdqz_params.xml"};

        std::vector<std::string> names = {"Ocean parameters",
                                          "Atmosphere parameters",
                                          "Sea ice parameters",
                                          "CoupledModel parameters",
                                          "Continuation parameters",
                                          "JDQZ parameters"};

        for (int i = 0; i != (int) files.size(); ++i)
        {
            params.push_back(Utils::obtainParams(files[i], names[i]));
        }

        params[OCEAN]->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");

        INFO('\n' << "Overwriting:");
        // Allow dominant parameterlists. Not that this trick ignores
        // any hierarchy. The Continuation and CoupledModel
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
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, RHS)
{
    coupledModel->computeRHS();
    std::shared_ptr<Combined_MultiVec> F = coupledModel->getRHS('V');
    for (int i = 0; i != F->Size(); ++i)
    {
        INFO(" submodel " << i << ": ||F|| = " << Utils::norm((*F)(i)));
    }

    EXPECT_LT(Utils::norm(F), 1e-7);
}

//------------------------------------------------------------------
TEST(CoupledModel, Newton)
{
    // Initialize the model state
    coupledModel->initializeState();

    // One step in a 'natural continuation'
    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    std::shared_ptr<Combined_MultiVec> solV =
        coupledModel->getSolution('V');

    solV->PutScalar(0.0);

    // set parameter
    coupledModel->setPar("Combined Forcing", 0.01);

    // try to converge
    int maxit = 10;
    std::shared_ptr<Combined_MultiVec> F = coupledModel->getRHS('V');
    std::shared_ptr<Combined_MultiVec> b = coupledModel->getRHS('C');
    int niter = 0;
    std::shared_ptr<Combined_MultiVec> x = coupledModel->getSolution('V');
    std::shared_ptr<Combined_MultiVec> y = coupledModel->getSolution('C');
    coupledModel->computeJacobian();
    coupledModel->computeRHS();

    for (; niter < maxit; ++niter)
    {
        coupledModel->computeJacobian();
        coupledModel->computeRHS();

        b = coupledModel->getRHS('C');
        CHECK_ZERO(b->Scale(-1.0));

        coupledModel->solve(b);  // J dx = - F

        stateV->Update(1.0, *x, 1.0); // x = x + dx;

        coupledModel->computeRHS();
        coupledModel->computeJacobian();

        coupledModel->applyMatrix(*x, *y);
        double normb = Utils::norm(b);
        y->Update(1.0, *b, -1.0);
        y->Scale(1. / normb);

        INFO("\n ||r|| / ||b|| = " << Utils::norm(y));
        INFO("        ||dx|| = " << Utils::norm(x) << " ");
        INFO("         ||F|| = " << Utils::norm(F) << "\n");

        for (int i = 0; i != b->Size(); ++i)
        {
            INFO(" submodel " << i << " ||F|| = " << Utils::norm((*b)(i)));
            INFO("   " << " dx = " << Utils::norm((*x)(i)));
            INFO("   " << " ||r|| / ||b|| = " << Utils::norm((*y)(i)));
        }
        if (Utils::norm(coupledModel->getRHS('V')) < 1e-8)
            break;
    }
    Utils::save(F, "F");
    EXPECT_LT(Utils::norm(coupledModel->getRHS('V')), 1e-8);
    EXPECT_LT(niter, maxit);
}

//------------------------------------------------------------------
TEST(CoupledModel, numericalJacobian)
{
    // only do this test for small problems in serial
    int nmax = 2e3;

    if ( (comm->NumProc() == 1) &&
         (coupledModel->getState('V')->GlobalLength() < nmax) )
    {
        bool failed = false;
        try
        {
            // get analytical jacobian blocks
            coupledModel->computeJacobian();
            coupledModel->dumpBlocks();

            INFO("compute njC");

            NumericalJacobian<std::shared_ptr<CoupledModel>,
                              std::shared_ptr<Combined_MultiVec> > njC;

            njC.setTolerance(1e-10);
            njC.seth(1e-6);
            njC.compute( coupledModel, coupledModel->getState('V') );

            std::string fnameJnC("JnC");

            INFO(" Printing Numerical Jacobian " << fnameJnC);

            njC.print(fnameJnC);

            // test individual elements
            NumericalJacobian<std::shared_ptr<CoupledModel>,
                              std::shared_ptr<Combined_MultiVec> >::CCS ccs;
            njC.fillCCS(ccs);

            EXPECT_NE(ccs.beg.back(), 0);

            std::shared_ptr<Combined_MultiVec> x = coupledModel->getState('C');

            testEntries(coupledModel, ccs, x);

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

    if (coupledModel->getState('V')->GlobalLength() > nmax)
    {
        std::cout << ("****Numerical Jacobian test cannot run for this problem size****\n");
        INFO("****Numerical Jacobian test cannot run for this problem size****");
    }
}


//-----------------------------------------------------------------
// 1st integral condition test for atmosphere
TEST(CoupledModel, AtmosphereIntegralCondition1)
{
    Teuchos::RCP<Epetra_Vector> intCoeff = atmos->getIntCoeff();
    Teuchos::RCP<Epetra_Vector> atmosX   = atmos->getState('C');

    double result = Utils::dot(intCoeff, atmosX);

    EXPECT_GT(Utils::norm(atmosX), 1e-7);
    INFO("  atmosphere state norm: " << Utils::norm(atmosX));
    INFO("  atmosphere integral condition on q: " << result);

    EXPECT_NEAR(result, 0.0, 1e-7);
}


//------------------------------------------------------------------
TEST(CoupledModel, AtmosphereEPfields)
{
    bool failed = false;
    try
    {
        Teuchos::RCP<Epetra_Vector> E = atmos->interfaceE();
        Teuchos::RCP<Epetra_Vector> P = atmos->interfaceP();

        Utils::print(E, "file.txt");
        Utils::print(P, "file.txt");

        double tol = 1e-12;

        int nnzE = Utils::nnz(E, tol);
        int nnzP = Utils::nnz(P, tol);
        EXPECT_EQ(nnzE, nnzP);

        int nelE = E->Map().NumMyElements();
        int nelP = P->Map().NumMyElements();

        EXPECT_EQ(nelE, nelP);

        // expecting data distributions to be equal
        for (int i = 0; i != nelE; ++i)
        {
            if ( std::abs((*E)[i]) > tol )
            {
                EXPECT_GT( std::abs((*P)[i]), tol );
            }
        }
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, EPIntegral)
{
    Teuchos::RCP<Epetra_Vector> intcoeff = atmos->getPIntCoeff();

    Teuchos::RCP<Epetra_Vector> E = atmos->interfaceE();
    Teuchos::RCP<Epetra_Vector> P = atmos->interfaceP();

    double integralE = Utils::dot(intcoeff, E);
    EXPECT_GT(std::abs(integralE), 0.0);

    double integralP = Utils::dot(intcoeff, P);
    EXPECT_GT(std::abs(integralP), 0.0);

    double totalArea;
    intcoeff->Norm1(&totalArea);

    std::cout << "integralP = " << integralP << std::endl;
    std::cout << "integralE = " << integralE << std::endl;
    std::cout << "totalArea = " << totalArea << std::endl;

    // not sure how strict this test should be
    EXPECT_NEAR( (integralP - integralE) / integralP, 0.0, 1e-7);
}

//------------------------------------------------------------------
// full continuation
TEST(CoupledModel, Continuation)
{
    bool failed = false;
    try
    {
        // One step in an arclength continuation
        // initialize state in model
        coupledModel->setPar("Combined Forcing", 0.0);
        coupledModel->initializeState();

        std::shared_ptr<Combined_MultiVec> solV =
            coupledModel->getSolution('V');

        solV->PutScalar(0.0);

        // Create continuation
        Continuation<std::shared_ptr<CoupledModel>>
            continuation(coupledModel, params[CONT]);

        // Run continuation
        int status = continuation.run();
        EXPECT_EQ(status, 0);

        // Dump blocks
        coupledModel->dumpBlocks();
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}


//------------------------------------------------------------------
TEST(CoupledModel, EPIntegral2)
{
    Teuchos::RCP<Epetra_Vector> intcoeff = atmos->getPIntCoeff();

    Teuchos::RCP<Epetra_Vector> E = atmos->interfaceE();
    Teuchos::RCP<Epetra_Vector> P = atmos->interfaceP();

    double integralE = Utils::dot(intcoeff, E);

    double integralP = Utils::dot(intcoeff, P);
    EXPECT_NEAR(integralE, integralP, 1e-8);

    // Compute integral of E-P
    E->Update(-1.0, *P, 1.0);
    double integralEP = Utils::dot(intcoeff, E);

    EXPECT_NEAR(std::abs(integralEP), 0.0, 1e-8);

    INFO(" int P   " << integralP);
    INFO(" int E   " << integralE);
    INFO(" int E-P " << integralEP);
}

#ifdef HAVE_JDQZPP

//------------------------------------------------------------------
TEST(CoupledModel, JDQZSolve)
{
    bool failed = false;
    try
    {
        INFO("Creating ComplexVector...");
        std::shared_ptr<Combined_MultiVec> x = coupledModel->getSolution('C');
        std::shared_ptr<Combined_MultiVec> y = coupledModel->getSolution('C');

        x->PutScalar(0.0);
        y->PutScalar(0.0);

        ComplexVector<Combined_MultiVec> z(*x, *y);
        ComplexVector<Combined_MultiVec> r(*x, *y); // residue
        ComplexVector<Combined_MultiVec> t(*x, *y); // tmp

        JDQZInterface<std::shared_ptr<CoupledModel>,
                      ComplexVector<Combined_MultiVec> > matrix(coupledModel, z);

        JDQZ<JDQZInterface<std::shared_ptr<CoupledModel>,
                           ComplexVector<Combined_MultiVec> > > jdqz(matrix, z);

        jdqz.setParameters(*params[EIGEN]);
        jdqz.printParameters();

        jdqz.solve();

        std::vector<ComplexVector<Combined_MultiVec> > eivec = jdqz.getEigenVectors();
        std::vector<std::complex<double> >             alpha = jdqz.getAlpha();
        std::vector<std::complex<double> >             beta  = jdqz.getBeta();

        EXPECT_GT(jdqz.kmax(), 0);
        for (int j = 0; j != jdqz.kmax(); ++j)
        {
            r.zero();
            t.zero();
            matrix.AMUL(eivec[j], r);
            r.scale(beta[j]);
            matrix.BMUL(eivec[j], t);
            r.axpy(-alpha[j], t);
            std::cout << "alpha: " << std::setw(30) << alpha[j]
                      << " beta: " << std::setw(15) << beta[j];
            std::cout << " alpha) / beta: " << std::setw(30)
                      << alpha[j] / beta[j];
            std::cout << " residue: " << r.norm() << std::endl;
            EXPECT_LT(r.norm(), 1e-7);
        }
    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

#endif

//------------------------------------------------------------------
// 2nd integral condition test
TEST(CoupledModel, AtmosphereIntegralCondition2)
{
    Teuchos::RCP<Epetra_Vector> intCoeff = atmos->getIntCoeff();
    Teuchos::RCP<Epetra_Vector> atmosX   = atmos->getState('C');

    double result = Utils::dot(intCoeff, atmosX);

    INFO("  atmosphere state norm: " << Utils::norm(atmosX));
    INFO("  atmosphere integral condition on q: " << result);

    EXPECT_NEAR(result, 0.0, 1e-7);
}

//------------------------------------------------------------------
// using the solution from the previous continuation we check the
// integrity of the Jacobian matrix
TEST(CoupledModel, SmallPerturbation)
{
    std::shared_ptr<Combined_MultiVec> x  = coupledModel->getState('V');
    std::shared_ptr<Combined_MultiVec> xp = coupledModel->getState('C');
    xp->Scale(0.01);                // perturbation
    double nrmxp = Utils::norm(xp);
    x->Update(1.0, *xp, 1.0);       // perturb state
    coupledModel->computeRHS();

    // temporary
    Combined_MultiVec tmp(*x);
    tmp.PutScalar(0.0);

    double nrmp = Utils::norm(coupledModel->getRHS('V')); // perturbed norm
    INFO("Perturbed norm: " << nrmp / nrmxp);

    x->Update(-1.0, *xp, 1.0);   // unperturb state
    coupledModel->computeRHS();
    coupledModel->computeJacobian();

    coupledModel->applyMatrix(*xp, tmp);
    tmp.Update(1.0,*coupledModel->getRHS('V'),1.0);

    double nrm = Utils::norm(&tmp); // perturbed norm
    INFO("Linearized norm: " << nrm / nrmxp);

    EXPECT_NEAR(nrm / nrmxp ,nrmp / nrmxp,1e-03);
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
