#include "TestDefinitions.H"
#include <limits>
//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    std::shared_ptr<Ocean>         ocean;
    std::shared_ptr<AtmospherePar> atmos;
    std::shared_ptr<CoupledModel>  coupledModel;
    RCP<Teuchos::ParameterList> oceanParams;
    RCP<Teuchos::ParameterList> atmosphereParams;
    RCP<Teuchos::ParameterList> coupledmodelParams;
    RCP<Teuchos::ParameterList> continuationParams;
    RCP<Teuchos::ParameterList> jdqzParams;
    RCP<Epetra_Comm>            comm;
}

//------------------------------------------------------------------
TEST(ParameterLists, Initialization)
{
    bool failed = false;
    try
    {
        // Create parameter object for Ocean
        oceanParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
        oceanParams->setName("Ocean parameters");

        // Create parameter object for Atmosphere
        atmosphereParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());
        atmosphereParams->setName("Atmosphere parameters");

        // Create parameter object for CoupledModel
        coupledmodelParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("coupledmodel_params.xml", coupledmodelParams.ptr());
        coupledmodelParams->setName("CoupledModel parameters");

        // Create parameter object for Continuation
        continuationParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params.xml", continuationParams.ptr());
        continuationParams->setName("Continuation parameters");

        // Create parameter object for JDQZ
        jdqzParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("jdqz_params.xml", jdqzParams.ptr());
        jdqzParams->setName("JDQZ parameters");

        INFO('\n' << "Overwriting:");
        // The Continuation and CoupledModel parameterlists overwrite settings
        Utils::overwriteParameters(oceanParams,        coupledmodelParams);
        Utils::overwriteParameters(atmosphereParams,   coupledmodelParams);

        Utils::overwriteParameters(oceanParams,        continuationParams);
        Utils::overwriteParameters(atmosphereParams,   continuationParams);
        Utils::overwriteParameters(coupledmodelParams, continuationParams);
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
        ocean = std::make_shared<Ocean>(comm, oceanParams);
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
        atmos = std::make_shared<AtmospherePar>(comm, atmosphereParams);
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
        coupledModel = std::make_shared<CoupledModel>(ocean,atmos,coupledmodelParams);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}


//------------------------------------------------------------------
TEST(CoupledModel, Newton)
{
    // One step in a 'natural continuation'
    // initialize state in model
    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    stateV->PutScalar(0.0);

    std::shared_ptr<Combined_MultiVec> solV =
        coupledModel->getSolution('V');

    solV->PutScalar(0.0);

    // set parameter
    coupledModel->setPar(0.005);

    // try to converge
    int maxit = 10;
    std::shared_ptr<Combined_MultiVec> b;
    int niter = 0;
    for (; niter != maxit; ++niter)
    {
        coupledModel->computeRHS();

        coupledModel->computeJacobian();

        b = coupledModel->getRHS('C');

        INFO(" ocean F  = " << Utils::norm(coupledModel->getRHS('V')->First()) );
        INFO(" atmos F  = " << Utils::norm(coupledModel->getRHS('V')->Second()) );

        CHECK_ZERO(b->Scale(-1.0));

        double normb = Utils::norm(b);

        coupledModel->solve(b);

        std::shared_ptr<Combined_MultiVec> x = coupledModel->getSolution('C');
        std::shared_ptr<Combined_MultiVec> y = coupledModel->getSolution('C');

        INFO(" ocean x  = " << Utils::norm(stateV->First()) );
        INFO(" atmos x  = " << Utils::norm(stateV->Second()) );
        INFO(" ocean dx = " << Utils::norm(x->First()) );
        INFO(" atmos dx = " << Utils::norm(x->Second()) );

        stateV->Update(1.0, *x, 1.0); // x = x + dx;

        coupledModel->applyMatrix(*x, *y);

        y->Update(1.0, *b, -1.0);
        y->Scale(1./normb);

        Utils::print(y, "residual");

        INFO(" ocean ||r|| / ||b||  = " << Utils::norm(y->First()));
        INFO(" atmos ||r|| / ||b||  = " << Utils::norm(y->Second()));
        INFO(" total ||r|| / ||b||  = " << Utils::norm(y));

        if (Utils::norm(coupledModel->getRHS('V')) < 0.1)
            break;
    }

    EXPECT_LT(Utils::norm(coupledModel->getRHS('V')), 0.1);
    EXPECT_LT(niter, 10);
    INFO("CoupledModel, Newton converged in " << niter << " iterations");
}

//------------------------------------------------------------------
// 1st integral condition test for atmosphere
TEST(CoupledModel, AtmosphereIntegralCondition1)
{
    Teuchos::RCP<Epetra_Vector> intCoeff = atmos->getIntCoeff();
    Teuchos::RCP<Epetra_Vector> atmosX   = atmos->getState('C');

    double result = Utils::dot(intCoeff, atmosX);

    INFO("  atmosphere state norm: " << Utils::norm(atmosX));
    INFO("  atmosphere integral condition on q: " << result);
    
    EXPECT_NEAR(result, 0.0, 1e-4);
}

//------------------------------------------------------------------
TEST(CoupledModel, AtmosphereEPfields)
{
    bool failed = false;
    try
    {
        atmos->computeEP();
        Teuchos::RCP<Epetra_Vector> E = atmos->getE();
        Teuchos::RCP<Epetra_Vector> P = atmos->getP();

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
    
    Teuchos::RCP<Epetra_Vector> intcoeff = atmos->getPrecipIntCo();
    
    Teuchos::RCP<Epetra_Vector> E = atmos->getE();
    Teuchos::RCP<Epetra_Vector> P = atmos->getP();
    
    double integralE = Utils::dot(intcoeff, E);
    EXPECT_GT(std::abs(integralE), 0.0);
                              
    double integralP = Utils::dot(intcoeff, P);
    EXPECT_GT(std::abs(integralP), 0.0);

    double totalArea;
    intcoeff->Norm1(&totalArea);

    std::cout << "integralP = " << integralP << std::endl;
    std::cout << "integralE = " << integralE << std::endl;
    std::cout << "totalArea = " << totalArea << std::endl;

    EXPECT_NEAR(integralP, integralE, 1e-7);
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
        std::shared_ptr<Combined_MultiVec> stateV =
            coupledModel->getState('V');
        
        stateV->PutScalar(1.234);
        
        std::shared_ptr<Combined_MultiVec> solV =
            coupledModel->getSolution('V');
        
        solV->PutScalar(0.0);
        
        // set initial parameter
        coupledModel->setPar(0.0);
        
        // Create continuation
        Continuation<std::shared_ptr<CoupledModel>,
                     Teuchos::RCP<Teuchos::ParameterList> >
            continuation(coupledModel, continuationParams);

        // Test continuation
        continuation.test();

        stateV->PutScalar(0.0);
        solV->PutScalar(0.0);
         
        // Run continuation        
        continuation.run();

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
    Teuchos::RCP<Epetra_Vector> intcoeff = atmos->getPrecipIntCo();
    
    Teuchos::RCP<Epetra_Vector> E = atmos->getE();
    Teuchos::RCP<Epetra_Vector> P = atmos->getP();
    
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

        jdqz.setParameters(*jdqzParams);
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

//------------------------------------------------------------------
TEST(CoupledModel, numericalJacobian)
{
    // only do this test for small problems in serial
    int nmax = 2e3;

    // get Jacobian blocks in the model
    coupledModel->dumpBlocks();
    
    if ( (comm->NumProc() == 1) &&
         (coupledModel->getState('V')->GlobalLength() < nmax) )
    {
        bool failed = false;
        try
        {
            INFO("compute njC");

            NumericalJacobian<std::shared_ptr<CoupledModel>,
                              std::shared_ptr<Combined_MultiVec> > njC;

            njC.setTolerance(1e-11);
            njC.seth(1e-12);
            njC.compute(coupledModel, coupledModel->getState('V'));

            std::string fnameJnC("JnC");

            INFO(" Printing Numerical Jacobian " << fnameJnC);

            njC.print(fnameJnC);

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
        INFO("****Numerical Jacobian test cannot run in parallel****");
    }

    if (coupledModel->getState('V')->GlobalLength() > nmax)
    {
        INFO("****Numerical Jacobian test cannot run for this problem size****");
    }
}


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
    // TESTINGv
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    // Get rid of possibly parallel objects for a clean ending.
    ocean        = std::shared_ptr<Ocean>();
    atmos        = std::shared_ptr<AtmospherePar>();
    coupledModel = std::shared_ptr<CoupledModel>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;
    
    MPI_Finalize();
    return out;
}
