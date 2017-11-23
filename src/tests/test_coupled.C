#include "TestDefinitions.H"
#include "NumericalJacobian.H"

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
TEST(CoupledModel, inspectState)
{
    bool failed = false;
    try
    {
        std::shared_ptr<Combined_MultiVec> state = coupledModel->getState('V');
        state->Random();

        int firstL  = state->First()->GlobalLength();
        int secondL = state->Second()->GlobalLength();
        int stateL  = state->GlobalLength();

        INFO(" global 1: " << firstL << " 2: " << secondL
             << " 1+2: " << stateL);

        EXPECT_EQ(firstL + secondL, stateL);

        firstL  = state->First()->MyLength();
        secondL = state->Second()->MyLength();
        stateL  = state->MyLength();

        INFO( " local 1: " << firstL << " 2: " << secondL
             << " 1+2: " << stateL );
        EXPECT_EQ(firstL + secondL, stateL);

        double firstNrm = Utils::norm(state->First());
        double secndNrm = Utils::norm(state->Second());
        double stateNrm = Utils::norm(state);

        EXPECT_NEAR(stateNrm, sqrt(pow(firstNrm,2) + pow(secndNrm,2)), 1e-7);

        INFO( " norm 1: " << firstNrm << " 2: " << secndNrm
              << " 1+2: " << stateNrm << std::endl );

        Utils::print(state->First(),  "rand_state_first");
        Utils::print(state->Second(), "rand_state_second");
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, computeJacobian)
{
    bool failed = false;
    try
    {
        // Set a small parameter
        coupledModel->setPar(0.1);

        // randomize state
        coupledModel->getState('V')->Random();

        coupledModel->computeJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> atmosJac = atmos->getJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> oceanJac = ocean->getJacobian();

        Utils::print(atmosJac, "atmosJac");
        Utils::print(oceanJac, "oceanJac");

        Utils::print(&atmosJac->ColMap(), "atmosJacColMap");
        Utils::print(&atmosJac->DomainMap(), "atmosJacDomainMap");

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
            INFO("compute njC");

            NumericalJacobian<std::shared_ptr<CoupledModel>,
                              std::shared_ptr<Combined_MultiVec> > njC;

            njC.setTolerance(1e-10);
            njC.seth(1);
            njC.compute(coupledModel, coupledModel->getState('V'));

            std::string fnameJnC("JnC");

            INFO(" Printing Numerical Jacobian " << fnameJnC);

            njC.print(fnameJnC);

            // --> todo test equality of element sums (see test_domain.C)
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
// We need this information from THCM
extern "C" _SUBROUTINE_(getdeps)(double*, double*, double*, double*);

TEST(CoupledModel, applyMatrix)
{
    bool failed = false;

    try
    {
        std::shared_ptr<Combined_MultiVec> x = coupledModel->getState('C');
        std::shared_ptr<Combined_MultiVec> y = coupledModel->getState('C');

        double normIn = Utils::norm(y);

        coupledModel->applyMatrix(*x, *y);

        double normOut = Utils::norm(y);

        EXPECT_NE(normIn, normOut);

        // Test the coupling blocks in the matrix separately

        CouplingBlock<std::shared_ptr<Ocean>,
                      std::shared_ptr<AtmospherePar> > C12(ocean, atmos);

        CouplingBlock<std::shared_ptr<AtmospherePar>,
                      std::shared_ptr<Ocean> > C21(atmos, ocean);

        Teuchos::RCP<Epetra_MultiVector> oceanVec = x->First();
        Teuchos::RCP<Epetra_MultiVector> atmosVec = x->Second();

        int n = ocean->getNdim();
        int m = ocean->getMdim();
        int l = ocean->getLdim();
        double value[3] = {1.234, 2.12, -23.5};

        // Test atmos -> ocean coupling
        for (int v = 0; v != 3; ++v)
        {
            atmosVec->PutScalar(value[v]);
            oceanVec->PutScalar(0.0);

            C12.applyMatrix(*atmosVec, *oceanVec);

            // Get ocean parameters
            double Ooa, Os, gamma, eta;
            FNAME(getdeps)(&Ooa, &Os, &gamma, &eta);

            // Test first surface element (temperature)
            int surfbT = FIND_ROW2(_NUN_, n, m, l, 0, 0, l-1, TT);

            INFO( "first surface element TT " << surfbT );

            int lid;
            double surfval;
            if (oceanVec->Map().MyGID(surfbT))
            {
                lid = oceanVec->Map().LID(surfbT);
                surfval = (*oceanVec)[0][lid];
                EXPECT_NEAR(-Ooa * value[v], surfval , 1e-7);
            }

            // Test first surface element (salinity)
            int surfbS = FIND_ROW2(_NUN_, n, m, l, 0, 0, l-1, SS);
            INFO( "first surface element SS " << surfbS );

            if (oceanVec->Map().MyGID(surfbS))
            {
                lid = oceanVec->Map().LID(surfbS);
                surfval = (*oceanVec)[0][lid];
                EXPECT_NEAR( (-eta * gamma - gamma ) *value[v], surfval , 1e-7);
            }
        }

        double qdim, nuq, eta, dqso, dqdt;
        atmos->getConstants(qdim, nuq, eta, dqso, dqdt);

        // Test ocean -> atmos coupling
        for (int v = 0; v != 3; ++v)
        {
            // Put values in temperature domain points. Atmosphere
            // only depends on sst.
            oceanVec->PutScalar(0.0);
            for (int i = TT-1; i < oceanVec->MyLength(); i += _NUN_)
            {
                (*oceanVec)[0][i] = value[v];
            }

            atmosVec->PutScalar(0.0);

            C21.applyMatrix(*oceanVec, *atmosVec);

            // assume single atmosphere layer
            l = 1;

            // test atmos temperature point
            int surfbT = FIND_ROW_ATMOS0(ATMOS_NUN_, n, m, l, 0, 0, l-1, ATMOS_TT_);

            int lid;
            double surfval;
            if (atmosVec->Map().MyGID(surfbT))
            {
                lid = atmosVec->Map().LID(surfbT);
                surfval = (*atmosVec)[0][lid];
                EXPECT_NEAR( value[v], surfval, 1e-7);
            }

            // Test atmos humidity point
            int surfbQ = FIND_ROW_ATMOS0(ATMOS_NUN_, n, m, l, 0, 0, l-1, ATMOS_QQ_);

            if (atmosVec->Map().MyGID(surfbQ))
            {
                lid = atmosVec->Map().LID(surfbQ);
                surfval = (*atmosVec)[0][lid];
                EXPECT_NEAR( dqdt * value[v], surfval, 1e-7);
            }

            // Test final element in range (results from sst integral)

            // First check if we have auxiliary unknowns:
            if (atmosVec->GlobalLength() > ATMOS_NUN_ * m * n * l)
            {
                Teuchos::RCP<Epetra_Vector> precipintco = atmos->getPrecipIntCo();
                double totalArea;
                precipintco->Norm1(&totalArea);

                Teuchos::RCP<Epetra_Vector> ones =
                    Teuchos::rcp(new Epetra_Vector(*precipintco));
                ones->PutScalar(1.0);

                double sstInt = Utils::dot(precipintco, ones);
                double intval = sstInt * (eta / totalArea) * (dqso / qdim);

                int last = FIND_ROW_ATMOS0(ATMOS_NUN_, n, m, l, n-1 , m-1, l-1, ATMOS_QQ_);
                if (atmosVec->Map().MyGID(last + 1))
                {
                    lid = atmosVec->Map().LID(last + 1);
                    surfval = (*atmosVec)[0][lid];
                    EXPECT_NEAR( intval * value[v], surfval, 1e-7);
                }
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
TEST(CoupledModel, View)
{
    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    std::shared_ptr<Combined_MultiVec> rhsV =
        coupledModel->getRHS('V');

    stateV->PutScalar(0.0);
    coupledModel->computeRHS();

    double norm1 = Utils::norm(rhsV);

    stateV->PutScalar(1.0);
    coupledModel->computeRHS();

    double norm2 = Utils::norm(rhsV);

    EXPECT_NE(norm1, norm2);
}

//------------------------------------------------------------------
TEST(CoupledModel, Synchronization)
{
    bool failed = false;

    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    try
    {
        stateV->First()->PutScalar(1.234);
        stateV->Second()->PutScalar(2.345);

        // Set a small parameter
        coupledModel->setPar(0.01);

        // At RHS computation the coupledModel synchronizes the states
        coupledModel->computeRHS();
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    // Obtain atmosphere temperature existing in ocean
    Teuchos::RCP<Epetra_Vector> oceanAtmosT  = ocean->getLocalAtmosT();

    double maxValue, minValue, meanValue;
    oceanAtmosT->MaxValue(&maxValue);
    EXPECT_EQ(maxValue, 2.345);

    oceanAtmosT->MinValue(&minValue);
    EXPECT_EQ(minValue, 2.345);

    if ( oceanAtmosT->Map().UniqueGIDs() )
    {
        oceanAtmosT->MeanValue(&meanValue);
        EXPECT_NEAR(meanValue, 2.345, 1e-7);
    }
    else
        INFO(" oceanAtmosT (overl.) GID's not unique, which is fine in parallel.");

    Teuchos::RCP<Epetra_Vector> atmosOceanT = atmos->getLocalOceanT();

    Utils::print(atmosOceanT, "atmosOceanT" + std::to_string(comm->MyPID()));

    atmosOceanT->MaxValue(&maxValue);
    EXPECT_EQ(maxValue, 1.234);
    atmosOceanT->MinValue(&minValue);
    EXPECT_EQ(minValue, 1.234);

    if (atmosOceanT->Map().UniqueGIDs())
    {
        atmosOceanT->MeanValue(&meanValue);
        EXPECT_NEAR(meanValue, 1.234, 1e-7);
    }
    else
        INFO(" atmosOceanT (overl.) GID's not unique, which is fine in parallel.");

    // Randomize combined state
    stateV->Random();

    try
    {
        // At RHS computation the coupledModel synchronizes the states
        coupledModel->computeRHS();
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    // Evaporation is calculated simultaneously in Ocean and in Atmosphere
    // during the RHS computation above.
    // Here we check whether they return the same vector.
    Teuchos::RCP<Epetra_Vector> atmosE = atmos->getE();
    Teuchos::RCP<Epetra_Vector> oceanE = ocean->getE();
    double nrmAtmosE = Utils::norm(atmosE);
    double nrmOceanE = Utils::norm(oceanE);
    EXPECT_NEAR(nrmAtmosE, nrmOceanE, 1e-7);

    // Precipitation should have been calculated in the atmosphere model
    Teuchos::RCP<Epetra_Vector> oceanAtmosP = ocean->getLocalAtmosP();

#ifdef GNU
    Utils::print(oceanAtmosT,  "oceanAtmosT" +
                 std::to_string(oceanAtmosT->Map().Comm().MyPID()) + ".txt");
    Utils::print(oceanAtmosP,  "oceanAtmosP" +
                 std::to_string(oceanAtmosP->Map().Comm().MyPID()) + ".txt");
#endif

    oceanAtmosP->MaxValue(&maxValue);
    EXPECT_GT(std::abs(maxValue), 0.0);
}


//------------------------------------------------------------------
// Test hashing functions of Combined_MultiVec and Utils
TEST(CoupledModel, Hashing)
{
    std::shared_ptr<Combined_MultiVec> x  = coupledModel->getState('C');
    std::size_t hash1 = x->hash(); INFO(" hash1 = " << hash1);
    std::size_t hash2 = x->hash(); INFO(" hash2 = " << hash2);
    EXPECT_EQ(hash1, hash2);

    x->Scale(1.0001);

    std::size_t hash3 = x->hash(); INFO(" hash3 = " << hash3);

    EXPECT_NE(hash2, hash3);

    x->PutScalar(1.00001);
    std::size_t hash4 = x->hash(); INFO(" hash4 = " << hash4);

    EXPECT_NE(hash3, hash4);

}

//------------------------------------------------------------------
TEST(CoupledModel, Solve)
{
    bool failed = false;
    try
    {
        std::shared_ptr<Combined_MultiVec> x = coupledModel->getState('C');
        std::shared_ptr<Combined_MultiVec> b = coupledModel->getState('C');
        x->PutScalar(1.0);
        coupledModel->applyMatrix(*x, *b);
        coupledModel->solve(b);
        std::shared_ptr<Combined_MultiVec> sol = coupledModel->getSolution('C');
    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    initializeEnvironment(argc, argv);
    if (outFile == Teuchos::null)
        throw std::runtime_error("ERROR: Specify output streams");

    ::testing::InitGoogleTest(&argc, argv);

    // -------------------------------------------------------
    // TESTING
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
