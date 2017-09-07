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
        coupledModel->computeJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> atmosJac = atmos->getJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> oceanJac = ocean->getJacobian();

        Utils::print(atmosJac, "atmosJac");
        Utils::print(oceanJac, "oceanJac");

        Utils::print(&atmosJac->ColMap(), "atmosJacColMap");
        Utils::print(&atmosJac->DomainMap(), "atmosJacDomainMap");

    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
// We need this information from THCM
extern "C" _SUBROUTINE_(getdeps)(double*, double*, double*);

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

        double value[3] = {1.234, 2.12, -23.5};

        // Test atmos -> ocean coupling
        for (int v = 0; v != 3; ++v)
        {
            atmosVec->PutScalar(value[v]);
            oceanVec->PutScalar(0.0);

            C12.applyMatrix(*atmosVec, *oceanVec);

            // Get ocean parameters
            double Ooa, Os, qdep;
            FNAME(getdeps)(&Ooa, &Os, &qdep);

            double max = 0;
            double el;
            for (int i = 0; i != oceanVec->MyLength(); ++i)
            {
                el = std::abs( (*oceanVec)[0][i] );
                max = (el > max) ? el : max;
            }

            EXPECT_NEAR( std::abs(-Ooa * value[v]), max, 1e-7);
        }

        // Test ocean -> atmos coupling
        for (int v = 0; v != 3; ++v)
        {
            oceanVec->PutScalar(0.0);
            for (int i = 4; i < oceanVec->MyLength(); i+=6)
            {
                (*oceanVec)[0][i] = value[v];
            }

            atmosVec->PutScalar(0.0);

            C21.applyMatrix(*oceanVec, *atmosVec);

            double max = 0;
            double el;
            for (int i = 0; i != atmosVec->MyLength(); ++i)
            {
                el = std::abs( (*atmosVec)[0][i] );
                max = (el > max) ? el : max;
            }
            EXPECT_NEAR( std::abs(value[v]), max, 1e-7);

            oceanVec->PutScalar(0.0);
            for (int i = 3; i < oceanVec->MyLength(); i+=6)
            {
                (*oceanVec)[0][i] = value[v];
            }

            atmosVec->PutScalar(0.0);
            C21.applyMatrix(*oceanVec, *atmosVec);

            max = 0;
            for (int i = 0; i != atmosVec->MyLength(); ++i)
            {
                el = std::abs( (*atmosVec)[0][i] );
                max = (el > max) ? el : max;
            }
            EXPECT_NE( std::abs(value[v]), max);

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
        coupledModel->setPar(0.005);

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
int main(int argc, char **argv)
{
    // Initialize the environment:
    initializeEnvironment(argc, argv);
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
