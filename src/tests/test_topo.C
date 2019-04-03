#include "TestDefinitions.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
	RCP<Ocean>  ocean;
	RCP<Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> > > topo;
    RCP<Epetra_Comm>  comm;
}

//------------------------------------------------------------------
TEST(Ocean, Initialization)
{
	bool failed = false;
	try
	{
		// Create parallel Ocean
		RCP<Teuchos::ParameterList> oceanParams =
			rcp(new Teuchos::ParameterList);
		updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
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
		RCP<Teuchos::ParameterList> topoParams = rcp(new Teuchos::ParameterList);
		updateParametersFromXmlFile("topo_params.xml", topoParams.ptr());
		topo = rcp(new Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> >
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
TEST(Topo, RHS)
{
    bool failed = false;
    try
    {
        topo->setPar(0.3);
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
TEST(Topo, SpinupContinuation)
{
	bool failed = false;
	try
	{
        // not sure if needed
        ocean->setPar(0.0);
        ocean->getState('V')->PutScalar(0.0);
        topo->getSolution('V')->PutScalar(0.0); // superfluous

		// Create parameter object for continuation
		RCP<Teuchos::ParameterList> continuationParams =
			rcp(new Teuchos::ParameterList);
		updateParametersFromXmlFile("spinup_continuation_params.xml",
									continuationParams.ptr());

		// Create spinup continuation
        Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
            continuation(ocean, continuationParams);

		// Run continuation
        INFO("--**-- Topo test: running spinup...");
		int status = continuation.run();
        EXPECT_EQ(status, 0);
        INFO("--**-- Topo test: running spinup... done");
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
	try
	{
        // not sure if needed
        topo->setPar(0.0);

		// Create parameter object for continuation
		RCP<Teuchos::ParameterList> continuationParams =
			rcp(new Teuchos::ParameterList);
		updateParametersFromXmlFile("topo_continuation_params.xml",
									continuationParams.ptr());

		// Create topo continuation
		Continuation<RCP<Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> > >,
					 RCP<Teuchos::ParameterList> >
			continuation(topo, continuationParams);

        // Run topo continuation
        int nMasks    = topo->nMasks();
        int startMask = topo->startMaskIdx();

        INFO(" Running topo cont...");// Run continuation

        for (int maskIdx = startMask; maskIdx != nMasks-1; maskIdx++)
        {
            topo->setMaskIndex(maskIdx);
            topo->initialize();

            topo->predictor();

            int status = continuation.run();
            EXPECT_EQ(status, 0);

        }
        INFO(" Running topo cont... done");


	}
	catch (...)
	{
		failed = true;
        throw;
	}

	EXPECT_EQ(failed, false);
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
