#include "TestDefinitions.H"

//------------------------------------------------------------------

namespace // local unnamed namespace (similar to static in C)
{	
	RCP<Ocean>                    ocean;
	std::shared_ptr<Atmosphere>   atmos;
	std::shared_ptr<CoupledModel> coupledModel;
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
		RCP<Teuchos::ParameterList> atmosphereParams =
			rcp(new Teuchos::ParameterList);
		updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());
		atmos = std::make_shared<Atmosphere>(atmosphereParams);
	}
	catch (...)
	{
		failed = true;
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
		RCP<Teuchos::ParameterList> coupledModelParams =
			rcp(new Teuchos::ParameterList);
		updateParametersFromXmlFile("coupledmodel_params.xml", coupledModelParams.ptr());
		coupledModel = std::make_shared<CoupledModel>(ocean,atmos,coupledModelParams);
	}
	catch (...)
	{
		failed = true;
	}
	
	EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Continuation)
{
	bool failed = false;
	try
	{
		// Create continuation
		RCP<Teuchos::ParameterList> continuationParams = rcp(new Teuchos::ParameterList);
		updateParametersFromXmlFile("continuation_params.xml", continuationParams.ptr());
		
		Continuation<std::shared_ptr<CoupledModel>, RCP<Teuchos::ParameterList> >
			continuation(coupledModel, continuationParams);

		continuation.run();
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
	ocean = Teuchos::null;
	
	comm->Barrier();
	std::cout << "TEST exit code proc #" << comm->MyPID()
			  << " " << out << std::endl;
	
	MPI_Finalize();
	return out;
}
