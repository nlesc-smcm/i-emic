#include "TestDefinitions.H"

// Testing the serial and parallel atmosphere

//------------------------------------------------------------------
namespace
{
	std::shared_ptr<Atmosphere>    atmos;
	std::shared_ptr<AtmospherePar> atmosPar;
}


//------------------------------------------------------------------
TEST(ALL, Initialization)
{
	bool failed = false;
	// Create atmosphere parameters
	RCP<Teuchos::ParameterList> atmosphereParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());

	try
	{
		atmos = std::make_shared<Atmosphere>(atmosphereParams);
	}
	catch (...)
	{
		failed = true;
	}
	EXPECT_EQ(failed, false);

	failed = false;
	try
	{
		atmosPar = std::make_shared<AtmospherePar>(comm, atmosphereParams);
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
	atmos     = std::shared_ptr<Atmosphere>();
	atmosPar  = std::shared_ptr<AtmospherePar>();
	
	comm->Barrier();
	std::cout << "TEST exit code proc #" << comm->MyPID()
			  << " " << out << std::endl;


	MPI_Finalize();
	return out;
}
