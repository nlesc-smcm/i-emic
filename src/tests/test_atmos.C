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
TEST(ALL, RHS)
{
	bool failed = false;
	try
	{
		// Check dimensions
		int dim1 = atmos->dim();
		int dim2 = atmosPar->dim();
		EXPECT_EQ(dim1, dim2);

		// Create random state from atmosPar standardMap state
		Teuchos::RCP<Epetra_Vector> state_par = atmosPar->getState();
		Teuchos::RCP<Epetra_BlockMap> gmap = Utils::AllGather(state_par->Map());
		Teuchos::RCP<Epetra_Import> imp =
			Teuchos::rcp(new Epetra_Import(*gmap, state_par->Map()));
		state_par->Random();
		Teuchos::RCP<Epetra_Vector> gathered =
			Teuchos::rcp(new Epetra_Vector(*gmap));

		gathered->Import(*state_par, *imp, Insert);

		std::shared_ptr< std::vector<double> > state =
			std::make_shared<std::vector<double> >(dim1, 0.0);
		
		gathered->ExtractCopy(&(*state)[0], dim1);

		atmos->setState(state);
		atmosPar->setState(state_par);
		
		// Compute RHS in both models
		atmos->computeRHS();
		atmosPar->computeRHS();
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
