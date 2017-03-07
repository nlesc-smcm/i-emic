#include "TestDefinitions.H"

// Testing the serial and parallel atmosphere

//------------------------------------------------------------------
namespace
{
	std::shared_ptr<Atmosphere>    atmos;
	std::shared_ptr<AtmospherePar> atmosPar;
}

//------------------------------------------------------------------
TEST(Atmosphere, Initialization)
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
	catch (std::exception const &e)
	{
		INFO("TEST(Atmosphere, Initialization) exception: " << e.what());
		failed = true;
	}
	catch (int error)
	{
		INFO("TEST(Atmosphere, Initialization) exception: error code = " << error);
		failed = true;
	}
	catch (...)
	{
		INFO("TEST(Atmosphere, Initialization) some exception thrown...");
		failed = true;
	}
	EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, RHS)
{
	bool failed = false;
	try
	{
		// Check dimensions
		int dim1 = atmos->dim();
		int dim2 = atmosPar->dim();
		EXPECT_EQ(dim1, dim2);

		// Create random state from atmosPar standardMap state
		Teuchos::RCP<Epetra_Vector> state_par = atmosPar->getState('C');

		state_par->SetSeed(rand() % 100);		
		INFO("TEST(Atmosphere, RHS): seed = " << state_par->Seed());

		// Randomize and scale standardMap state
		state_par->Random();
		
		// Gather state_par into std::vector
		std::shared_ptr<std::vector<double> > state =
			getGatheredVector(state_par);
		
		// Put the full state in the serial atmosphere
		atmos->setState(state);

		// Put the standardMap parallel state into the parallel atmosphere
		atmosPar->setState(state_par);
		
		// Compute RHS in both models
		atmos->computeRHS();
		atmosPar->computeRHS();

		// Obtain RHS from both models
		Teuchos::RCP<Epetra_Vector> rhsPar = atmosPar->getRHS('V');
		std::shared_ptr<std::vector<double> > rhsSer = atmos->getRHS('V');

		// Check norms parallel and serial rhs
		EXPECT_NEAR(norm(rhsPar), norm(rhsSer), 1e-7);
		INFO("TEST(Atmosphere, RHS): serial norm = " << norm(rhsSer));
		INFO("TEST(Atmosphere, RHS): parall norm = " << norm(rhsPar));
	}
	catch (...)
	{
		failed = true;		
	}
	EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, State)
{
	bool failed = false;
	try
	{
		// Obtain state from both models
		Teuchos::RCP<Epetra_Vector> statePar = atmosPar->getState('V');
		std::shared_ptr<std::vector<double> > stateSer = atmos->getState('V');

		// Check norms parallel and serial rhs
		EXPECT_NEAR(norm(statePar), norm(stateSer), 1e-7);
	}
	catch(...)
	{
		failed = true;
	}
	EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, SurfaceTemperature)
{
	bool failed = false;
	try
	{
		// Obtain map for surface valeus
		Teuchos::RCP<Epetra_BlockMap> surfmap = atmosPar->getStandardSurfaceMap();

		// Create random vector from surface map
		Teuchos::RCP<Epetra_Vector> surfvals = Teuchos::rcp(new Epetra_Vector(*surfmap));
		surfvals->SetSeed(rand() % 100);		
		INFO("TEST(Atmosphere, SurfaceTemperature): seed = " << surfvals->Seed());

		surfvals->Random();
		surfvals->Scale(1e1);

		// Create gather map
		Teuchos::RCP<Epetra_BlockMap> gmap =
			Utils::AllGather(*surfmap);

		// Create import strategy
		Teuchos::RCP<Epetra_Import> gimp =
			Teuchos::rcp(new Epetra_Import(*gmap, *surfmap));

		// Create gather vector
		Teuchos::RCP<Epetra_Vector> gathered = Teuchos::rcp(new Epetra_Vector(*gmap));

		// Gather <surfvals> into <gathered>
		gathered->Import(*surfvals, *gimp, Insert);

		// Create vector for serial atmosphere
		int size = surfmap->NumGlobalElements();
		std::shared_ptr<std::vector<double> > fullsurfvals =
			std::make_shared<std::vector<double> >(size);

		// Extract gathered 
		gathered->ExtractCopy(&(*fullsurfvals)[0], size);

		// Put the randomized ocean surface values in both models
		atmos->setOceanTemperature(*fullsurfvals);
		atmosPar->setOceanTemperature(surfvals);

		// Compute RHS in both models
		atmos->computeRHS();
		atmosPar->computeRHS();

		// Obtain RHS from both models
		Teuchos::RCP<Epetra_Vector> rhsPar = atmosPar->getRHS('V');
		std::shared_ptr<std::vector<double> > rhsSer = atmos->getRHS('V');

		// Check norms parallel and serial rhs
		EXPECT_NEAR(norm(rhsPar), norm(rhsSer), 1e-7);
		INFO("TEST(Atmosphere, SurfaceTemperature): serial norm = " << norm(rhsSer));
		INFO("TEST(Atmosphere, SurfaceTemperature): parall norm = " << norm(rhsPar));

	}
	catch (...)
	{
		failed = true;
	}
	EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, Jacobian)
{
	bool failed = false;
	try
	{
		atmos->computeJacobian();
		atmosPar->computeJacobian();

		// Obtain vector and randomize
		Teuchos::RCP<Epetra_Vector> x = atmosPar->getState('C');
		x->Random();
		x->Scale(1e1);

		// Gather randomized vector
		std::shared_ptr<std::vector<double> > gx = getGatheredVector(x);

		// Check the gather routine
		EXPECT_NEAR(norm(gx), norm(x), 1e-7);

		std::shared_ptr<std::vector<double> > outSer =
			std::make_shared<std::vector<double> >(gx->size(), 0.0);

		// apply serial Jacobian
		atmos->applyMatrix(gx, outSer);

		Teuchos::RCP<Epetra_Vector> outPar = atmosPar->getState('C');
		outPar->PutScalar(0.0);

		// apply parallel Jacobian
		atmosPar->applyMatrix(*x, *outPar);

		// check results
		EXPECT_NEAR(norm(outSer), norm(outPar), 1e-7);

		INFO("TEST(Atmosphere, apply Jacobian): ||serial out|| = " << norm(outSer));
		INFO("TEST(Atmosphere, apply Jacobian): ||parall out|| = " << norm(outPar));
	}
	catch (std::exception const &e)
	{
		INFO("TEST(Atmosphere, Initialization) exception: " << e.what());
		failed = true;
	}
	catch (int error)
	{
		INFO("TEST(Atmosphere, Initialization) exception: error code = " << error);
		failed = true;
	}
	catch (...)
	{
		INFO("TEST(Atmosphere, Initialization) some exception thrown...");
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
