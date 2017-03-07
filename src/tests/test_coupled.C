#include "TestDefinitions.H"

//------------------------------------------------------------------

namespace // local unnamed namespace (similar to static in C)
{	
	std::shared_ptr<Ocean>         ocean;
	std::shared_ptr<AtmospherePar> atmos;
	std::shared_ptr<CoupledModel>  coupledModel;
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
		
		ocean = std::make_shared<Ocean>(comm, oceanParams);	
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
		
		atmos = std::make_shared<AtmospherePar>(comm, atmosphereParams);
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
	}
	EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
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
		
	}
	catch (...)
	{
		failed = true;
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

	// set parameter
	coupledModel->setPar(0.0001);

	// try to converge
	int maxit = 10;
	for (int i = 0; i != maxit; ++i)
	{
		coupledModel->computeRHS();
		coupledModel->computeJacobian();
	
		std::shared_ptr<Combined_MultiVec> b = coupledModel->getRHS('C');
		b->Scale(-1.0);
		
		double normb = Utils::norm(b);
		
		coupledModel->solve(b);

		std::shared_ptr<Combined_MultiVec> x = coupledModel->getSolution('C');		
		std::shared_ptr<Combined_MultiVec> y = coupledModel->getSolution('C');

		INFO(" ||x||  = " << Utils::norm(stateV) );
		stateV->Update(1.0, *x, 1.0); // x = x + dx;
		
		coupledModel->applyMatrix(*x, *y);

		y->Update(1.0, *b, -1.0);
		y->Scale(1./normb);

		Utils::print(y, "residual");
	
		INFO(" ocean ||r|| / ||b||  = " << Utils::norm(y->First()));
		INFO(" atmos ||r|| / ||b||  = " << Utils::norm(y->Second()));
		INFO(" total ||r|| / ||b||  = " << Utils::norm(y));

		INFO(" ||F|| = " << normb);
		
	}
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
