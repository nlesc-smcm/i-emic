//=======================================================================
// Main continuation of the ocean model
//=======================================================================

#include "RunDefinitions.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
void run(RCP<Epetra_Comm> Comm);

//------------------------------------------------------------------
int main(int argc, char **argv)
{
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);

	// run the ocean model
	run(Comm);

	// print the profile
	if (Comm->MyPID() == 0)
		printProfile(profile);
	
    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	MPI_Finalize();	
}

//------------------------------------------------------------------
void run(RCP<Epetra_Comm> Comm)
{
	TIMER_START("Total time...");

	//------------------------------------------------------------------
	// Check if outFile is specified
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");	
	
	// Create parameter object for Ocean
	RCP<Teuchos::ParameterList> oceanParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
    oceanParams->setName("Ocean parameters");
    
 	// Create parameter object for continuation
	RCP<Teuchos::ParameterList> continuationParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("continuation_params.xml", continuationParams.ptr());
    continuationParams->setName("Continuation parameters");

    // Let the continuation parameters dominate over ocean parameters
    INFO("Overwriting:");
    Utils::overwriteParameters(oceanParams, continuationParams);

	// Create parallelized Ocean object
	RCP<Ocean> ocean = Teuchos::rcp(new Ocean(Comm, oceanParams));
	
	// Create continuation
	Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
		continuation(ocean, continuationParams);
	
	// Run continuation
	continuation.run();

	// Calculate eigenvalues and eigenvectors
	// JDQZMAT<Ocean> matrix;
	// JDQZ<JDQZMAT<Ocean> > jdqz(matrix);

	//------------------------------------------------------------------
	TIMER_STOP("Total time...");		
}

