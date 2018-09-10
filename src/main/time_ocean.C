//=======================================================================
// Thetastepper with the ocean model
//=======================================================================

#include "RunDefinitions.H"
#include "ThetaStepper.H"
#include "Theta.H"       

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
void runOceanModel(RCP<Epetra_Comm> Comm);

//------------------------------------------------------------------
int main(int argc, char **argv)
{
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);

	// run the ocean model
	runOceanModel(Comm);

    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	MPI_Finalize();
}

//------------------------------------------------------------------
void runOceanModel(RCP<Epetra_Comm> Comm)
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

    double comb =
        oceanParams->sublist("THCM").
        sublist("Starting Parameters").get("Combined Forcing", 0.0);

    if (std::abs(comb) < 1e-7)
    {
        WARNING("Nothing will happen without any forcing: par(comb) = "
                << comb, __FILE__, __LINE__);
    }

	// Create parallelized Theta<Ocean> object
    Teuchos::RCP<Theta<Ocean> > oceanTheta =
        Teuchos::rcp(new Theta<Ocean>(Comm, oceanParams));

    // Create parameter object for time stepping
	RCP<Teuchos::ParameterList> timeParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("timestepper_params.xml", timeParams.ptr());
    timeParams->setName("time stepper parameters");

    // Create ThetaStepper
    ThetaStepper<Teuchos::RCP<Theta<Ocean> >, Teuchos::RCP<Teuchos::ParameterList> >
        stepper(oceanTheta, timeParams);

    // Run ThetaStepper
    stepper.run();

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile(profile);

	TIMER_STOP("Total time...");
}

