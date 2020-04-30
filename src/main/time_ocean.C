//=======================================================================
// Thetastepper with the ocean model
//=======================================================================

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "GlobalDefinitions.H"

#include "Ocean.H"
#include "TransientFactory.H"

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
    RCP<Teuchos::ParameterList> oceanParams =
        Utils::obtainParams("ocean_params.xml", "Ocean");

    Utils::obtainParams(oceanParams, "solver_params.xml", "Belos Solver");

    double comb =
        oceanParams->sublist("THCM").
        sublist("Starting Parameters").get("Combined Forcing", 0.0);

    if (std::abs(comb) < 1e-7)
    {
        WARNING("Nothing will happen without any forcing: par(comb) = "
                << comb, __FILE__, __LINE__);
    }

    // Create Ocean object
    Teuchos::RCP<Ocean> ocean =
        Teuchos::rcp(new Ocean(Comm, oceanParams));

    // Create parameter object for time stepping
    RCP<Teuchos::ParameterList> timeParams =
        Utils::obtainParams("timestepper_params.xml", "time stepper parameters");

    // Create time stepper
    auto stepper = TransientFactory(ocean, timeParams);

    // Run time stepper
    int status = stepper->run();
    if (status != 0)
        ERROR("Timestepper failed", __FILE__, __LINE__);

    TIMER_STOP("Total time...");

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile();
}

