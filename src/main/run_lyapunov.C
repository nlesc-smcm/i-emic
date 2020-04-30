//=======================================================================
// Main continuation of the Lyapunov model
//=======================================================================

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <memory>

#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Continuation.H"
#include "Ocean.H"
#include "LyapunovModel.H"

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

    // Create parameter object for continuation
    RCP<Teuchos::ParameterList> continuationParams =
        Utils::obtainParams("continuation_params.xml", "Continuation");

    // Add the JDQZ parameters
    Utils::obtainParams(continuationParams, "jdqz_params.xml", "JDQZ");

    // Let the continuation parameters dominate over ocean parameters
    Utils::overwriteParameters(oceanParams, continuationParams);

    // Create parallelized Ocean object
    RCP<LyapunovModel<Ocean> > ocean = Teuchos::rcp(
        new LyapunovModel<Ocean>(Comm, oceanParams));

    // Create continuation
    Continuation<RCP<LyapunovModel<Ocean>>> continuation(ocean, continuationParams);

    // Run continuation
    int status = continuation.run();
    assert(status == 0);

    TIMER_STOP("Total time...");

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile();
}

