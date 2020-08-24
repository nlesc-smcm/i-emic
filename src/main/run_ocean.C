//=======================================================================
// Main continuation of the ocean model
//=======================================================================

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <memory>

#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Continuation.H"
#include "Ocean.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
void runOceanModel(RCP<Epetra_Comm> Comm);

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    Teuchos::CommandLineProcessor clp(false);
    clp.setDocString("This program performs an ocean-only run of I-EMIC\n");

    bool describe_parameters = false;
    clp.setOption("describe-parameters", "disable-describe-parameters", &describe_parameters,
                  "Describe all possible parameter values");

    Teuchos::CommandLineProcessor::EParseCommandLineReturn ret = clp.parse(argc, argv);

    if (ret == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
        return 0;
    if (ret != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
        return -1;

    if (describe_parameters)
    {
        Teuchos::ParameterList::PrintOptions opts;
        opts.showDoc(true);
        Teuchos::ParameterList params;
        params.sublist("Ocean") = Ocean::getDefaultInitParameters();
        params.sublist("Continuation") = Continuation<RCP<Ocean>>::getDefaultInitParameters();
        params.print(std::cout, opts);
        return 0;
    }

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

#ifdef HAVE_JDQZPP
    // Add the JDQZ parameters
    Utils::obtainParams(continuationParams, "jdqz_params.xml", "JDQZ");
#endif

    // Let the continuation parameters dominate over ocean parameters
    Utils::overwriteParameters(oceanParams, continuationParams);

    // Create parallelized Ocean object
    RCP<Ocean> ocean = Teuchos::rcp(new Ocean(Comm, oceanParams));

    // Create continuation
    Continuation<RCP<Ocean>> continuation(ocean, continuationParams);

    // Run continuation
    int status = continuation.run();
    if (status != 0)
        ERROR("Continuation failed", __FILE__, __LINE__);

    TIMER_STOP("Total time...");

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile();
}
