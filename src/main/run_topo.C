//=======================================================================
// Continuation of topographies in an ocean model
//=======================================================================

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "GlobalDefinitions.H"

#include "Continuation.H"
#include "Ocean.H"
#include "Topo.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    //  - MPI
    //  - output files
    //  - returns Trilinos' communicator Epetra_Comm
    RCP<Epetra_Comm> comm = initializeEnvironment(argc, argv);

    if (outFile == Teuchos::null)
        throw std::runtime_error("ERROR: Specify output streams");

    // Create parameters parallel Ocean
    RCP<Teuchos::ParameterList> oceanParams =
        rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());

    // Create parameters for topography continuation class
    RCP<Teuchos::ParameterList> topoParams = rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("topo_params.xml", topoParams.ptr());

    // Create parameter object for continuation
    RCP<Teuchos::ParameterList> continuationParams =
        rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("continuation_params.xml",
                                continuationParams.ptr());

    // Create ocean
    Teuchos::RCP<Ocean> ocean = Teuchos::rcp(new Ocean(comm, oceanParams));

    // Create topography continuation class
    RCP<Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> > > topo =
        rcp(new Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> >
            (ocean, topoParams));

    // Create continuation
    Continuation<RCP<Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> > >,
                 RCP<Teuchos::ParameterList> >
        continuation(topo, continuationParams);

    // Run
    int nMasks    = topo->nMasks();
    int startMask = topo->startMaskIdx();

    for (int maskIdx = startMask; maskIdx != nMasks-1; maskIdx++)
    {
        topo->setMaskIndex(maskIdx);
        topo->initialize();

        TIMER_START("  TOPO:  Predictor");
        topo->predictor();
        TIMER_STOP ("  TOPO:  Predictor");

        TIMER_START("  TOPO:  Homotopy Continuation");
        int status = continuation.run();
        TIMER_STOP ("  TOPO:  Homotopy Continuation");
        if (status == 0)
            ERROR("Continuation failed", __FILE__, __LINE__);
        
        topo->setPar(1.0);
        
        if (argc == 1)
            topo->postProcess();
        else
            std::cout << " arguments given, no usage yet" << std::endl;
    }

    ocean = Teuchos::null;
    topo  = Teuchos::null;

    // print the profile
    if (comm->MyPID() == 0)
        printProfile();

    comm->Barrier();
    MPI_Finalize();
}
