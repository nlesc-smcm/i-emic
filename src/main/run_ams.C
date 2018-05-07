//=======================================================================
// Main ams of the ocean model
//=======================================================================

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Ocean.H"
#include "AMS.H"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

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

    // Create parameter object for ams
    RCP<Teuchos::ParameterList> amsParams = rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("ams_params.xml", amsParams.ptr());
    amsParams->setName("AMS parameters");

    // Let the ams parameters dominate over ocean parameters
    Utils::overwriteParameters(oceanParams, amsParams);

    std::string stateA = amsParams->get("State A", "");
    std::string stateB = amsParams->get("State B", "");
    oceanParams->set("Input file", stateA);

    RCP<Ocean> ocean = Teuchos::rcp(new Ocean(Comm, oceanParams));

    RCP<Epetra_Vector> sol1 = ocean->getState('C');
    Utils::load(sol1, stateA);
    RCP<Epetra_Vector> sol2 = ocean->getState('C');
    Utils::load(sol2, stateB);

    // EpetraExt::MultiVectorToMatrixMarketFile("sol1.mtx", *sol1);
    // EpetraExt::MultiVectorToMatrixMarketFile("sol2.mtx", *sol2);

    // Create ams
    AMS<RCP<Ocean> > ams(ocean, amsParams, sol1, sol2);

    // ocean->computeJacobian();
    // EpetraExt::RowMatrixToMatrixMarketFile("A.mtx", *ocean->getJacobian());
    // EpetraExt::MultiVectorToMatrixMarketFile("M.mtx", *ocean->getDiagB());

    ocean->computeForcing();
    EpetraExt::RowMatrixToMatrixMarketFile("F.mtx", *ocean->getForcing());

    ams.run();

    TIMER_STOP("Total time...");

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile();
}

