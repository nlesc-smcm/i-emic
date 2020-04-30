//=======================================================================
// Main ams of the ocean model
//=======================================================================

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Ocean.H"
#include "TransientFactory.H"

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
    RCP<Teuchos::ParameterList> oceanParams =
        Utils::obtainParams("ocean_params.xml", "Ocean");

    Utils::obtainParams(oceanParams, "solver_params.xml", "Belos Solver");

    // Create parameter object for ams
    RCP<Teuchos::ParameterList> amsParams =
        Utils::obtainParams("ams_params.xml", "AMS parameters");

    // Let the ams parameters dominate over ocean parameters
    Utils::overwriteParameters(oceanParams, amsParams);

    std::string stateA = amsParams->get("State A", "");
    std::string stateB = amsParams->get("State B", "");
    std::string stateD = amsParams->get("Unstable state", "");
    oceanParams->set("Input file", stateA);
    oceanParams->set("Load state", true);

    RCP<Ocean> ocean = Teuchos::rcp(new Ocean(Comm, oceanParams));

    RCP<Epetra_Vector> sol1 = ocean->getState('C');
    Utils::load(sol1, stateA);
    RCP<Epetra_Vector> sol2 = ocean->getState('C');
    Utils::load(sol2, stateB);
    RCP<Epetra_Vector> sol3 = Teuchos::null;
    if (stateD != "")
    {
        sol3 = ocean->getState('C');
        Utils::load(sol3, stateD);
    }

    bool writeMatrices = amsParams->get("write matrices", false);

    if (writeMatrices)
    {
        EpetraExt::MultiVectorToMatrixMarketFile("sol1.mtx", *sol1);
        EpetraExt::MultiVectorToMatrixMarketFile("sol2.mtx", *sol2);
        if (sol3 != Teuchos::null)
            EpetraExt::MultiVectorToMatrixMarketFile("sol3.mtx", *sol3);

        *ocean->getState('V') = *sol1;
        ocean->computeJacobian();
        EpetraExt::RowMatrixToMatrixMarketFile("A1.mtx", *ocean->getJacobian());
        *ocean->getState('V') = *sol2;
        ocean->computeJacobian();
        EpetraExt::RowMatrixToMatrixMarketFile("A2.mtx", *ocean->getJacobian());
        EpetraExt::MultiVectorToMatrixMarketFile("M.mtx", *ocean->getMassMat());

        ocean->computeForcing();
        EpetraExt::RowMatrixToMatrixMarketFile("F.mtx", *ocean->getForcing());
    }

    // Create ams
    auto ams = TransientFactory(ocean, amsParams, sol1, sol2, sol3);

    ams->run();

    TIMER_STOP("Total time...");

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile();
}

