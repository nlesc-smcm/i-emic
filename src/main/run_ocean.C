//=======================================================================
// Main continuation of the ocean model
//=======================================================================

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <memory>

#include "GlobalDefinitions.H"
#include "Utils.H"

#include "ComplexVector.H"
#include "JDQZInterface.H"
#include "jdqz.hpp"

#include "Continuation.H"
#include "Ocean.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
using JDQZsolver = JDQZ<JDQZInterface<Teuchos::RCP<Ocean>,
                                      ComplexVector<Epetra_Vector> > >;

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

    // Create parameter object for continuation
    RCP<Teuchos::ParameterList> continuationParams = rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("continuation_params.xml", continuationParams.ptr());
    continuationParams->setName("Continuation parameters");

    // Create parameter object for JDQZ
    RCP<Teuchos::ParameterList> jdqzParams =
        Utils::obtainParams("jdqz_params.xml", "JDQZ parameters");

    // Let the continuation parameters dominate over ocean parameters
    Utils::overwriteParameters(oceanParams, continuationParams);

    // Create parallelized Ocean object
    RCP<Ocean> ocean = Teuchos::rcp(new Ocean(Comm, oceanParams));

    // Create continuation
    Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
        continuation(ocean, continuationParams);

    // Create JDQZ generalized eigenvalue solver:
    // 1) Create a vector with complex arithmetic based on an Epetra_Vector
    Epetra_Vector t = *ocean->getSolution('C');
    t.PutScalar(0.0);
    ComplexVector<Epetra_Vector> z(t);

    // 2) Create JDQZ (mass) matrix and preconditioning interface
    JDQZInterface<Teuchos::RCP<Ocean>,
                  ComplexVector<Epetra_Vector> > matrix(ocean, z);

    // 3) Create JDQZ solver
    std::shared_ptr<JDQZsolver> jdqz = std::make_shared<JDQZsolver>(matrix, z);
    jdqz->setParameters(*jdqzParams);

    // Couple JDQZ to continuation
    continuation.setEigenSolver(jdqz);

    // Run continuation
    int status = continuation.run();
    if (status != 0)
        ERROR("Continuation failed", __FILE__, __LINE__);

    TIMER_STOP("Total time...");

    // print the profile
    if (Comm->MyPID() == 0)
    {
        printProfile();
        jdqz->printProfile("jdqz_profile");
    }
}
