//=======================================================================
// Main continuation of the coupled model
//=======================================================================

#include "RunDefinitions.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

using JDQZsolver = JDQZ<JDQZInterface<std::shared_ptr<CoupledModel>, 
                                      ComplexVector<Combined_MultiVec> > >;

void runCoupledModel(RCP<Epetra_Comm> Comm);

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    //  - MPI
    //  - output files
    //  - returns Trilinos' communicator Epetra_Comm
    RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);

    runCoupledModel(Comm);

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile(profile);

    //--------------------------------------------------------
    // Finalize MPI
    //--------------------------------------------------------
    MPI_Finalize();
}

//------------------------------------------------------------------
void runCoupledModel(RCP<Epetra_Comm> Comm)
{
    TIMER_START("Total time...");

    TIMER_START("Total initialization");
    //------------------------------------------------------------------
    // Check if outFile is specified
    if (outFile == Teuchos::null)
        throw std::runtime_error("ERROR: Specify output streams");

    // First we create a bunch of parameterlists
    
    // Create parameter object for Ocean
    RCP<Teuchos::ParameterList> oceanParams =
        obtainParams("ocean_params.xml", "Ocean parameters");

    // Create parameter object for Atmosphere
    RCP<Teuchos::ParameterList> atmosphereParams =
        obtainParams("atmosphere_params.xml", "Atmosphere parameters"); 

    // Create parameter object for CoupledModel
    RCP<Teuchos::ParameterList> coupledmodelParams =
        obtainParams("coupledmodel_params.xml", "CoupledModel parameters"); 

    // Create parameter object for Continuation
    RCP<Teuchos::ParameterList> continuationParams =
        obtainParams("continuation_params.xml", "Continuation parameters"); 

    // Create parameter object for JDQZ
    RCP<Teuchos::ParameterList> jdqzParams =
        obtainParams("jdqz_params.xml", "JDQZ parameters"); 

    // The Continuation and CoupledModel parameterlists overwrite settings
    Utils::overwriteParameters(oceanParams,        coupledmodelParams);
    Utils::overwriteParameters(atmosphereParams,   coupledmodelParams);

    Utils::overwriteParameters(oceanParams,        continuationParams);
    Utils::overwriteParameters(atmosphereParams,   continuationParams);
    Utils::overwriteParameters(coupledmodelParams, continuationParams);
    
    // Create parallelized Ocean object
    std::shared_ptr<Ocean> ocean = std::make_shared<Ocean>(Comm, oceanParams);

    // Create parallelized Atmosphere object
    std::shared_ptr<AtmospherePar> atmos =
        std::make_shared<AtmospherePar>(Comm, atmosphereParams);

    // Create CoupledModel
    std::shared_ptr<CoupledModel> coupledModel =
        std::make_shared<CoupledModel>(ocean, atmos, coupledmodelParams);

    // Create JDQZ generalized eigenvalue solver
    Combined_MultiVec t = *coupledModel->getSolution('C');
    t.PutScalar(0.0);
    ComplexVector<Combined_MultiVec> z(t);
    
    JDQZInterface<std::shared_ptr<CoupledModel>,
                  ComplexVector<Combined_MultiVec> > matrix(coupledModel, z);
    
    std::shared_ptr<JDQZsolver> jdqz = std::make_shared<JDQZsolver>(matrix, z);
    jdqz->setParameters(*jdqzParams);
    
    // Create Continuation
    Continuation<std::shared_ptr<CoupledModel>, RCP<Teuchos::ParameterList> >
        continuation(coupledModel, continuationParams);

    // Couple JDQZ to continuation
    continuation.setEigenSolver(jdqz);

    TIMER_STOP("Total initialization");

    // Run continuation
    continuation.run();
    
    TIMER_STOP("Total time...");
}
