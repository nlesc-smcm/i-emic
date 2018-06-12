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

    // First we create a bunch of parameterlists based on a bunch of
    // xml files.
    std::vector<string> files = {"ocean_params.xml",
                                 "atmosphere_params.xml",
                                 "seaice_params.xml",
                                 "coupledmodel_params.xml",
                                 "continuation_params.xml",
                                 "jdqz_params.xml"};

    std::vector<string> names = {"Ocean parameters",
                                 "Atmosphere parameters",
                                 "Sea ice parameters",
                                 "CoupledModel parameters",
                                 "Continuation parameters",
                                 "JDQZ parameters"};
    
    enum Ident { OCEAN, ATMOS, SEAICE, COUPLED, CONT, EIGEN };

    std::vector<Teuchos::RCP<Teuchos::ParameterList> > params;

    for (int i = 0; i != (int) files.size(); ++i)
        params.push_back(obtainParams(files[i], names[i]));

    INFO('\n' << "Overwriting:");
    // Allow dominant parameterlists. Not that this trick uses a
    // 'flattened' hierarchy. The Continuation and CoupledModel
    // parameterlists are allowed to overwrite settings.
    Utils::overwriteParameters(params[OCEAN],  params[COUPLED]);
    Utils::overwriteParameters(params[ATMOS],  params[COUPLED]);
    Utils::overwriteParameters(params[SEAICE], params[COUPLED]);

    Utils::overwriteParameters(params[OCEAN],  params[CONT]);
    Utils::overwriteParameters(params[ATMOS],  params[CONT]);
    Utils::overwriteParameters(params[SEAICE], params[CONT]);

    Utils::overwriteParameters(params[COUPLED], params[CONT]);
    
    // Create parallelized Ocean object
    std::shared_ptr<Ocean> ocean = std::make_shared<Ocean>(Comm, params[OCEAN]);

    // Create parallelized Atmosphere object
    std::shared_ptr<Atmosphere> atmos =
        std::make_shared<Atmosphere>(Comm, params[ATMOS]);

    // Create parallelized Atmosphere object
    std::shared_ptr<SeaIce> seaice =
        std::make_shared<SeaIce>(Comm, params[SEAICE]);

    // Create CoupledModel
    std::shared_ptr<CoupledModel> coupledModel =
        std::make_shared<CoupledModel>(ocean, atmos, seaice, params[COUPLED]);

    // Create JDQZ generalized eigenvalue solver
    Combined_MultiVec t = *coupledModel->getSolution('C');
    t.PutScalar(0.0);
    ComplexVector<Combined_MultiVec> z(t);
    
    JDQZInterface<std::shared_ptr<CoupledModel>,
                  ComplexVector<Combined_MultiVec> > matrix(coupledModel, z);
    
    std::shared_ptr<JDQZsolver> jdqz = std::make_shared<JDQZsolver>(matrix, z);
    jdqz->setParameters(*params[EIGEN]);
    
    // Create Continuation
    Continuation<std::shared_ptr<CoupledModel>, RCP<Teuchos::ParameterList> >
        continuation(coupledModel, params[CONT]);

    // Couple JDQZ to continuation
    continuation.setEigenSolver(jdqz);

    TIMER_STOP("Total initialization");

    // Run continuation
    continuation.run();

    // print the profile
    if (Comm->MyPID() == 0)
    {
        printProfile(profile);
        jdqz->printProfile("jdqz_profile");
    }    
    
    TIMER_STOP("Total time...");
}
