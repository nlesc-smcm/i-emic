//=======================================================================
// Main continuation of the coupled model
//=======================================================================

#include <Teuchos_RCP.hpp>

#include <string>
#include <vector>
#include <memory>

#include "GlobalDefinitions.H"

#include "ComplexVector.H"
#include "JDQZInterface.H"
#include "jdqz.hpp"

#include "Continuation.H"
#include "Ocean.H"
#include "Atmosphere.H"
#include "SeaIce.H"
#include "CoupledModel.H"
#include "ThetaStepper.H"
#include "Theta.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
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
    std::vector<std::string> files = {"ocean_params.xml",
                                      "atmosphere_params.xml",
                                      "seaice_params.xml",
                                      "coupledmodel_params.xml",
                                      "continuation_params.xml",
                                      "jdqz_params.xml",
                                      "timestepper_params.xml"};

    std::vector<std::string> names = {"Ocean parameters",
                                      "Atmosphere parameters",
                                      "Sea ice parameters",
                                      "CoupledModel parameters",
                                      "Continuation parameters",
                                      "JDQZ parameters",
                                      "timestepper_params.xml"};
    
    enum Ident { OCEAN, ATMOS, SEAICE, COUPLED, CONT, EIGEN, TIME};

    std::vector<Teuchos::RCP<Teuchos::ParameterList> > params;

    for (int i = 0; i != (int) files.size(); ++i)
        params.push_back(Utils::obtainParams(files[i], names[i]));

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

    // check params[OCEAN]
    double comb =
        params[OCEAN]->sublist("THCM").
        sublist("Starting Parameters").get("Combined Forcing", 0.0);
    
    if (std::abs(comb) < 1e-7)
    {
        WARNING("Nothing will happen without any forcing: par(comb) = "
                << comb, __FILE__, __LINE__);
    }
    
    // Create parallelized Theta<Ocean> object
    std::shared_ptr<Theta<Ocean> > ocean =
        std::make_shared<Theta<Ocean> >(Comm, params[OCEAN]);

    // Create parallelized Theta<Atmosphere> object
    std::shared_ptr<Theta<Atmosphere> > atmos =
        std::make_shared<Theta<Atmosphere> >(Comm, params[ATMOS]);

    // Create parallelized Theta<Atmosphere> object
    std::shared_ptr<Theta<SeaIce> > seaice =
        std::make_shared<Theta<SeaIce> >(Comm, params[SEAICE]);

    // Create CoupledModel
    std::shared_ptr<CoupledModel> coupledModel =
        std::make_shared<CoupledModel>(ocean, atmos, seaice, params[COUPLED]);

    // Create ThetaStepper
    ThetaStepper<std::shared_ptr<CoupledModel>, Teuchos::RCP<Teuchos::ParameterList> >
        stepper(coupledModel, params[TIME]);


    TIMER_STOP("Total initialization");

    // run timestepper
    int status = stepper.run();
    assert(status == 0);

    TIMER_STOP("Total time...");

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile();
}
