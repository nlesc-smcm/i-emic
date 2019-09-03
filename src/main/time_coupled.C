//=======================================================================
// Main continuation of the coupled model
//=======================================================================

#include <Teuchos_RCP.hpp>

#include <string>
#include <vector>
#include <memory>

#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Ocean.H"
#include "Atmosphere.H"
#include "SeaIce.H"
#include "CoupledModel.H"
#include "CoupledThetaModel.H"
#include "TransientFactory.H"

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
    
    // Create parallelized ThetaModel<Ocean> object
    std::shared_ptr<ThetaModel<Ocean> > ocean =
        std::make_shared<ThetaModel<Ocean> >(Comm, params[OCEAN], params[TIME]);

    // Create parallelized ThetaModel<Atmosphere> object
    std::shared_ptr<ThetaModel<Atmosphere> > atmos =
        std::make_shared<ThetaModel<Atmosphere> >(Comm, params[ATMOS], params[TIME]);

    // Create parallelized ThetaModel<Atmosphere> object
    std::shared_ptr<ThetaModel<SeaIce> > seaice =
        std::make_shared<ThetaModel<SeaIce> >(Comm, params[SEAICE], params[TIME]);

    // Create CoupledModel
    std::shared_ptr<CoupledModel> coupledModel =
        std::make_shared<CoupledModel>(ocean, atmos, seaice, params[COUPLED]);

    // Create CoupledThetaModel
    std::shared_ptr<CoupledThetaModel> coupledThetaModel =
        std::make_shared<CoupledThetaModel>(*coupledModel, params[TIME]);

    // Create time stepper
    auto stepper = Teuchos::rcp(
        new AdaptiveTransient<decltype(coupledThetaModel)>(coupledThetaModel, params[TIME]));

    TIMER_STOP("Total initialization");

    // run timestepper
    int status = stepper->run();
    if (status != 0)
        ERROR("Timestepper failed", __FILE__, __LINE__);

    TIMER_STOP("Total time...");

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile();
}
