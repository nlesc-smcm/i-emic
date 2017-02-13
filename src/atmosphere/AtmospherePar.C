#include "AtmospherePar.H"
#include "AtmosphereDefinitions.H"

//==================================================================
// Constructor
AtmospherePar::AtmospherePar(Teuchos::RCP<Epetra_Comm> comm, ParameterList params)
	:
	params_          (params),
	comm_            (comm),
	n_               (params->get("Global Grid-Size n", 16)),
	m_               (params->get("Global Grid-Size m", 16)),
	l_               (params->get("Global Grid-Size l", 1)),
	periodic_        (params->get("Periodic", false)),
	inputFile_       (params->get("Input file", "atmos_input.h5")),
	outputFile_      (params->get("Output file", "atmos_output.h5")),
	loadState_       (params->get("Load state", false)),
	saveState_       (params->get("Save state", false))
{
	INFO("AtmospherePar: constructor...");
	// Define degrees of freedom
	dof_ = ATMOS_NUN_;
	
	// Define domain
	xmin_ = params->get("Global Bound xmin", 286.0) * PI_ / 180.0;
	xmax_ = params->get("Global Bound xmax", 350.0) * PI_ / 180.0;
	ymin_ = params->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
	ymax_ = params->get("Global Bound ymax", 74.0)  * PI_ / 180.0;

	// Create domain object
	domain_ = Teuchos::rcp(new TRIOS::Domain(n_, m_, l_, dof_,
											 xmin_, xmax_, ymin_, ymax_,
											 periodic_, 1.0, comm_));

	// Compute 2D decomposition
	domain_->Decomp2D();

	// Obtain local dimensions
	double xminloc = domain_->XminLoc();
	double xmaxloc = domain_->XmaxLoc();
	double yminloc = domain_->YminLoc();
	double ymaxloc = domain_->YmaxLoc();

	// Obtain local grid dimensions
	int nloc = domain_->LocalN();
	int mloc = domain_->LocalM();
	int lloc = domain_->LocalL();

	// Obtain overlapping and non-overlapping maps	
	assemblyMap_ = domain_->GetAssemblyMap();
	standardMap_ = domain_->GetStandardMap();

	// Create overlapping and non-overlapping states
	state_      = Teuchos::rcp(new Epetra_Vector(*standardMap_));
	rhs_        = Teuchos::rcp(new Epetra_Vector(*standardMap_));
	localState_ = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
	localRHS_   = Teuchos::rcp(new Epetra_Vector(*assemblyMap_)); 

	// Create local Atmosphere object
	// Periodicity is handled by Atmosphere if there is a single
	// core in the x-direction. 
	Teuchos::RCP<Epetra_Comm> xComm = domain_->GetProcRow(0);
	bool perio = (periodic_ && xComm->NumProc() == 1);
	atmos_ = std::make_shared<Atmosphere>(nloc, mloc, lloc, perio,
										  xminloc, xmaxloc, yminloc, ymaxloc,
										  params_);
	
	INFO("AtmospherePar: constructor done");
}

//==================================================================
void AtmospherePar::computeRHS()
{
	INFO("AtmosepherePar: computeRHS...");
	atmos_->computeRHS();
	INFO("AtmosepherePar: computeRHS done");
}
