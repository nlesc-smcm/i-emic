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
	dim_ = n_ * m_ * l_ * dof_;
	
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

	// Create assembly state
	domain_->Solve2Assembly(*state_, *localState_);
	
	// local problem size
	int numMyElements = assemblyMap_->NumMyElements();

	std::shared_ptr<std::vector<double> > localState =
		std::make_shared<std::vector<double> >(numMyElements, 0.0);

	localState_->ExtractCopy(&(*localState)[0], numMyElements);

	atmos_->setState(localState);

	// compute local rhs and check bounds
	atmos_->computeRHS();
 	std::shared_ptr<std::vector<double> > rhs = atmos_->getRHS('V');
	
	if ((int) rhs->size() != numMyElements)
	{
		ERROR("RHS incorrect size", __FILE__, __LINE__);
	}

	// obtain view
	double *rhs_tmp;
	localRHS_->ExtractView(&rhs_tmp);

	// fill view
	for (int i = 0; i != numMyElements; ++i)
	{
		rhs_tmp[i] = (*rhs)[i];
	}

	// set datamember
	domain_->Assembly2Solve(*localRHS_, *rhs_);

#ifdef DEBUGGING_NEW
	std::ofstream file1, file2;
	file1.open("rhs_epetra" + std::to_string(comm_->MyPID()) + ".txt");
	file2.open("rhs_atmos" + std::to_string(comm_->MyPID()) + ".txt");
	rhs_->Print(file1);
	for (auto &it: *rhs)
		file2 << it << " ";
	file1.close();
	file2.close();
	double nrm;
	rhs_->Norm2(&nrm);
	INFO("AtmospherePar rhs norm: " << nrm);
#endif
		
	INFO("AtmospherePar: computeRHS done");
}

//==================================================================
void AtmospherePar::idealized()
{
	// initialize local rhs with idealized values
	atmos_->idealized();

	// local problem size
	int numMyElements = assemblyMap_->NumMyElements();

	// obtain view of assembly state
	double *state_tmp;
	localState_->ExtractView(&state_tmp);

	// obtain local state and check bounds
	std::shared_ptr<std::vector<double> > state = atmos_->getState('V');
	if ((int) state->size() != numMyElements)
	{
		ERROR("state incorrect size", __FILE__, __LINE__);
	}

	// fill assembly view with local state
	for (int i = 0; i != numMyElements; ++i)
	{
		state_tmp[i] = (*state)[i];
	}
	
	// set solvemap state
	domain_->Assembly2Solve(*localState_, *state_);

#ifdef DEBUGGING_NEW
	std::ofstream file1,file2;
	file1.open("state_epetra" + std::to_string(comm_->MyPID()) + ".txt");
	file2.open("state_atmos" + std::to_string(comm_->MyPID()) + ".txt");
	state_->Print(file1);
	for (auto &it: *state)
		file2 << it << " ";

	file1.close();
	file2.close();
	double nrm;
	state_->Norm2(&nrm);
	INFO("AtmospherePar idealized state norm: " << nrm);
#endif
	
}
