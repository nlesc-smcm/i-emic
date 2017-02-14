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

	// Obtain special maps
	// depth-averaged, single unknown for ocean surface temperature
	standardSurfaceMap_ = domain_->CreateStandardMap(1, true);
	assemblySurfaceMap_ = domain_->CreateAssemblyMap(1, true);

	// Create overlapping and non-overlapping vectors
	state_      = Teuchos::rcp(new Epetra_Vector(*standardMap_));
	rhs_        = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    sst_       	= Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

	localState_ = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
	localRHS_   = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
	localSST_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));

	// create graph
	createMatrixGraph();
		
	// jac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *matrixGraph_));

	// Periodicity is handled by Atmosphere if there is a single
	// core in the x-direction. 
	Teuchos::RCP<Epetra_Comm> xComm = domain_->GetProcRow(0);
	bool perio = (periodic_ && xComm->NumProc() == 1);
	
	// Create local Atmosphere object
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
 	std::shared_ptr<std::vector<double> > localRHS = atmos_->getRHS('V');
	
	if ((int) localRHS->size() != numMyElements)
	{
		ERROR("RHS incorrect size", __FILE__, __LINE__);
	}

	// obtain view
	double *rhs_tmp;
	localRHS_->ExtractView(&rhs_tmp);

	// fill view
	for (int i = 0; i != numMyElements; ++i)
	{
		rhs_tmp[i] = (*localRHS)[i];
	}

	// set datamember
	domain_->Assembly2Solve(*localRHS_, *rhs_);
		
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
	
}

//==================================================================
Teuchos::RCP<Epetra_Vector> AtmospherePar::getVector(char mode, Teuchos::RCP<Epetra_Vector> vec)
{
	if (mode == 'C') // copy
	{
    	Teuchos::RCP<Epetra_Vector> copy = Teuchos::rcp(new Epetra_Vector(*vec));
		return copy;
    }
	else if (mode == 'V')
		return vec;
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}
	
}


//==================================================================
void AtmospherePar::setOceanTemperature(Teuchos::RCP<Epetra_Vector> in)
{
	if (!(in->Map().SameAs(*standardSurfaceMap_)))
	{
		ERROR("Map of ocean surface input vector not same as surface map", __FILE__, __LINE__);
	}

	// assign to our own datamember
	sst_ = in;

	// create assembly 
	domain_->Solve2Assembly(*sst_, *localSST_);

	// local vector size
	int numMyElements = assemblySurfaceMap_->NumMyElements();

	std::shared_ptr<std::vector<double> > localSST =
		std::make_shared<std::vector<double> >(numMyElements, 0.0);

	localSST_->ExtractCopy(&(*localSST)[0], numMyElements);
	atmos_->setOceanTemperature(*localSST);	
}


//==================================================================
void AtmospherePar::computeJacobian()
{
	// compute jacobian in local atmosphere
	atmos_->computeJacobian();

	// obtain CRS matrix from local atmosphere
	std::shared_ptr<Atmosphere::CRSMat> localJac =
		atmos_->getJacobian();

	// max nonzeros per row
	const int maxnnz = ATMOS_NUN_ * ATMOS_NP_ + 1;

	// indices array
	int indices[maxnnz];
	
	// values array
	double values[maxnnz];
	
	int numMyElements = assemblyMap_->NumMyElements();
	assert(numMyElements == (int) (*localJac)["beg"].size() - 1);

	int index, numentries;
	for (int i = 0; i < numMyElements; ++i)
	{
		if (!domain_->IsGhost(i))
		{
			index = (*localJac)["beg"][i]; // beg contains 1-based indices!
			numentries = (*localJac)["beg"][i+1] - index;
			for (int j = 0; j < numentries; ++j)
			{
				indices[j] = assemblyMap_->GID((*localJac)["jco"][index-1+j]);
				values[j] = (*localJac)["co"][index-1+j];
			}
			
//			int ierr = jac_->ReplaceGlobalValues(assemblyMap_->GID(i),
//													 numentries,
//													 values, indices);
		}
	}
}

//==================================================================
// Very similar to the THCM function
// Create a graph to initialize the Jacobian
void AtmospherePar::createMatrixGraph()
{
	int n    = domain_->LocalN();
	int m    = domain_->LocalM();
	int l    = domain_->LocalL();
	int ndim = standardMap_->NumMyElements();
	int *numEntriesPerRow = new int[ndim];

	for (int k = 1; k <= l; k++)
		for (int j = 1; j <= m; j++)
			for (int i = 1; i <= n; i++)
			{
				// get 1-based row from local atmos, convert to 0-based
				int lidU = atmos_->find_row(i, j, k, ATMOS_TT_) - 1;

				// get global id of row
				int gidU = assemblyMap_->GID(lidU);

				if (standardMap_->MyGID(gidU)) // otherwise: ghost cell
				{
					// obtain local id
					int lid0 = standardMap_->LID(gidU) - 1;
					for (int xx = 1; xx <= dof_; ++xx)
					{
						// better safe than sorry
						numEntriesPerRow[lid0+xx] = ATMOS_NP_;
					}
				}
			}
	
}
