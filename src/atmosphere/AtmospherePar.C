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
	sol_        = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    sst_       	= Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

	localState_ = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
	localRHS_   = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
	localSol_   = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
	localSST_   = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));

	// create graph
	createMatrixGraph();
		
	jac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *matrixGraph_));

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
	// set all entries to zero
	CHECK_ZERO(jac_->PutScalar(0.0));

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

	// check size
	int numMyElements = assemblyMap_->NumMyElements();
	assert(numMyElements == (int) (*localJac)["beg"].size() - 1);

	// loop over local elements
	int index, numentries;
	for (int i = 0; i < numMyElements; ++i)
	{
		// ignore ghost rows
		if (!domain_->IsGhost(i, ATMOS_NUN_))
		{
			// obtain indices and values from CRS container
			index = (*localJac)["beg"][i]; // beg contains 1-based indices!
			numentries = (*localJac)["beg"][i+1] - index;
			for (int j = 0; j < numentries; ++j)
			{
				indices[j] = assemblyMap_->GID((*localJac)["jco"][index-1+j] - 1);
				values[j]  = (*localJac)["co"][index-1+j];
			}

			// put values in Jacobian
			int ierr = jac_->ReplaceGlobalValues(assemblyMap_->GID(i),
												 numentries,
												 values, indices);
			// debugging
			if (ierr != 0)
			{
				for (int ii = 0; ii < numentries; ++ii)
				{
					std::cout << "proc" << comm_->MyPID()
							  << " entries: (" << indices[ii]
							  << " " << values[ii] << ")" <<  std::endl;
					
					std::cout << "proc" << comm_->MyPID() << " "
							  << assemblyMap_->GID(i) << std::endl;
					
					INFO(" debug info: " << indices[ii] << " " << values[ii]);
				}
				
				INFO(" GRID: "<< assemblyMap_->GID(i));
				INFO(" number of entries: " << numentries);
				INFO(" numMyElements: " << numMyElements);
				INFO(" is ghost " << domain_->IsGhost(i, ATMOS_NUN_));
				INFO(" maxnnz: " << maxnnz);

				CHECK_ZERO(jac_->ExtractGlobalRowCopy(assemblyMap_->GID(i),
													  maxnnz, numentries,
													  values, indices));
				INFO("\noriginal row: ");
				INFO("number of entries: "<<numentries);

				for (int ii = 0; ii < numentries; ++ii)
				{
					std::cout << "proc" << comm_->MyPID()
							  << " entries: (" << indices[ii]
							  << " " << values[ii] << ")" <<  std::endl;
					INFO(" debug info: " << indices[ii] << " " << values[ii]);
				}

				INFO ("Error in ReplaceGlobalValues: " << ierr);
				ERROR("Error in ReplaceGlobalValues", __FILE__, __LINE__);
			}
		}
	}

	// Finalize matrix
	CHECK_ZERO(jac_->FillComplete());
}

//==================================================================
void AtmospherePar::applyMatrix(Teuchos::RCP<Epetra_Vector> const &in,
								Teuchos::RCP<Epetra_Vector> &out)
{
	jac_->Apply(*in, *out);
}

//==================================================================
void AtmospherePar::applyPrecon(Teuchos::RCP<Epetra_Vector> const &in,
								Teuchos::RCP<Epetra_Vector> &out)
{
	// Todo
}

//==================================================================
void AtmospherePar::solve(Teuchos::RCP<Epetra_Vector> const &b)
{
	solveSubDomain(b);
}

//==================================================================
void AtmospherePar::solveSubDomain(Teuchos::RCP<Epetra_Vector> const &b)
{
	if (!(rhs->Map().SameAs(*standardMap_)))
	{
		ERROR("AtmospherePar::solve, map of b not same as standard map",
			  __FILE__, __LINE__);
	}

	// obtain assembly vector, put b in sol
	domain_->Solve2Assembly(*b, *localSol_);

	// local vector size
	int numMyElements = assemblyMap_->NumMyElements();

	// create local vector
	std::shared_ptr<std::vector<double> > localSol =
		std::make_shared<std::vector<double> >(numMyElements, 0.0);

	// extract assembly into local vector
	localB_->ExtractCopy(&(*localSol)[0], numMyElements);

	// solve using serial model
	atmos_->solve(localSol);

	// obtain solution
	localSol = atmos_->getSolution('C');

	// obtain view
	double *sol_tmp;
	localSol_->ExtractView(&sol_tmp);

	// fill view
	for (int i = 0; i != numMyElements; ++i)
	{
		sol_tmp[i] = (*localSol)[i];
	}

	// obtain the solution for the standard map
	domain_->Assembly2Solve(*localSol_, *sol_);
}

//==================================================================
// Very similar to the THCM function
// Create a graph to initialize the Jacobian
void AtmospherePar::createMatrixGraph()
{
	// We know from the discretization that we have at
	// most 5 dependencies in each row. --> This will change, obviously, when extending
	// the atmosphere model. If the number of dependencies varies greatly per row
	// we need to supply them differently (variable). 
	int maxDeps = 5;
	matrixGraph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy, *standardMap_, maxDeps, false));

	// Here we start specifying the indices
	int indices[maxDeps];

	// Get global domain size
	int N = domain_->GlobalN();
	int M = domain_->GlobalM();
	int L = domain_->GlobalL();

	// Get our local range in all directions
	// 0-based
	int I0 = domain_->FirstRealI();
	int J0 = domain_->FirstRealJ();
	int K0 = domain_->FirstRealK();
	int I1 = domain_->LastRealI();
	int J1 = domain_->LastRealJ();
	int K1 = domain_->LastRealK();

	int pos; // position in indices array, not really useful here
	for (int k = K0; k <= K1; ++k)								 
		for (int j = J0; j <= J1; ++j)
			for (int i = I0; i <= I1; ++i)
			{
				// Obtain row corresponding to i,j,k,TT, using 0-based find_row
				int gidU = FIND_ROW_ATMOS0(ATMOS_NUN_, N, M, L, i, j, k, ATMOS_TT_);
				int gid0 = gidU - 1; // used as offset

				pos = 0;

				// Specify dependencies, see Atmosphere::discretize()
				// ATMOS_TT_: 5-point stencil
				insert_graph_entry(indices, pos, i, j, k, ATMOS_TT_, N, M, L);
				insert_graph_entry(indices, pos, i-1, j, k, ATMOS_TT_, N, M, L);
				insert_graph_entry(indices, pos, i+1, j, k, ATMOS_TT_, N, M, L);
				insert_graph_entry(indices, pos, i, j-1, k, ATMOS_TT_, N, M, L);
				insert_graph_entry(indices, pos, i, j+1, k, ATMOS_TT_, N, M, L);

				// Insert dependencies in matrixGraph
				CHECK_ZERO(matrixGraph_->InsertGlobalIndices(gid0+ATMOS_TT_, pos, indices));
			}
	
	// Finalize matrixgraph
	CHECK_ZERO(matrixGraph_->FillComplete());
}

//=============================================================================
// Copied from THCM, adjusted for AtmospherePar
void AtmospherePar::insert_graph_entry(int* indices, int& pos,
									   int i, int j, int k, int xx,
									   int N, int M, int L) const
{
	int ii = i; // if x-boundary is periodic i may be out of bounds.
	// ii will be adjusted in that case:
	if (domain_->IsPeriodic())
    {
		ii = MOD((double)i, (double)N);
    }
	if ((ii>=0) && (j>=0) && (k>=0) &&
		(ii< N) && (j< M) && (k< L) )
	{
		// find index using 0-based find_row
        indices[pos++] = FIND_ROW_ATMOS0(ATMOS_NUN_, N, M, L, ii, j, k, xx);
	}
}
