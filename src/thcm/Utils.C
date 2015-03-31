#include "Utils.H"
//========================================================================================
int Utils::SplitBox(int nx, int ny, int nz,
					int nparts, int& ndx, int& ndy, int& ndz,
					int sx, int sy, int sz)
{
	// Factor the number of processors into two dimensions. (nprocs = npN*npM)
	double rmin = 1e100;
	int ret = 1;

	int npx = nx / sx;
	int npy = ny / sy;
	int npz = nz / sz;

	std::string s1 = Teuchos::toString(nx) + "x" +
		Teuchos::toString(ny) + "x" + Teuchos::toString(nz);
	std::string s2 = Teuchos::toString(npx) + "x" +
		Teuchos::toString(npy) + "x" + Teuchos::toString(npz);

	// check all possibilities:
	for (int t1 = 1; t1 <= nparts; t1++)
		for (int t2 = 1; t2 <= (int)(nparts / t1); t2++)
		{
			int t3 = (int)(nparts/(t1*t2));

			if (t1 * t2 * t3 == nparts)
			{
				std::string s3 = Teuchos::toString(t1) + "x" +
					Teuchos::toString(t2) + "x" + Teuchos::toString(t3);
				int my_nx = nx / t1;
				int my_ny = ny / t2;
				int my_nz = nz / t3;
				if ((my_nx * t1 != nx) || (my_ny * t2 != ny) || (my_nz * t3 != nz))
				{
					DEBUG("Can't partition a "+s1+" domain into "+s3+" parts.");
					continue;
				}
				int my_npx = npx / t1;
				int my_npy = npy / t2;
				int my_npz = npz / t3;
				if ((my_npx * sx != my_nx) ||
					(my_npy * sy != my_ny) || (my_npz * sz != my_nz))
				{
					DEBUG("Can't partition "+s2+" domains onto "+s3+" processors.");
					continue;
				}

				double r1 = std::abs((double)nx / (double)t1 - (double)ny / (double)t2);
				double r2 = std::abs((double)nx / (double)t1 - (double)nz / (double)t3);
				double r3 = std::abs((double)ny / (double)t2 - (double)nz / (double)t3);
				double r = r1 + r2 + r3;

				if (r < rmin)
				{
					rmin = r;
					ndx  = t1;
					ndy  = t2;
					ndz  = t3;
					ret  = 0;
				}
			}
		}
	return ret;
}
//========================================================================================
void Utils::ind2sub(int nx, int ny, int nz, int dof, 
					int idx, int& i, int& j, int& k, int& var)
{
#ifdef TESTING    
	if ((idx < 0) || (idx >= (nx * ny * nz * dof)))
	{
		std::cerr << "dim=[" << nx << "," << ny << ","
				  << nz << "], dof=" << dof << ": ind=" <<idx << std::endl;
		Error("ind2sub: Index out of range!",__FILE__,__LINE__);
	}
#endif      
	int rem = idx;
	var = MOD(rem, dof);
	rem = (rem - var) / dof;
	i   = MOD(rem , nx);
	rem = (rem - i) / nx;
	j   = MOD(rem, ny);
	rem = (rem - j) / ny;
	k   = MOD(rem, nz);
}
//========================================================================================
void Utils::ind2sub(int nx, int ny, int nz, int idx, int& i, int& j, int& k)
{
	int dummy;
	ind2sub(nx,ny,nz,1,idx,i,j,k,dummy);
	return;
}
//========================================================================================
//! converts cartesian subscripts to linear index
int Utils::sub2ind(int nx, int ny, int nz, int dof, int i, int j, int k, int var)
{
#ifdef TESTING    
	std::string msg1 = "sub2ind: ";
	std::string msg3 = " out of range ";
	if ((i < 0) || (i >= nx))
	{
		std::string msg2 = "i-Index "+Teuchos::toString(i);
		std::string msg4 = "[0,"+Teuchos::toString(nx)+"]";
		ERROR(msg1+msg2+msg3+msg4,__FILE__,__LINE__);
	}
	if ((j < 0) || ( j >= ny))
	{
		std::string msg2 = "j-Index "+Teuchos::toString(j);
		std::string msg4 = "[0,"+Teuchos::toString(ny)+"]";
		ERROR(msg1+msg2+msg3+msg4,__FILE__,__LINE__);
	}
	if ((k < 0) || (k >= nz))
	{
		std::string msg2 = "k-Index "+Teuchos::toString(j);
		std::string msg4 = "[0,"+Teuchos::toString(nz)+"]";
		ERROR(msg1+msg2+msg3+msg4,__FILE__,__LINE__);
	}
	if ((var<0)||(var>=dof))
	{
		std::string msg2 = "var-Index "+Teuchos::toString(j);
		std::string msg4 = "[0,"+Teuchos::toString(dof)+"]";
		ERROR(msg1+msg2+msg3+msg4,__FILE__,__LINE__);
	}
#endif      
	return ((k*ny+j)*nx+i)*dof+var;
}

//========================================================================================
Teuchos::RCP<Epetra_Map> Utils::CreateMap(int nx, int ny, int nz,
										  int dof, int indexbase,
										  const Epetra_Comm& comm,
										  int numActiveProcs)
{
	if (numActiveProcs == -1) numActiveProcs = comm.NumProc();
	// create a parallel map. We first figure out where in the domain we are   
	int np  = std::min(comm.NumProc(),nx*ny*nz);
	np      = std::min(np,numActiveProcs);
	int pid = comm.MyPID();

	int npX, npY, npZ;
	int pidX  = -1, pidY = -1, pidZ = -1;
	int offX  = -1, offY = -1, offZ = -1;
	int nXloc = 0, nYloc = 0, nZloc = 0;
	
	bool active = (pid < np);
	
	Utils::SplitBox(nx, ny, nz, np, npX, npY, npZ);
	
	if (active)
    {
		Utils::ind2sub(npX, npY, npZ, pid, pidX, pidY, pidZ);
		// dimension of subdomain
		nXloc = (int) (nx / npX);
		nYloc = (int) (ny / npY);
		nZloc = (int) (nz / npZ);

		// offsets for local->global index conversion
		offX = pidX * (int) (nx / npX);
		offY = pidY * (int) (ny / npY);
		offZ = pidZ * (int) (nz / npZ);
    
		// distribute remaining points among first few cpu's
		int remX = nx % npX;
		int remY = ny % npY;
		int remZ = nz % npZ;

		if (pidX < remX) nXloc++;
		if (pidY < remY) nYloc++;
		if (pidZ < remZ) nZloc++;

		for (int i=0;i<std::min(remX,pidX);i++) offX++;
		for (int i=0;i<std::min(remY,pidY);i++) offY++;
		for (int i=0;i<std::min(remZ,pidZ);i++) offZ++;
    } //active?
	if (indexbase != 0)
    {
		// this could easily be implemented but I didn't need it up to now.
		ERROR("only index base 0 is implemented right now",
			  __FILE__, __LINE__);
    }

	return CreateMap(offX, offX + nXloc - 1,                                  
					 offY, offY + nYloc - 1,
					 offZ, offZ + nZloc - 1,
					 0, nx-1, 0, ny-1, 0, nz-1,
					 dof, comm);
      
}
//========================================================================================
Teuchos::RCP<Epetra_Map> Utils::CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,
										  int I0, int I1, int J0, int J1, int K0, int K1,
										  int dof,
										  const Epetra_Comm& comm)
{
	Teuchos::RCP<Epetra_Map> result = Teuchos::null;
	DEBUG("Utils::CreateMap ");
	DEBUG("["<<i0<<".."<<i1<<"]");
	DEBUG("["<<j0<<".."<<j1<<"]");
	DEBUG("["<<k0<<".."<<k1<<"]");
		  
	int n = std::max(i1-i0+1,0); int N=I1-I0+1;
	int m = std::max(j1-j0+1,0); int M=J1-J0+1;
	int l = std::max(k1-k0+1,0); int L=K1-K0+1;
		  
	DEBVAR(N);
	DEBVAR(M);
	DEBVAR(L);
		  
	int NumMyElements = n*m*l*dof;
	int NumGlobalElements = -1; // note that there may be overlap
	int *MyGlobalElements = new int[NumMyElements];
		  
	int pos = 0;
	for (int k=k0; k<=k1; k++)
		for (int j=j0; j<=j1; j++)
			for (int i=i0; i<=i1; i++)
				for (int var=0;var<dof;var++)
				{
					MyGlobalElements[pos++] = 
						Utils::sub2ind(N,M,L,dof,i,j,k,var);
				}
	result = Teuchos::rcp(new Epetra_Map(NumGlobalElements,
										 NumMyElements,MyGlobalElements,0,comm));
	delete [] MyGlobalElements;
	return result;
}
//========================================================================================
Teuchos::RCP<Epetra_Map> Utils::CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,        
												int I0, int I1, int J0, int J1, int K0, int K1,
												const Epetra_Comm& comm)
{
	Teuchos::RCP<Epetra_Map> result = Teuchos::null;
    
	DEBUG("Utils::CreateMap ");
	DEBUG("["<<i0<<".."<<i1<<"]");
	DEBUG("["<<j0<<".."<<j1<<"]");
	DEBUG("["<<k0<<".."<<k1<<"]");
      
	int n = i1-i0+1; int N=I1-I0+1;
	int m = j1-j0+1; int M=J1-J0+1;
	int l = k1-k0+1; int L=K1-K0+1;
      
	DEBVAR(M);
	DEBVAR(N);
	DEBVAR(L);
      
	int NumMyElements = n*m*l;
	int NumGlobalElements = -1; // note that there may be overlap
	int *MyGlobalElements = new int[NumMyElements];
      
	int pos = 0;
	for (int k=k0; k<=k1; k++)
        for (int j=j0; j<=j1; j++)
			for (int i=i0; i<=i1; i++)
				MyGlobalElements[pos++] = k*N*M + j*N + MOD((double)i,(double)N);

	result = Teuchos::rcp(new Epetra_Map(NumGlobalElements,
										 NumMyElements,MyGlobalElements,0,comm));
	delete [] MyGlobalElements;
	return result;
}

//========================================================================================
Teuchos::RCP<Epetra_CrsMatrix> Utils::Gather(const Epetra_CrsMatrix& mat, int root)
{
	const Epetra_Map& rowmap_dist = mat.RowMap();
	// we take the domain map as the colmap is potentially overlapping
	const Epetra_Map& colmap_dist = mat.DomainMap();
	// gather the row map
	Teuchos::RCP<Epetra_Map> rowmap = Gather(rowmap_dist,root);
	// gather the col map
	Teuchos::RCP<Epetra_Map> colmap = Gather(colmap_dist,root);
    
	// we only guess the number of row entries, this routine is not performance critical
	// as it should only be used for debugging anyway
	int num_entries = mat.NumGlobalNonzeros() / mat.NumGlobalRows();
	Teuchos::RCP<Epetra_CrsMatrix> gmat =
		Teuchos::rcp(new Epetra_CrsMatrix(Copy,*rowmap, *colmap, num_entries) );
      
	Teuchos::RCP<Epetra_Import> import =
		Teuchos::rcp(new Epetra_Import(*rowmap,rowmap_dist) );
      
	CHECK_ZERO(gmat->Import(mat,*import,Insert));
	CHECK_ZERO(gmat->FillComplete());
	gmat->SetLabel(mat.Label());
	return gmat;
}
//========================================================================================
Teuchos::RCP<Epetra_MultiVector> Utils::Gather(const Epetra_MultiVector& vec, int root)
{
	const Epetra_BlockMap& map_dist = vec.Map();
	Teuchos::RCP<Epetra_BlockMap> map = Gather(map_dist,root);
	Teuchos::RCP<Epetra_MultiVector> gvec = 
        Teuchos::rcp(new Epetra_MultiVector(*map,vec.NumVectors()));
	Teuchos::RCP<Epetra_Import> import =
		Teuchos::rcp(new Epetra_Import(*map,map_dist) );
	CHECK_ZERO(gvec->Import(vec,*import,Insert));
	gvec->SetLabel(vec.Label());
	return gvec;      
}
//========================================================================================
// create "Gather" map from "Solve" map
Teuchos::RCP<Epetra_BlockMap> Utils::Gather(const Epetra_BlockMap& map, int root)
{
    int ElementSize = map.ElementSize();
#ifdef TESTING
    if (ElementSize != 1)
	{
		DEBVAR(ElementSize);
		ERROR("this is possibly not implemented correctly!",
			  __FILE__, __LINE__);
		ElementSize=1;
	}
#endif
    int NumMyElements = map.NumMyElements();
    int NumGlobalElements = map.NumGlobalElements();
    const Epetra_Comm& Comm = map.Comm();
    int *MyGlobalElements = new int[NumMyElements];
    int *AllGlobalElements = NULL;	
    for (int i = 0; i < NumMyElements; i++)
	{
		MyGlobalElements[i] = map.GID(i);
	}    
    if (Comm.MyPID() == root)
	{
		AllGlobalElements = new int[NumGlobalElements];
	}
	if (Comm.NumProc()>1)
	{
#ifdef HAVE_MPI		
		const Epetra_MpiComm MpiComm = dynamic_cast<const Epetra_MpiComm&>(Comm);
		int *counts, *disps;
		counts = new int[Comm.NumProc()];
		disps  = new int[Comm.NumProc()+1];
		MPI_Gather(&NumMyElements,1,MPI_INTEGER,
				   counts,1,MPI_INTEGER,root,MpiComm.GetMpiComm());
		
		if (Comm.MyPID()==root)
		{
			disps[0]=0;
			for (int p=0;p<Comm.NumProc();p++)
			{
				disps[p+1] = disps[p]+counts[p];
			}
		}
		
		MPI_Gatherv(MyGlobalElements, NumMyElements,MPI_INTEGER, 
					AllGlobalElements, counts,disps, MPI_INTEGER, root, MpiComm.GetMpiComm());
		delete [] counts;
		delete [] disps;                
#else
		ERROR("No MPI but still parallel??? We don't do that.", __FILE__, __LINE__);
#endif
	}
	else
	{
		for (int i = 0; i < NumMyElements; i++)
			AllGlobalElements[i] = MyGlobalElements[i];
	}  
    if (Comm.MyPID()!=root) 
	{
		NumMyElements=0;
	}
    else
	{
		NumMyElements = NumGlobalElements;
		Teuchos::ArrayView<int> view(AllGlobalElements, NumGlobalElements);
		std::sort(view.begin(), view.end());
		Teuchos::ArrayView<int>::iterator new_end = std::unique(view.begin(), view.end());
		NumMyElements = std::distance(view.begin(), new_end);
		NumGlobalElements = NumMyElements;
	}
	CHECK_ZERO(Comm.Broadcast(&NumGlobalElements, 1, 0));
	// build the new (gathered) map
	Teuchos::RCP<Epetra_BlockMap> gmap =
		Teuchos::rcp(new Epetra_BlockMap(NumGlobalElements, NumMyElements, AllGlobalElements, 
										 ElementSize, map.IndexBase(), Comm) );
	if (Comm.MyPID()==root)
	{      
		delete [] AllGlobalElements;
	}
    delete [] MyGlobalElements;
    return gmap;    
}
//========================================================================================		
// create "Gather" map from "Solve" map
Teuchos::RCP<Epetra_Map> Utils::Gather(const Epetra_Map& map, int root)
{
	int NumMyElements = map.NumMyElements();
	int NumGlobalElements = map.NumGlobalElements();
	const Epetra_Comm& Comm = map.Comm();
    
	int *MyGlobalElements = new int[NumMyElements];
	int *AllGlobalElements = NULL;

	for (int i = 0; i < NumMyElements; i++)
		MyGlobalElements[i] = map.GID(i);
	    
	if (Comm.MyPID() == root)
		AllGlobalElements = new int[NumGlobalElements];

	if (Comm.NumProc()>1)
	{    
#ifdef HAVE_MPI    
		const Epetra_MpiComm MpiComm = dynamic_cast<const Epetra_MpiComm&>(Comm);
		int *counts, *disps;
		counts = new int[Comm.NumProc()];
		disps  = new int[Comm.NumProc()+1];
		MPI_Gather(&NumMyElements, 1, MPI_INTEGER, 
				   counts, 1, MPI_INTEGER, root, MpiComm.GetMpiComm());
		
		if (Comm.MyPID() == root)
		{
			disps[0] = 0;
			for (int p = 0; p < Comm.NumProc(); p++)
				disps[p+1] = disps[p] + counts[p];
		}

		MPI_Gatherv(MyGlobalElements,  NumMyElements,MPI_INTEGER, 
					AllGlobalElements, counts,disps, MPI_INTEGER, root, MpiComm.GetMpiComm());
		delete [] counts;
		delete [] disps;                
#else
		ERROR("No MPI but still parallel??? We don't do that.",__FILE__,__LINE__);
#endif
	}
	else
		for (int i=0; i < NumMyElements; i++)
			AllGlobalElements[i] = MyGlobalElements[i];

	if (Comm.MyPID() != root) 
		NumMyElements = 0;
	else
	{
		NumMyElements=NumGlobalElements;
		std::sort(AllGlobalElements,AllGlobalElements+NumGlobalElements);
	}
	// build the new (gathered) map
	Teuchos::RCP<Epetra_Map> gmap =
		Teuchos::rcp(new Epetra_Map(NumGlobalElements, NumMyElements, AllGlobalElements, 
									map.IndexBase(), Comm));
    
	if (Comm.MyPID() == root)
		delete [] AllGlobalElements;
	
	delete [] MyGlobalElements;
    
	return gmap;    
}
//========================================================================================
// distribute a gathered vector among processors
Teuchos::RCP<Epetra_MultiVector> Utils::Scatter
(const Epetra_MultiVector& vec, const Epetra_BlockMap& distmap)
{
    Teuchos::RCP<Epetra_MultiVector> dist_vec =
		Teuchos::rcp(new Epetra_MultiVector(distmap,vec.NumVectors()));
    Teuchos::RCP<Epetra_Import> import =
		Teuchos::rcp(new Epetra_Import(vec.Map(),distmap));
    CHECK_ZERO(dist_vec->Export(vec,*import,Insert));
    return dist_vec;
}
//========================================================================================
