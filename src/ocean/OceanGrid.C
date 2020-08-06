/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#include "THCMdefs.H"
#include "TRIOS_Domain.H"

#include "OceanGrid.H"

#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"


#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

extern "C" {
// these are defined in usrc.F90 and transform data from a 1D array to 6 3D arrays
    _MODULE_SUBROUTINE_(m_thcm_utils,usol1d)(double *data,
                                             double* u, double* v, double* w,
                                             double* p, double* T, double* S);

// get a copy of the local landm array
    _MODULE_SUBROUTINE_(m_thcm_utils,get_landm)(MaskType* landm);

// z-integration of u-velocity
    _MODULE_SUBROUTINE_(m_thcm_utils,depth_int_u)(double* u, double* us);

// TODO: this should not be used!!!
    _SUBROUTINE_(solu)(double *data, double* u, double* v, double* w, double* p, double* T, double* S);
    _MODULE_SUBROUTINE_(m_thcm_utils,compute_psim)(double *vs, double *psim);

// get grid arrays
    _MODULE_SUBROUTINE_(m_usr,get_grid_data)(double* x, double* y, double* z,
                                             double* xu,double* yv,double* zw);

}

OceanGrid::OceanGrid(Teuchos::RCP<TRIOS::Domain> dom)
{
    DEBUG("Enter OceanGrid constructor...");
    domain = dom;
    n = domain->LocalN();
    m = domain->LocalM();
    l = domain->LocalL();

    dx = (domain->Xmax()-domain->Xmin())/(domain->GlobalN());
    dy = (domain->Ymax()-domain->Ymin())/(domain->GlobalM());
    dz = (domain->Zmax()-domain->Zmin())/(domain->GlobalL());

    importVector = Teuchos::rcp(new Epetra_Vector(*(domain->GetAssemblyMap())));
    // allocate memory and set pointers
    U_=allocateGrid(0,n,0,m,0,l+1);
    V_=allocateGrid(0,n,0,m,0,l+1);
    W_=allocateGrid(0,n+1,0,m+1,0,l);
    P_=allocateGrid(0,n+1,0,m+1,0,l+1);
    T_=allocateGrid(0,n+1,0,m+1,0,l+1);
    S_=allocateGrid(0,n+1,0,m+1,0,l+1);
    Rho_=allocateGrid(0,n+1,0,m+1,0,l+1);

    LandMask_ = new MaskType[(m+2)*(n+2)*(l+2)];

    F90NAME(m_thcm_utils,get_landm)(LandMask_);

    PsiM_ = new double[(m+1)*(l+1)];
    for (int i=0; i<(m+1)*(l+1);i++) PsiM_[i]=0.0;
    PsiB_ = new double[(n+1)*(m+1)];
    for (int i=0; i<(m+1)*(n+1);i++) PsiB_[i]=0.0;

    recompute_PsiM_   = true;
    recompute_PsiB_   = true;
    recompute_MaxVel_ = true;

    // get a communicator with all processes in my 'row' of the
    // processor grid.

    // communication in x-direction:
    DEBUG("OceanGrid: get x-comm");
    xComm = domain->GetProcRow(0);

    // communication in y-direction:
    DEBUG("OceanGrid: get y-comm");
    yComm = domain->GetProcRow(1);

    // create maps for the coordinate vectors
    // we construct them such that the distribution is the same
    // as for our grid cells
    DEBUG("Create x/y/z maps...");
    int nglob = domain->GlobalN();
    xMap = Teuchos::rcp(new Epetra_Map(nglob,0,*xComm));
    int mglob = domain->GlobalM();
    yMap = Teuchos::rcp(new Epetra_Map(mglob,0,*yComm));
    // since the z-direction is not split up we use a replicated map
    zMap = Teuchos::rcp(new Epetra_LocalMap(l,0,*(domain->GetComm())));

    DEBVAR(*xMap);
    DEBVAR(*yMap);
    DEBVAR(*zMap);

    xc_ = Teuchos::rcp(new Epetra_Vector(*xMap));
    xc_->SetLabel("x-coordinate of cell centers");
    yc_ = Teuchos::rcp(new Epetra_Vector(*yMap));
    yc_->SetLabel("y-coordinate of cell centers");
    zc_ = Teuchos::rcp(new Epetra_Vector(*zMap));
    zc_->SetLabel("z-coordinate of cell centers");

    // create xMap0, yMap0.
    // these maps are for the _nodes_, which gives an overlapping decomposition
    int n0 = xMap->NumMyElements()+1;
    int *my_inds = new int[n0];
    for (int i=1;i<n0;i++) my_inds[i]=xMap->GID(i-1)+1;
    my_inds[0] = xMap->MinMyGID();
    int numEl=-1;
    xMap0=Teuchos::rcp(new Epetra_Map(numEl,n0,my_inds,0,*xComm));
    delete [] my_inds;

    int m0 = yMap->NumMyElements()+1;
    my_inds = new int[m0];
    for (int i=1;i<m0;i++) my_inds[i]=yMap->GID(i-1)+1;
    my_inds[0] = yMap->MinMyGID();
    numEl=-1;
    yMap0=Teuchos::rcp(new Epetra_Map(numEl,m0,my_inds,0,*yComm));
    delete [] my_inds;
    int l0 = l+1;
    zMap0 = Teuchos::rcp(new Epetra_LocalMap(l0,0,*(domain->GetComm())));

    xu_ = Teuchos::rcp(new Epetra_Vector(*xMap0));
    xu_->SetLabel("x-coordinate of grid nodes");
    yv_ = Teuchos::rcp(new Epetra_Vector(*yMap0));
    yv_->SetLabel("y-coordinate of grid nodes");
    zw_ = Teuchos::rcp(new Epetra_Vector(*zMap0));
    zw_->SetLabel("z-coordinate of grid nodes");

    // get arrays from fortran. These contain overlap, so we can't
    // put them directly into the newly created vectors.
    double *x, *y, *z, *xu, *yv, *zw;
    x = new double[n];   y = new double[m];   z = new double[l];
    xu= new double[n+1]; yv= new double[m+1]; zw= new double[l+1];

    DEBUG("get grid arrays from Fortran...");

    F90NAME(m_usr,get_grid_data)(x,y,z,xu,yv,zw);

    int i0 = domain->FirstRealI()-domain->FirstI();
    int j0 = domain->FirstRealJ()-domain->FirstJ();
    int k0 = domain->FirstRealK()-domain->FirstK();

    for (int i=0;i<n0-1;i++) (*xc_)[i] = x[i0+i];
    for (int j=0;j<m0-1;j++) (*yc_)[j] = y[j0+j];
    for (int k=0;k<l0-1;k++) (*zc_)[k] = z[k0+k];

    // note: there are either two or no ghost cells,
    //       but there is always a grid node 0 (even
    //       at global domain bounds). If there is a
    //       ghost layer we need to pick index 1 as
    //       first node.
    i0 = std::max(0,i0-1);
    j0 = std::max(0,j0-1);
    k0 = std::max(0,k0-1);

    for (int i=0;i<n0;i++) (*xu_)[i] = xu[i0+i];
    for (int j=0;j<m0;j++) (*yv_)[j] = yv[j0+j];
    for (int k=0;k<l0;k++) (*zw_)[k] = zw[k0+k];


    delete [] x; delete [] xu;
    delete [] y; delete [] yv;
    delete [] z; delete [] zw;
    DEBUG("OceanGrid constructor done");
}

OceanGrid::~OceanGrid()
{
    INFO("Destroy OceanGrid...");
    delete [] U_;
    delete [] V_;
    delete [] W_;
    delete [] P_;
    delete [] T_;
    delete [] S_;
    delete [] Rho_;
    delete [] PsiM_;
    delete [] PsiB_;
    delete [] LandMask_;
}



// put data from vector in grid format.
//! input should be based on the 'Solve' map
void OceanGrid::ImportData(const Epetra_Vector& input)
{
    DEBUG("Import vector into grid...");
    // import boundaries into ghost-nodes
    DEBUG("(a) get ghost-nodes...");
    domain->Solve2Assembly(input,*importVector);
    double *data;
    DEBUG("(b) extract pointer...");
    CHECK_ZERO(importVector->ExtractView(&data));

    // puts the data into the locations (1:n,1:m,1:l) and
    // adds some dummy-data and boundary conditions in the boundary cells
    DEBUG("(c) call Fortran routine...");
    F90NAME(m_thcm_utils,usol1d)(data,U_,V_,W_,P_,T_,S_);

    recompute_PsiM_ = true;
    recompute_PsiB_ = true;
    recompute_MaxVel_ = true;

    DEBUG("done!");
//    DEBVAR(input);
//    DEBVAR(*this);
}

void OceanGrid::ImportLandMask(const Epetra_IntVector& mask)
{
    // we require the input to be distributed already, including
    // overlap and all, so we can simply extract the data
    CHECK_ZERO(mask.ExtractCopy(LandMask_));
}

//! put data from grid in a vector.
//! input should be based on the 'Solve' map
void OceanGrid::ExportData(Epetra_Vector& output)
{
    double *data;
    CHECK_ZERO(importVector->ExtractView(&data));
    // puts the data from the locations (1:n,1:m,1:l)
    // into the assembly vector
    //TODO: this should not be used, first make a C-compatible '1D' version!
    ERROR("Not Implemented!",__FILE__,__LINE__);
    FNAME(solu)(data,U_,V_,W_,P_,T_,S_);
    // import boundaries into ghost-nodes
    domain->Assembly2Solve(*importVector,output);
}

//  private:

// memory allocation function. This just
// allocates a block of memory large enough
// to hold the required data. Indexing has to be
// done by the access operators u() etc.
double* OceanGrid::allocateGrid(int imin, int imax,
                                int jmin, int jmax,
                                int kmin, int kmax)
{
    DEBUG("Allocate Grid-Array: ["<<imin<<".."<<imax<<"]x["<<jmin<<".."<<jmax<<"]x["<<kmin<<".."<<kmax<<"]");
    //total memory required
    int dim = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);
    //allocate big block
    double *G = new double[dim];
    if (G==NULL)
    {
        INFO("failed to allocate "<<dim*8*1e-6<<" MB of memory for a grid function!");
        ERROR("Memory allocation error!",__FILE__,__LINE__);
    }
    for (int i=0;i<dim;i++) G[i]=0.0;
    return G;
}




void OceanGrid::recomputePsiM()
{
    double sum;

    DEBUG("OceanGrid: compute Psi_m");

    // integrate y-velocity v in x-direction,
    // resulting in an average v-field in the y-z plane

    // vs(0:m,1:l)
    double *vs = new double[(m+1)*l];

// k is 1-based in this indexing function for vs:
#define IND_VS(j,k) (m+1)*((k)-1)+(j)

    // perform local integration
    // we perform the local integration into PsiM, which
    // is abused as a buffer here (PsiM is bigger than vs, so we
    // can do this). We omit ghost-cells to get the global
    // integral correct:

    // note that the grid requires 1-based indexing
    int nghostleft = domain->FirstRealI()-domain->FirstI();
    int nghostright = domain->LastI()-domain->LastRealI();
    int imin = 1+nghostleft;
    int imax = n - nghostright;
    DEBVAR(nghostleft);
    DEBVAR(nghostright);
    DEBVAR(imin);
    DEBVAR(imax);
    for (int k=1;k<=l;k++)
        for (int j=0;j<=m;j++)
        {
            sum = 0.0;
            for (int i=imin;i<=imax;i++)
            {
                sum += v(i,j,k);
            }
            PsiM_[IND_VS(j,k)] = sum*dx;
        }
    // perform an allreduce operation over all processes in the
    // x-direction, the resulting array is the global integral
    // this is allowed because all subdomains in the x-direction
    // have the same (m,l). The values in the ghost-nodes have
    // been computed as well (at j=0 and m).
    CHECK_ZERO(xComm->SumAll(PsiM_,vs,(m+1)*l));

    // clean away the rubbish we dumped in this aray:
    for (int i=0;i<(m+1)*l;i++) PsiM_[i]=0.0;

    // now we can compute the stream-function Psi

    // this is sequential as the z-direction is not split up:
    F90NAME(m_thcm_utils,compute_psim)(vs, PsiM_);

    // compute the local maximum and minimum
    double PsiMinL=psiM(0,0);
    double PsiMaxL = PsiMinL;
    for (int j=0; j<=m; j++)
        for (int k=0; k<=l; k++)
        {
            PsiMaxL = std::max(psiM(j,k),PsiMaxL);
            PsiMinL = std::min(psiM(j,k),PsiMinL);
        }
    // global maximum
    Teuchos::RCP<Epetra_Comm> comm = domain->GetComm();
    CHECK_ZERO(comm->MaxAll(&PsiMaxL,&PsimMax_,1));
    CHECK_ZERO(comm->MinAll(&PsiMinL,&PsimMin_,1));

    DEBVAR(PsiMinL);
    DEBVAR(PsiMaxL);
    DEBVAR(PsimMin_);
    DEBVAR(PsimMax_);

    delete [] vs;
    recompute_PsiM_=false;
}

void OceanGrid::recomputePsiB()
{
    DEBUG("OceanGrid: compute Psi_b");

    // integrate x-velocity u in z-direction,
    // resulting in an average u-field in the x-y plane

    // us(0:n,0:m)
    double *us = new double[(m+1)*(n+1)];

#define IND_US(i,j) (n+1)*(j)+(i)

    // this is parallel since the z-direction is not split up:
    DEBUG("depth-integrate U for computation of Psi_B");
    F90NAME(m_thcm_utils,depth_int_u)(U_,us);

    // number of ghost nodes at top and bottom of domain
    int nghosttop = domain->LastJ()-domain->LastRealJ();
    int nghostbottom = domain->FirstRealJ()-domain->FirstJ();

    int imin = 0;
    int imax = n;
    int jmin = 1+nghostbottom;
    int jmax = m-nghosttop;

    // now we can compute the stream-function Psi

    // first we compute the barotropic streamfunction (\int_y us dy)
    // on the subdomain, then we perform a scan-sum to make it global.

    DEBUG("compute local Psi_B");


    // this is inspired by THCM file mstream.f, subroutine bstream
    for (int i=imin; i<=imax; i++) psiB(i,jmin-1) = 0.0;

    for (int i=imin; i<=imax; i++)
        for (int j=jmin; j<=jmax; j++)
        {
            psiB(i,j)=0.5*(us[IND_US(i,j-1)]+us[IND_US(i,j)])*dy + psiB(i,j-1);
        }

    // now, on a process in processor row J, we have the integral from
    // j0 to j in local row j, j0 may be a ghost node. We need to add
    // to all psi-values in the subdomain the row j0-1 from proc p-1.

#ifdef HAVE_MPI
    DEBUG("Psi_B: MPI part");
    // The Comm->ScanSum function can't be used here, we do it ourselves:
    Teuchos::RCP<Epetra_MpiComm> mpi = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(yComm);
    DEBVAR(*mpi);
    if (mpi!=Teuchos::null) // might be a serial comm, I guess
    {
        MPI_Comm ycomm = mpi->Comm(); // get actual MPI communicator

        // get info on the processor column
        int yprocs,yrank;
        MPI_Comm_rank(ycomm,&yrank);
        MPI_Comm_size(ycomm,&yprocs);

        int count = n+1;
        double *buf = new double[count];
        MPI_Datatype type = MPI_DOUBLE;
        int tag = 0;
        if (yrank==0)
        {
            if (yprocs>1)
            {
                for (int i=0;i<=n;i++)
                {
                    buf[i] = psiB(i,jmax);
                }
                MPI_Send(buf,count,type,1,tag,ycomm);
            }
        }
        else
        {
            // we sequentially send data to the next proc in the column
            // this could be done more efficiently but this shouldn't be
            // a crucial routine when it comes to time.
            for (int proc=1; proc<yprocs; proc++)
            {
                if (yrank==proc)
                {
                    MPI_Recv(buf,count,type,yrank-1,tag,ycomm,MPI_STATUS_IGNORE);
                    for (int j=jmin; j<=jmax; j++)
                        for (int i=imin;i<=imax;i++)
                        {
                            psiB(i,j) += buf[i];
                        }
                    if (yrank<yprocs-1)
                    {
                        for (int i=imin;i<=imax;i++)
                        {
                            buf[i] = psiB(i,jmax);
                        }
                        MPI_Send(buf,count,type,yrank+1,tag,ycomm);
                    }
                }
            }// for all procs
        }// not the root process
        delete [] buf;
    }// not an MPI communicator
#endif

    delete [] us;
#ifdef DEBUGGING
    std::cout << "\nPsiM(0:m,0:l): "<<std::endl;
    for (int k=l; k>=0; k--)
    {
        for (int j=0; j<=m; j++)
        {
            std::cout << psiM(j,k) << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << "\nPsiB(0:n,0:l): "<<std::endl;
    for (int j=m; j>=0; j--)
    {
        for (int i=0; i<=n; i++)
        {
            std::cout << psiB(i,j) << "\t";
        }
        std::cout << std::endl;
    }
#endif

    // compute the local maximum and minimum
    DEBUG("Compute max/min of Psi_B");
    double PsiMinL=psiB(0,0);
    double PsiMaxL = PsiMinL;
    for (int j=jmin; j<=jmax; j++)
        for (int i=imin; i<=imax; i++)
        {
            PsiMaxL = std::max(psiB(i,j),PsiMaxL);
            PsiMinL = std::min(psiB(i,j),PsiMinL);
        }
    // global maximum
    Teuchos::RCP<Epetra_Comm> comm = domain->GetComm();
    CHECK_ZERO(comm->MaxAll(&PsiMaxL,&PsibMax_,1));
    CHECK_ZERO(comm->MinAll(&PsiMinL,&PsibMin_,1));

    //DEBVAR(PsiMinL);
    //DEBVAR(PsiMaxL);
    //DEBVAR(PsibMin_);
    //DEBVAR(PsibMax_);

    recompute_PsiB_=false;
}

void OceanGrid::recomputeMaxVel(void)
{
    // local maxima
    double maxU=0;
    double maxV=0;
    double maxW=0;
    // compute local max of u,v,w
    for (int k=1;k<=l;k++)
        for (int j=1;j<=m;j++)
            for (int i=1;i<=n;i++)
            {
                maxU=std::max(maxU,std::abs((u(i,j,k)+u(i-1,j,k)+u(i-1,j-1,k)+u(i,j-1,k))*0.25));
                maxV=std::max(maxV,std::abs((v(i,j,k)+v(i-1,j,k)+v(i-1,j-1,k)+v(i,j-1,k))*0.25));
                maxW=std::max(maxW,std::abs((w(i,j,k)+w(i,j,k-1))*0.5));
            }
    // compute global maxima
    CHECK_ZERO(domain->GetComm()->MaxAll(&maxU,&MaxU,1));
    CHECK_ZERO(domain->GetComm()->MaxAll(&maxV,&MaxV,1));
    CHECK_ZERO(domain->GetComm()->MaxAll(&maxW,&MaxW,1));
    DEBUG("Maximum velocities: ");
    DEBVAR(MaxU);
    DEBVAR(MaxV);
    DEBVAR(MaxW);
    recompute_MaxVel_=false;
}


// output to file stream
std::ostream& OceanGrid::print(std::ostream& os) const
{
    os << "OceanGrid"<<std::endl;
    os << "n = "<<n<<std::endl;
    os << "m = "<<m<<std::endl;
    os << "l = "<<l<<std::endl;
    os << std::endl;
    os << *xc_<<std::endl;
    os << *yc_<<std::endl;
    os << *zc_<<std::endl;
    os << *xu_<<std::endl;
    os << *yv_<<std::endl;
    os << *zw_<<std::endl;
    os << "LandMask(0:n+1,0:m+1,0:l+1): "<<std::endl;
    for (int k=0;k<=l+1;k++)
    {
        os << "k="<<k<<std::endl;
        for (int j=m+1; j>=0; j--)
        {
            for (int i=0; i<=n+1; i++)
            {
                os << landm(i,j,k);
            }
            os << std::endl;
        }
        os << std::endl;
    }
    os << "U(0:n,0:m,0:l+1): "<<std::endl;
    for (int k=0;k<=l+1;k++)
    {
        os << "k="<<k<<std::endl;
        for (int j=m; j>=0; j--)
        {
            for (int i=0; i<=n; i++)
            {
                os << u(i,j,k) << "\t";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    os << "\nV(0:n,0:m,0:l+1): "<<std::endl;
    for (int k=0;k<=l+1;k++)
    {
        os << "k="<<k<<std::endl;
        for (int j=m; j>=0; j--)
        {
            for (int i=0; i<=n; i++)
            {
                os << v(i,j,k) << "\t";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    os << "\nW(0:n,0:m,0:l): "<<std::endl;
    for (int k=0;k<=l;k++)
    {
        os << "k="<<k<<std::endl;
        for (int j=m; j>=0; j--)
        {
            for (int i=0; i<=n; i++)
            {
                os << w(i,j,k) << "\t";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    os << "\np(0:n+1,0:m+1,0:l+1): "<<std::endl;
    for (int k=0;k<=l+1;k++)
    {
        os << "k="<<k<<std::endl;
        for (int j=m+1; j>=0; j--)
        {
            for (int i=0; i<=n+1; i++)
            {
                os << p(i,j,k) << "\t";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    os << "\nT(0:n+1,0:m+1,0:l+1): "<<std::endl;
    for (int k=0;k<=l+1;k++)
    {
        os << "k="<<k<<std::endl;
        for (int j=m+1; j>=0; j--)
        {
            for (int i=0; i<=n+1; i++)
            {
                os << T(i,j,k) << "\t";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    os << "\nS(0:n+1,0:m+1,0:l+1): "<<std::endl;
    for (int k=0;k<=l+1;k++)
    {
        os << "k="<<k<<std::endl;
        for (int j=m+1; j>=0; j--)
        {
            for (int i=0; i<=n+1; i++)
            {
                os << S(i,j,k) << "\t";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    if (std::abs(PsimMax_)>0)
    {
        os << "\nPsiM(0:m,0:l): "<<std::endl;
        for (int k=l; k>=0; k--)
        {
            for (int j=0; j<=m; j++)
            {
                os << psiM(j,k) << "\t";
            }
            os << std::endl;
        }
    }
    if (std::abs(PsibMax_)>0)
    {
        os << "\nPsiB(0:n,0:l): "<<std::endl;
        for (int j=m; j>=0; j--)
        {
            for (int i=0; i<=n; i++)
            {
                os << psiB(i,j) << "\t";
            }
            os << std::endl;
        }
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const OceanGrid& G)
{
    return G.print(os);
}


