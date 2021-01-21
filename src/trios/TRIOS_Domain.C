/**********************************************************************
 * Copyright by s Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

#include "TRIOS_Domain.H"

#include <fstream>
#include <vector>

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"

#include "Utils.H"
#include "THCMdefs.H"
#include "GlobalDefinitions.H"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include <mpi.h>
#endif

namespace TRIOS
{
    Grid::Grid() : Grid(0.0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
    {}

    Grid::Grid(double qz, int N, int M, int L, double xmin, double xmax,
               double ymin, double ymax, double hdim)
        : n_(N), m_(M), l_(L)
        , xmin_(xmin), xmax_(xmax)
        , ymin_(ymin), ymax_(ymax)
        , zmin_(-1.0), zmax_(0), hdim_(hdim)
        , x_(N), y_(M), z_(L), xu_(N+1), yv_(M+1), zw_(L+1)
        , qz_(qz)
        , depth_(L+1, 0.0), cumulativeDepth_(L, 0.0)
    {
        double dx = (xmax_ - xmin_) / n_;
        double dy = (ymax_ - ymin_) / m_;
        double dz = (zmax_ - zmin_) / l_;

        xu_[0] = xmin_;
        yv_[0] = ymin_;
        zw_[0] = zmin_;

        for (int i = 0; i != n_; ++i)
        {
            x_[i]    = xmin_ + (i + 0.5) * dx;
            xu_[i+1] = xmin_ + (i + 1.0) * dx;
        }

        for (int j = 0; j != m_; ++j)
        {
            y_[j]    = ymin_ + (j + 0.5) * dy;
            yv_[j+1] = ymin_ + (j + 1.0) * dy;
        }

        for (int k = 0; k != l_; ++k)
        {
            z_[k]    = fz(zmin_ + (k + 0.5) * dz);
            zw_[k+1] = fz(zmin_ + (k + 1.0) * dz);
            depth_[k] = dfdz(zmin_ + (k + 0.5) * dz);

            for (int i = 0; i <= k; i++)
            {
                cumulativeDepth_[i] += depth_[k];
            }
        }
    }

    Grid Grid::SubGrid(int Noff, int Moff, int nloc, int mloc, int lloc) const
    {
        double dx = (xmax_ - xmin_) / n_;
        double dy = (ymax_ - ymin_) / m_;

        double xminLoc = xmin_ + Noff*dx;
        double xmaxLoc = xmin_ + (Noff+nloc)*dx;
        double yminLoc = ymin_ + Moff*dy;
        double ymaxLoc = ymin_ + (Moff+mloc)*dy;

        return Grid(qz_, nloc, mloc, lloc, xminLoc, xmaxLoc, yminLoc, ymaxLoc, hdim_);
    }

    double Grid::getXposEdge(int x) const
    {
        if (x < 0 || x > n_) {
            ERROR("X index out of bounds!",__FILE__,__LINE__);
        }

        return xu_[x];
    }

    double Grid::getXposCenter(int x) const
    {
        if (x < 0 || x >= n_) {
            ERROR("X index out of bounds!",__FILE__,__LINE__);
        }

        return x_[x];
    }

    double Grid::getYposEdge(int y) const
    {
        if (y < 0 || y > m_) {
            ERROR("Y index out of bounds!",__FILE__,__LINE__);
        }

        return yv_[y];
    }

    double Grid::getYposCenter(int y) const
    {
        if (y < 0 || y >= m_) {
            ERROR("Y index out of bounds!",__FILE__,__LINE__);
        }

        return y_[y];
    }

    double Grid::getZposEdge(int z) const
    {
        if (z < 0 || z >= l_) {
            ERROR("Z index out of bounds!",__FILE__,__LINE__);
        }

        double dz = (zmax_ - zmin_) / l_;
        return -hdim_ * dz * cumulativeDepth_[z];
    }

    double Grid::getZposCenter(int z) const
    {
        if (z < 0 || z >= l_) {
            ERROR("Z index out of bounds!",__FILE__,__LINE__);
        }

        double dz = (zmax_ - zmin_) / l_;
        double result = cumulativeDepth_[z+1] + 0.5 * depth_[z];
        return -hdim_ * dz * result;
    }

    double Grid::fz(double z)
    {
        if (qz_ > 1.0)
            return -1.0 + tanh(qz_ * (z + 1.0)) / tanh(qz_);
        else
            return z;
    }

    double Grid::dfdz(double z)
    {
        double ch = cosh(qz_ * (z+1));
        if (qz_ > 1.0)
            return qz_ / (tanh(qz_) * ch * ch);
        else
            return 1.0 + (1.0 - qz_) * (1.0 - (2.0 * z));
    };

    /* Constructor
       input: dimensions of the global box:
       - N: east-west
       - M: north-south
       - L: z-direction
    */
    Domain::Domain(int N, int M, int L, int dof,
                   double Xmin, double Xmax, double Ymin, double Ymax,
                   bool Periodic, double Hdim, double qz, Teuchos::RCP<Epetra_Comm> Comm,
                   int aux)
        :
        comm(Comm),
        n(N), m(M), l(L),
        periodic(Periodic),
        dof_(dof),
        aux_(aux),
        gridGlb_(qz, N, M, L, Xmin, Xmax, Ymin, Ymax, Hdim)
    {
        int dim = m * n * l * dof_ + aux_;
        int *MyGlobalElements = new int[dim];

        for (int i = 0; i < dim; i++)
            MyGlobalElements[i] = i;

        ColMap = Teuchos::rcp(new Epetra_Map(dim, dim, MyGlobalElements, 0, *comm));
        delete [] MyGlobalElements;
    }

    // Destructor
    Domain::~Domain()
    {
        // destructor handled by Teuchos::rcp's
    }

    //=============================================================================
    // decompose the domain for a 2D processor array.
    // The xy-directions are split up.
    void Domain::Decomp2D()
    {
        int nprocs = comm->NumProc();
        int pid    = comm->MyPID();

        npL = 1;

        // Factor the number of processors into two dimensions. (nprocs = npN * npM)

        int t1 = nprocs;
        int t2 = 1;
        npM    = t1;
        npN    = t2;

        double r; //remainder
        double r_min = 100;

        while (t1 > 0)
        {
            t2 = (int)(nprocs / t1);
            r  =  std::abs(m/t1 - n/t2);
            if (t1 * t2 == nprocs && r <= r_min)
            {
                r_min = r;
                npM   = t1;
                npN   = t2;
            }
            t1--;
        }

        INFO("\n+++ 2D Domain decomposition +++");
        INFO(" factoring, np = " << nprocs);
        INFO("  n = "   << n);
        INFO("  m = "   << m);
        INFO("  npN = " << npN);
        INFO("  npM = " << npM << std::endl);

        // find out where in the domain we are situated. the
        // subdomains are numbered in a row-major 'matrix' fashion
        //
        // P0 P1 P2  |
        //           m,pidM
        //           |
        // P3 P4 P5  v
        //   --n->
        //   pidN
        //
        // note that this corresponds to the column-major ordering in
        // Fortran, i.e. i (n-direction) is the fastest index, and k
        // (l-dir.) the slowest

        pidL = 0; // for reasons of generality
        pidN = pid % npN;
        pidM = (pid - pidN) / npN;

        // dimension of actual subdomain (without ghost-nodes)
        mloc0 = (int) (m / npM);
        nloc0 = (int) (n / npN);
        lloc0 = l;

        // offsets for local->global index conversion
        Loff0 = 0;
        Moff0 = pidM * (int) (m / npM);
        Noff0 = pidN * (int) (n / npN);

        // distribute remaining points among first few cpu's
        int remM = m%npM;
        int remN = n%npN;

        if (pidM<remM) mloc0++;
        if (pidN<remN) nloc0++;
        for (int i=0;i<std::min(remM,pidM);i++) Moff0++;
        for (int i=0;i<std::min(remN,pidN);i++) Noff0++;

        //subdomain dimensions/offsets including ghost-nodes (will be
        // added further down)
        mloc=mloc0;
        nloc=nloc0;
        lloc=lloc0;

        Loff = Loff0;
        Moff = Moff0;
        Noff = Noff0;

        // we need two layers of overlap (ghost-nodes) between subdomains:
        //  _________............+
        // |      :  |           :
        // | SD11 :  |  SD12     :
        // |      :  |           :
        // +---------+-----------+
        // (for each variable, obviously)
        //
        // We need two because THCM assumes one layer of 'LAND' cells
        // at the domain boundaries. These are discarded only when
        // assembling global distributed vectors/matrices.

        // if there is only one subdomain in the
        // x-direction, periodicity is left to THCM
        Teuchos::RCP<Epetra_Comm> xcomm = this->GetProcRow(0);
        xparallel = (xcomm->NumProc()>1);

        if (pidM > 0)
        { mloc+=numGhosts; Moff-=numGhosts;}
        if (pidM < npM-1)
        { mloc+=numGhosts;}
        if ((pidN > 0) || (periodic && xparallel))
        { nloc+=numGhosts; Noff-=numGhosts;}
        if ((pidN < npN-1) || (periodic && xparallel))
        { nloc+=numGhosts;}

        // in the case of periodic boundary conditions the offsets may now be
        // negative or the local domain may exceed the global one. when
        // using nloc and noff, we therefore have to take mod(i,nglob)
        CommonSetup();
    }

    void Domain::CommonSetup()
    {
        INFO("processor position: (N,M,L) = ("<<pidN<<","<<pidM<<","<<pidL<<")");
        INFO("subdomain offsets: "<<Noff0<<","<<Moff0<<","<<Loff0);
        INFO("grid dimension on subdomain: "<<nloc0<<"x"<<mloc0<<"x"<<lloc0);
        INFO("  +++ including ghost nodes: +++");
        INFO("  subdomain offsets: "<<Noff<<","<<Moff<<","<<Loff);
        INFO("  grid dimension on subdomain: "<<nloc<<"x"<<mloc<<"x"<<lloc);
        INFO("  +++   auxiliary unknowns:  +++  " << aux_);

        // create the maps:
        StandardMap = CreateStandardMap(dof_);
        AssemblyMap = CreateAssemblyMap(dof_);

        // create surface single unknown maps:
        StandardSurfaceMap = CreateStandardMap(1, true);
        AssemblySurfaceMap = CreateAssemblyMap(1, true);

        // no load-balancing object available, yet (has to be set by user)
        SolveMap = StandardMap;

        // finally make the Import/Export objects (transfer function
        // between the two maps)
        as2std =
            Teuchos::rcp(new Epetra_Import(*AssemblyMap,
                                           *StandardMap));
        as2std_surf =
            Teuchos::rcp(new Epetra_Import(*AssemblySurfaceMap,
                                           *StandardSurfaceMap));

        std2sol = Teuchos::null;

        // Create local grid (including ghost nodes)
        gridLoc_ = gridGlb_.SubGrid(Noff, Moff, nloc, mloc, lloc);
    }

    // find out wether a particular local index is on a ghost node
    bool Domain::IsGhost(int ind, int nun_) const
    {
        int i,j,k,xx;
        int row = ind;

        bool result = false;

        // find out where the point is
        Utils::ind2sub(nloc,mloc,lloc,nun_,row,i,j,k,xx);

        // ghost nodes at periodic boundary only if more than one proc in x-direction
        bool perio = periodic && xparallel;

        for (int ii = 0; ii < numGhosts; ii++)
        {
            result = result||(i==ii && ((pidN>0)||perio));
            result = result||(i==nloc-1-ii && ((pidN<npN-1)||perio));
            result = result||(j==ii && pidM>0);
            result = result||(j==mloc-1-ii && pidM<npM-1);
            result = result||(k==ii && pidL>0);
            result = result||(k==lloc-1-ii && pidL<npL-1);
        }

        return result;
    }

    // public map creation function (can only create a limited range of maps)
    Teuchos::RCP<Epetra_Map> Domain::CreateSolveMap(int nun_, bool depth_av) const
    {
        Teuchos::RCP<Epetra_Map> M=Teuchos::null;
        if (UseLoadBalancing()==false) // no load-balancing? use our own decomposition
        {
            M=CreateStandardMap(nun_,depth_av);
        }
        else
        {
            ERROR("not implemented",__FILE__,__LINE__);
        }
        return M;
    }

    //=============================================================================
    Teuchos::RCP<Epetra_Map> Domain::CreateStandardMap(int nun_, bool depth_av) const
    {
        // Add auxiliary unknowns at final processor,
        int root    = comm->NumProc() - 1;
        bool addAux = (comm->MyPID() == root) ? true : false;

        Teuchos::RCP<Epetra_Map> M = Teuchos::null;
        if (depth_av)
        {
            M = CreateMap(Noff0, Moff0,0,nloc0, mloc0, 1, nun_);
        }
        else
        {
            M = CreateMap(Noff0, Moff0,Loff0,nloc0, mloc0, lloc0, nun_, addAux);
        }
        return M;
    }

    //==========================================================================
    // Public map creation function (can only create a limited range of maps)
    Teuchos::RCP<Epetra_Map> Domain::CreateAssemblyMap(int nun_, bool depth_av) const
    {
        // Add auxiliary unknowns at every processor
        bool addAux = true;

        Teuchos::RCP<Epetra_Map> M = Teuchos::null;
        if (depth_av) // create a map with a single vertical layer: 'depth-averaged'
        {
            M = CreateMap(Noff, Moff, 0, nloc, mloc, 1, nun_);
        }
        else
        {
            M = CreateMap(Noff, Moff, Loff, nloc, mloc, lloc, nun_, addAux);
        }
        return M;
    }

    //==========================================================================
    // This version is private and very general
    Teuchos::RCP<Epetra_Map> Domain::CreateMap(int noff_, int moff_, int loff_,
                                               int nloc_, int mloc_, int lloc_,
                                               int nun_,  bool addAux) const
    {
        // If a map is requested for all unknowns, we append the auxiliary unknowns as well.
        // Otherwise, we create the original map.
        int NumMyElements;
        if (addAux)
            NumMyElements = mloc_ * nloc_ * lloc_ * nun_ + aux_;
        else
            NumMyElements = mloc_ * nloc_ * lloc_ * nun_;


        int NumGlobalElements = -1; // Let Epetra figure it out herself

        // Make lists of all global elements that belong to my subdomain
        // (with and without ('0') ghost nodes)
        int* MyGlobalElements = new int[NumMyElements];
        int p;

        // construct list of all entries
        p = 0;
        for (int k = loff_; k < loff_ + lloc_; k++)
        {
            for (int j = moff_; j < moff_ + mloc_; j++)
            {
                for (int i = noff_; i < noff_ + nloc_; i++)
                {
                    // this is to handle periodic bc in the x-direction:
                    int ii = MOD((double) i, (double) n);
                    for (int xx = 1; xx <= nun_; xx++) //nun_ unknowns per element
                    {
                        MyGlobalElements[p] = FIND_ROW2(nun_, n, m, l, ii, j, k, xx);
                        p++;
                    }//xx
                }//i
            }//j
        }//k

        // Add auxiliary unknowns, every subdomain needs them, they're great
        if (addAux)
        {
            for (int aa = 1; aa <= aux_; ++aa)
            {
                MyGlobalElements[p] = FIND_ROW2(nun_, n, m, l, n-1, m-1, l-1, nun_) + aa;
                p++;
            }
        }

        // A little integrity check
        assert(p == NumMyElements);

        // create the map
        Teuchos::RCP<Epetra_Map> M =
            Teuchos::rcp(new Epetra_Map(NumGlobalElements,
                                        NumMyElements, MyGlobalElements, 0, *comm));
        delete [] MyGlobalElements;
        return M;
    }

    //==========================================================================
    // create communication groups
    Teuchos::RCP<Epetra_Comm> Domain::GetProcRow(int dim)
    {
#ifndef HAVE_MPI
        return comm; // can't split anything in sequential mode
#else //{
        Teuchos::RCP<Epetra_MpiComm> mpi_comm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(comm);
        if (mpi_comm==Teuchos::null) ERROR("Bad Communicator encountered!",__FILE__,__LINE__);
        MPI_Comm old_comm = mpi_comm->GetMpiComm();
        MPI_Comm row_comm;
        MPI_Group old_group, row_group;

        int nproc_row, disp, offset;

        if (npL>1) ERROR("This function is not implemented for 3D proc arrays!",
                         __FILE__,__LINE__);


        if (dim==0)
        {
            nproc_row=npN;
            disp=1;
            offset = pidM*npN;
        }
        else if (dim==1)
        {
            nproc_row=npM;
            disp=npN;
            offset = pidN;
        }
        else
        {
            INFO("Bad dimension to be extracted from proc array: "<<dim);
            ERROR("Cannot split dimension",__FILE__,__LINE__);
        }
        int *row_ranks = new int[nproc_row];


//  Create vector of ranks for row processes:
        DEBUG("My position: ("<<pidN<<", "<<pidM<<")");
        DEBUG("extract communicator for dim="<<dim);
        DEBUG("Creating sub-communicator consisting of:");
        for (int k = 0; k < nproc_row; k++)
        {
            row_ranks[k] = offset+k*disp;
            DEBUG(k<<": "<<row_ranks[k]);
        }
        /* ----------------------- */
        /* create the global group */
        /* ----------------------- */
        int ierr=MPI_Comm_group(old_comm, &old_group);

        /* --------------------- */
        /* Create the row group  */
        /* --------------------- */
        ierr = MPI_Group_incl(old_group, nproc_row, row_ranks, &row_group);
        if (ierr!=0) ERROR("MPI call 'Group_incl' failed!",__FILE__,__LINE__);

        /* ---------------------------------------------------- */
        /* Create the new communicator; MPI assigns the context */
        /* ---------------------------------------------------- */

        ierr = MPI_Comm_create(old_comm, row_group, &row_comm);
        if (ierr!=0) ERROR("MPI call 'Comm_create' failed!",__FILE__,__LINE__);

        // finally create the Epetra comm object
        Teuchos::RCP<Epetra_Comm> new_comm = Teuchos::rcp(new Epetra_MpiComm(row_comm));

        delete [] row_ranks;
        return new_comm;
#endif //}
    }

/////////////////////////////////////////////////////////////////////////////////
// Data Migration functions                                                    //
/////////////////////////////////////////////////////////////////////////////////

    //
    int Domain::Assembly2Standard
    (const Epetra_Vector& source, Epetra_Vector& target) const
    {
#ifdef DEBUGGING_NEW
        if (!(source.Map().SameAs(*AssemblyMap) &&
              target.Map().SameAs(*StandardMap)))
        {
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        CHECK_ZERO(target.Export(source,*as2std,Zero));
        return 0;
    }

    //
    int Domain::Standard2Assembly
    (const Epetra_Vector& source, Epetra_Vector& target) const
    {
#ifdef DEBUGGING_NEW
        if (!(source.Map().SameAs(*StandardMap) &&
              target.Map().SameAs(*AssemblyMap)))
        {
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        CHECK_ZERO(target.Import(source,*as2std,Insert));
        return 0;
    }

    int Domain::Assembly2StandardSurface
    (const Epetra_Vector& source, Epetra_Vector& target) const
    {
#ifdef DEBUGGING_NEW
        if (!(source.Map().SameAs(*AssemblySurfaceMap) &&
              target.Map().SameAs(*StandardSurfaceMap)))
        {
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        CHECK_ZERO(target.Export(source,*as2std_surf,Zero));
        return 0;
    }

    //
    int Domain::Standard2AssemblySurface
    (const Epetra_Vector& source, Epetra_Vector& target) const
    {
#ifdef DEBUGGING_NEW
        if (!(source.Map().SameAs(*StandardSurfaceMap) &&
              target.Map().SameAs(*AssemblySurfaceMap)))
        {
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        CHECK_ZERO(target.Import(source,*as2std_surf,Insert));
        return 0;
    }

    //
    int Domain::Standard2Solve
    (const Epetra_Vector& source, Epetra_Vector& target) const
    {
#ifdef DEBUGGING_NEW
        if (!(source.Map().SameAs(*StandardMap)&&target.Map().SameAs(*SolveMap)))
        {
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        if (UseLoadBalancing()==false)
        {
            target = source;
        }
        else
        {
            CHECK_ZERO(target.Export(source,*std2sol,Insert));
        }
        return 0;
    }

    //
    int Domain::Solve2Standard
    (const Epetra_Vector& source, Epetra_Vector& target) const
    {
#ifdef DEBUGGING_NEW
        if (!(source.Map().SameAs(*SolveMap)&&target.Map().SameAs(*StandardMap)))
        {
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        if (UseLoadBalancing()==false)
        {
            target = source;
        }
        else
        {
            CHECK_ZERO(target.Import(source,*std2sol,Insert));
        }
        return 0;
    }

    //
    int Domain::Solve2Assembly
    (const Epetra_Vector& source, Epetra_Vector& target) const
    {
#ifdef DEBUGGING_NEW
        if (!(source.Map().SameAs(*SolveMap) && target.Map().SameAs(*AssemblyMap)))
        {
            INFO("source.Map() != SolveMap " << !source.Map().SameAs(*SolveMap));
            INFO("target.Map() != AssemblyMap " << !target.Map().SameAs(*AssemblyMap));
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        if (UseLoadBalancing()==false)
        {
            this->Standard2Assembly(source,target);
        }
        else
        {
            Epetra_Vector tmp(*StandardMap);
            this->Solve2Standard(source,tmp);
            this->Standard2Assembly(tmp,target);
        }
        return 0;
    }

    //
    int Domain::Assembly2Solve
    (const Epetra_Vector& source, Epetra_Vector& target) const
    {
//        DEBUG("Enter Assembly2Solve");
#ifdef DEBUGGING_NEW
        if (!(source.Map().SameAs(*AssemblyMap)&&target.Map().SameAs(*SolveMap)))
        {
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        if (UseLoadBalancing()==false)
        {
            this->Assembly2Standard(source,target);
        }
        else
        {
            Epetra_Vector tmp(*StandardMap);
            this->Assembly2Standard(source,tmp);
            this->Standard2Solve(tmp,target);
        }
        return 0;
    }



    //
    int Domain::Standard2Solve
    (const Epetra_CrsMatrix& source, Epetra_CrsMatrix& target) const
    {
#ifdef DEBUGGING_NEW
        if (!(source.RowMap().SameAs(*StandardMap)&&target.RowMap().SameAs(*SolveMap)))
        {
            ERROR("Invalid Transfer Function called!",__FILE__,__LINE__);
        }
#endif
        if (UseLoadBalancing()==false)
        {
            target = source;
        }
        else
        {
            CHECK_ZERO(target.Export(source,*std2sol,Insert));
        }
        return 0;
    }
}//namespace
