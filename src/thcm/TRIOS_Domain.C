/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

#include "TRIOS_Domain.H"
#include <fstream>

#include "TRIOS_Macros.H"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#define CLASSNAME "Domain"
#include "TRIOS_ftrace.H"

namespace TRIOS {

      /* Constructor
          input: dimensions of the global box:
            - N: east-west
            - M: north-south
            - L: z-direction
      */
Domain::Domain(int N, int M, int L, int dof, 
                       double Xmin, double Xmax, double Ymin, double Ymax, 
                       bool Periodic, double Hdim, Teuchos::RCP<Epetra_Comm> Comm)
  : m(M),n(N),l(L),xmin(Xmin),xmax(Xmax),ymin(Ymin),ymax(Ymax),
    dof_(dof),
    periodic(Periodic),zmin(-Hdim),zmax(0),
    comm(Comm)
  {
  ENTER("Domain")
  //TODO: check if we want zmin=-Hdim or -1 (as in THCM)
  int dim = m*n*l*dof_;
  int *MyGlobalElements = new int[dim];
  for (int i=0;i<dim;i++) MyGlobalElements[i]=i;
  ColMap = Teuchos::rcp(new Epetra_Map(dim,dim,MyGlobalElements,0,*comm) );
  delete [] MyGlobalElements;
  LEAVE("Domain");
  }
      
// Destructor
Domain::~Domain()
  {
  ENTER("~Domain");
  // destructor handled by Teuchos::rcp's
  LEAVE("~Domain");
  }
      
                  
// decompose the domain for a 2D processor array.
// The xy-directions are split up.
void Domain::Decomp2D()
  {
  ENTER("Decomp2D");
  int nprocs = comm->NumProc();
  int pid = comm->MyPID();

////////////////////////////////////////////////////////////////////

  npL = 1;

// Factor the number of processors into two dimensions. (nprocs = npN*npM)

  int t1 = nprocs;
  int t2 = 1;
  npM=t1;
  npN=t2;

  double r;//remainder
  double r_min = 100;

  while (t1>0)
    {
    t2=(int)(nprocs/t1);
    r = std::abs(m/t1-n/t2);
    if (t1*t2==nprocs && r<=r_min)
      {
      r_min=r;
      npM=t1;
      npN=t2;
      }
//    (*info) << t1 << " " << t2 << " " << r << std::endl;
    t1--;
    }

  INFO("+++ Domain decomposition +++");
  INFO("factoring, np = "<<nprocs);
  INFO(" n = "<< n);
  INFO(" m = "<< m);
  INFO(" npN = "<< npN);
  INFO(" npM = "<< npM << std::endl);

                        
  // find out where in the domain we are situated.
  // the subdomains are numbered in a row-major 'matrix' fashion
  //                  
  // P0 P1 P2  |      
  //           m,pidM 
  //           |      
  // P3 P4 P5  v      
  //  --n-->          
  //   pidN           
  //                  
  
  // note that this corresponds to the column-major ordering in Fortran,
  // i.e. i (n-direction) is the fastest index, and k (l-dir.) the slowest
        
  pidL = 0; // for reasons of generality
  pidN = pid%npN;
  pidM = (pid-pidN)/npN;

  // dimension of actual subdomain (without ghost-nodes)
        
  mloc0 = (int)(m/npM);
  nloc0 = (int)(n/npN);
  lloc0 = l;
  
  // offsets for local->global index conversion
  Loff0 = 0;
  Moff0 = pidM*(int)(m/npM);
  Noff0 = pidN*(int)(n/npN);

  // distribute remaining points among first few cpu's
  int remM=m%npM;
  int remN=n%npN;
  
  if (pidM<remM) mloc0++;
  if (pidN<remN) nloc0++;
  for (int i=0;i<std::min(remM,pidM);i++) Moff0++;
  for (int i=0;i<std::min(remN,pidN);i++) Noff0++;

  //subdomain dimensions/offsets including ghost-nodes
  // (will be added further down)
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
  
  Teuchos::RCP<Epetra_Comm> xcomm = this->GetProcRow(0);
  xparallel = (xcomm->NumProc()>1); // if there is only one subdomain in the 
                                        // x-direction, periodicity is left to THCM

  if (pidM>0)      { mloc+=num_ghosts; Moff-=num_ghosts;}
  if (pidM<npM-1) { mloc+=num_ghosts;}
  if ((pidN>0)||(periodic&&xparallel))      { nloc+=num_ghosts; Noff-=num_ghosts;}
  if ((pidN<npN-1)||(periodic&&xparallel)) { nloc+=num_ghosts;}

// in the case of periodic boundary conditions the offsets may now be 
// negative or the local domain may exceed the global one. when       
// using nloc and noff, we therefore have to take mod(i,nglob)

  CommonSetup();		
  LEAVE("Decomp2D");
  }

// 3D domain decomposition (for Navier-Stokes)
// This implementation is more modern than the Decomp2D function
// because we have lots of functionality from HYMLS we can use now.
// NOTE: we do not include overlap if this function is used, that 
// would have to be added if we have a parallel code using 3D decom-
// positions.
void Domain::Decomp3D()
  {
  ENTER("Decomp3D");
  int nprocs = comm->NumProc();
  int pid = comm->MyPID();
  HYMLS::Tools::SplitBox(n,m,l,nprocs,npN,npM,npL);
  HYMLS::Tools::ind2sub(npN,npM,npL,pid,pidN,pidM,pidL);

  Loff0 = pidL*(int)(l/npL);
  Moff0 = pidM*(int)(m/npM);
  Noff0 = pidN*(int)(n/npN);
       
 Noff = Noff0; 
 Moff = Moff0; 
 Loff = Loff0;

  mloc0 = (int)(m/npM);
  nloc0 = (int)(n/npN);
  lloc0 = (int)(l/npL);
 
 nloc = nloc0;
 mloc = mloc0;
 lloc = lloc0;

  CommonSetup();  
  LEAVE("Decomp3D");
  }

  void Domain::CommonSetup()
    {
    ENTER("CommonSetup");

INFO("processor position: (N,M,L) = ("<<pidN<<","<<pidM<<","<<pidL<<")");
INFO("subdomain offsets: "<<Noff0<<","<<Moff0<<","<<Loff0);
INFO("grid dimension on subdomain: "<<nloc0<<"x"<<mloc0<<"x"<<lloc0);
DEBUG("+++ including ghost nodes: +++");
DEBUG("subdomain offsets: "<<Noff<<","<<Moff<<","<<Loff);
DEBUG("grid dimension on subdomain: "<<nloc<<"x"<<mloc<<"x"<<lloc);


  // create the maps:
       
  StandardMap = CreateStandardMap(dof_);
  AssemblyMap = CreateAssemblyMap(dof_);
  // no load-balancing object available, yet (has to be set by user)
  SolveMap = StandardMap;


#ifdef DEBUGGING
comm->Barrier();
DEBUG("create importers...");
#endif

  // finally make the Import/Export objects (transfer function
  // between the two maps)
  as2std = Teuchos::rcp(new Epetra_Import(*AssemblyMap,*StandardMap));
  std2sol = Teuchos::null;
/*
  INFO("importer: "<<*as2std);
  Epetra_Vector test_as(*AssemblyMap);
  Epetra_Vector test_so(*SolveMap);
  
  for (int i=0;i<test_as.MyLength();i++)
    {
    test_as[i]=i;
    }
  this->Assembly2Solve(test_as,test_so);
  INFO("assembly vector: "<<test_as);
  INFO("solve vector: "<<test_so);
  this->Solve2Assembly(test_so,test_as);
  INFO("new assembly vector: "<<test_as);
  

//  DEBUG("Importer: "<<std::endl)
//  DEBUG( (*Importer) );
*/

  // determine the physical bounds of the subdomain
  // (must be passed to THCM)
  
  // grid constants as computed in 'grid.f':
  double dx = (xmax-xmin)/n;
  double dy = (ymax-ymin)/m;
  double dz = (zmax-zmin)/l;
        
  xmin_loc = xmin + Noff*dx;
  xmax_loc = xmin + (Noff+nloc)*dx;
  ymin_loc = ymin + Moff*dy;
  ymax_loc = ymin + (Moff+mloc)*dy;
  zmin_loc = zmin + Loff*dz;
  zmax_loc = zmin + (Loff+lloc)*dz;

DEBVAR(xmin)
DEBVAR(xmax)
DEBVAR(ymin)
DEBVAR(ymax)
DEBVAR(xmin_loc)
DEBVAR(xmax_loc)
DEBVAR(ymin_loc)
DEBVAR(ymax_loc)
  LEAVE("CommonSetup");
    }

  // find out wether a particular local index is on a ghost node
  bool Domain::IsGhost(int ind,int nun_) const
    {
    int i,j,k,xx;
    int row = ind; 

    bool result = false;


//    DEBUG("GHOST CHECK: li (C++) = "<<ind);

    // find out where the point is
    HYMLS::Tools::ind2sub(nloc,mloc,lloc,nun_,row,i,j,k,xx);
    
    // ghost nodes at periodic boundary only if more than one proc in x-direction
    bool perio = periodic&&xparallel;

//    DEBUG("GHOST CHECK: i,j,k,xx (C++) = "<<i<<", "<<j<<", "<<k<<", "<<xx);

    for (int ii=0;ii<num_ghosts;ii++)
      {
      result = result||(i==ii && ((pidN>0)||perio));
      result = result||(i==nloc-1-ii && ((pidN<npN-1)||perio));
      result = result||(j==ii && pidM>0);
      result = result||(j==mloc-1-ii && pidM<npM-1);
      result = result||(k==ii && pidL>0);
      result = result||(k==lloc-1-ii && pidL<npL-1);
      }

//    DEBVAR(result);
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
    HYMLS::Tools::Error("not implemented",__FILE__,__LINE__);
//    int l_=depth_av?1:l;
//    M = loadbal->createBalancedMap(l_,nun_);
    }
  return M;
  }

Teuchos::RCP<Epetra_Map> Domain::CreateStandardMap(int nun_, bool depth_av) const
  {
  Teuchos::RCP<Epetra_Map> M = Teuchos::null;
  if (depth_av)
    {
    M = CreateMap(Noff0, Moff0,0,nloc0, mloc0, 1, nun_);
    }
  else
    {
    M = CreateMap(Noff0, Moff0,Loff0,nloc0, mloc0, lloc0, nun_);
    }
  return M;
  }

// public map creation function (can only create a limited range of maps)
Teuchos::RCP<Epetra_Map> Domain::CreateAssemblyMap(int nun_, bool depth_av) const
  {
  Teuchos::RCP<Epetra_Map> M=Teuchos::null;
  if (depth_av)
    {
    M = CreateMap(Noff, Moff,0,nloc, mloc, 1, nun_);
    }
  else
    {
    M = CreateMap(Noff, Moff,Loff,nloc, mloc, lloc, nun_);
    }
  return M;
  }
  
// this version is private and very general
Teuchos::RCP<Epetra_Map> Domain::CreateMap(int noff_, int moff_, int loff_,
                              int nloc_, int mloc_, int lloc_, int nun_) const
  {
  int NumMyElements = mloc_*nloc_*lloc_*nun_;
  DEBVAR(NumMyElements)  
  int NumGlobalElements = -1; // Let Epetra figure it out herself
                
        
  // Make lists of all global elements that belong to my subdomain 
  // (with and without ('0') ghost nodes)
  int* MyGlobalElements = new int[NumMyElements];
  int p;
        
  // construct list of all entries
  p = 0;
  for (int k=loff_;k<loff_+lloc_;k++)
    {
    for (int j=moff_;j<moff_+mloc_;j++)
      {
      for (int i=noff_;i<noff_+nloc_;i++)
        {
        // this is to handle periodic bc in the x-direction:
        int ii = MOD((double)i,(double)n);
        for (int xx=1;xx<=nun_;xx++) //6 unknowns per element
          {
//          DEBUG(ii<<" "<<j<<" "<<k<<" => "<<FIND_ROW2(nun_,n,m,l,i,j,k,xx));
          MyGlobalElements[p] = FIND_ROW2(nun_,n,m,l,ii,j,k,xx);
          p++;
          }//xx
        }//i
      }//j
    }//k



  // create the map
       
   Teuchos::RCP<Epetra_Map> M = Teuchos::rcp(new Epetra_Map(NumGlobalElements,
          NumMyElements,MyGlobalElements,0,*comm));
          
  DEBVAR(NumMyElements)  
   delete [] MyGlobalElements;
   return M;


}


// create communication groups
Teuchos::RCP<Epetra_Comm> Domain::GetProcRow(int dim)
  {
#ifndef HAVE_MPI
  return comm; // can't split anything in sequential mode
#else //{
  Teuchos::RCP<Epetra_MpiComm> mpi_comm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(comm);
  if (mpi_comm==Teuchos::null) HYMLS::Tools::Error("Bad Communicator encountered!",__FILE__,__LINE__);
  MPI_Comm old_comm = mpi_comm->GetMpiComm();
  MPI_Comm row_comm;
  MPI_Group old_group, row_group;

    int nproc_row, disp, offset;

    if (npL>1) HYMLS::Tools::Error("This function is not implemented for 3D proc arrays!",
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
      HYMLS::Tools::Error("Cannot split dimension",__FILE__,__LINE__);
      }
    int *row_ranks = new int[nproc_row];


//  Create vector of ranks for row processes: 
    DEBUG("My position: ("<<pidN<<", "<<pidM<<")");
    DEBUG("extract communicator for dim="<<dim);
    DEBUG("Creating sub-communicator consisting of:")
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
  ierr=MPI_Group_incl(old_group, nproc_row, row_ranks, &row_group);
  if (ierr!=0) HYMLS::Tools::Error("MPI call 'Group_incl' failed!",__FILE__,__LINE__);

  /* ---------------------------------------------------- */
  /* Create the new communicator; MPI assigns the context */
  /* ---------------------------------------------------- */
  
  ierr = MPI_Comm_create(old_comm, row_group, &row_comm);
  if (ierr!=0) HYMLS::Tools::Error("MPI call 'Comm_create' failed!",__FILE__,__LINE__);
  
  // finally create the Epetra comm object
  Teuchos::RCP<Epetra_Comm> new_comm = Teuchos::rcp(new Epetra_MpiComm(row_comm),false);
  
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
//        DEBUG("Enter Assembly2Standard");
#ifdef TESTING
        if (!(source.Map().SameAs(*AssemblyMap)&&target.Map().SameAs(*StandardMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
          }
#endif
        CHECK_ZERO(target.Export(source,*as2std,Zero));
//        DEBUG("Leave Assembly2Standard");
        return 0;
        }


      //
      int Domain::Standard2Assembly
        (const Epetra_Vector& source, Epetra_Vector& target) const
        {
//        DEBUG("Enter Standard2Assembly");
#ifdef TESTING
        if (!(source.Map().SameAs(*StandardMap)&&target.Map().SameAs(*AssemblyMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
          }
#endif
        CHECK_ZERO(target.Import(source,*as2std,Insert));
//        DEBUG("Leave Standard2Assembly");
        return 0;
        }

      //
      int Domain::Standard2Solve
        (const Epetra_Vector& source, Epetra_Vector& target) const
        {
//        DEBUG("Enter Standard2Solve");
#ifdef TESTING
        if (!(source.Map().SameAs(*StandardMap)&&target.Map().SameAs(*SolveMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
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
//        DEBUG("Leave Standard2Solve");
        return 0;
        }

      //
      int Domain::Solve2Standard
        (const Epetra_Vector& source, Epetra_Vector& target) const
        {
//        DEBUG("Enter Solve2Standard");
#ifdef TESTING
        if (!(source.Map().SameAs(*SolveMap)&&target.Map().SameAs(*StandardMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
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
//        DEBUG("Leave Solve2Standard");
        return 0;
        }

      //
      int Domain::Solve2Assembly
         (const Epetra_Vector& source, Epetra_Vector& target) const
        {
//        DEBUG("Enter Solve2Assembly");
#ifdef TESTING
        if (!(source.Map().SameAs(*SolveMap)&&target.Map().SameAs(*AssemblyMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
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
//        DEBUG("Leave Solve2Assembly");
        return 0;
        }


      //
      int Domain::Assembly2Solve
        (const Epetra_Vector& source, Epetra_Vector& target) const
        {
//        DEBUG("Enter Assembly2Solve");
#ifdef TESTING
        if (!(source.Map().SameAs(*AssemblyMap)&&target.Map().SameAs(*SolveMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
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
//        DEBUG("Leave Assembly2Solve");
        return 0;
        }
        
        

      //
      int Domain::Standard2Solve
        (const Epetra_CrsMatrix& source, Epetra_CrsMatrix& target) const
        {
//        DEBUG("Enter Standard2Solve");
#ifdef TESTING
        if (!(source.RowMap().SameAs(*StandardMap)&&target.RowMap().SameAs(*SolveMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
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
//        DEBUG("Enter Leave2Solve");
        return 0;
        }
}//namespace
