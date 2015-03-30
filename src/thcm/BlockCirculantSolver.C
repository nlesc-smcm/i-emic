#include "BlockCirculantSolver.H"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_RowMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "math.h"
#include "HYMLS_MatrixUtils.H"
#include "Amesos.h"
#include "EpetraExt_BlockVector.h"

#include "fstream" //de header om bestanden te kunnen maken.
#include "iostream"//de header voor cin.get();

//#ifdef HAVE_FFTW

#include "Komplex_LinearProblem.h"
//#endif

//#define FILE_DEBUG 1

#ifndef INFO
#define INFO(s) std::cout << s << std::endl;
#endif

#ifndef DEBUG
#define DEBUG(s)
#endif

using Teuchos::RCP;
using Teuchos::rcp;


//! Constructor
BlockCirculantSolver::
BlockCirculantSolver(const EpetraExt::BlockCrsMatrix& cloneMatrix_,
         const EpetraExt::BlockVector& cloneVector_,
         const Teuchos::RCP<const Epetra_CrsMatrix> jacobianBlock_, 
	 const Teuchos::RCP<const Epetra_CrsMatrix> massBlock_,
	 const Teuchos::RCP<const EpetraExt::MultiComm> globalComm_,
	 double dt_, Teuchos::ParameterList& params_) :
  jacobianBlock(jacobianBlock_),
  massBlock(massBlock_),
  globalComm(globalComm_),
  dt(dt_), params(params_),
  jacobian(cloneMatrix_),
  nullSpace_(Teuchos::null)
{
  DEBUG("Enter BlockCirculantSolver constructor");

  this->buildTimeComm();

  nullSpace_ = params.get("Null Space", nullSpace_);

  localComm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm> (rcp(&(globalComm->SubDomainComm()),false));

  numBlocks=globalComm->NumTimeSteps();
  blocksOnDomain=globalComm->NumTimeStepsOnDomain();
  blockSize=jacobianBlock->NumMyRows();

  globalMap=Teuchos::rcp(&(cloneMatrix_.RowMap()),false);
  standardMap=Teuchos::rcp(&(jacobianBlock->RowMap()),false);


  nullSpace_ = params.get("Null Space", nullSpace_);
  
  // allocate arrays
  rhsRe.resize(blocksOnDomain);
  lhsRe.resize(blocksOnDomain);
  rhsIm.resize(blocksOnDomain);
  lhsIm.resize(blocksOnDomain);

  Solver.resize(blocksOnDomain);

  complexProblem.resize(blocksOnDomain);

  // create all the objects: 
  
  Teuchos::ParameterList& solverList
        = params.sublist("Ifpack");

  std::string solver_type = solverList.get("amesos: solver type","Amesos_Klu");
  Amesos Factory;

  for (int i=0;i<blocksOnDomain;i++)
    {
    rhsRe[i]=Teuchos::rcp(new Epetra_Vector(*standardMap));
    rhsIm[i]=Teuchos::rcp(new Epetra_Vector(*standardMap));
    lhsRe[i]=Teuchos::rcp(new Epetra_Vector(*standardMap));
    lhsIm[i]=Teuchos::rcp(new Epetra_Vector(*standardMap));
    complexProblem[i] = rcp(new Komplex_LinearProblem(1.0,0.0,*jacobianBlock,0.0,1.0,*massBlock, *(lhsRe[i]),*(lhsIm[i]),*(rhsRe[i]),*(rhsIm[i])));
    Solver[i] = Teuchos::rcp(Factory.Create(solver_type.c_str(),*(complexProblem[i]->KomplexProblem())));
    Solver[i]->SetParameters(solverList);
    }
// at this stage the preconditioner is not computed, yet.
// A doesn't contain the circulant matrix and the LU decomposition
// is not yet computed. All of that is done in ComputePreconditioner()

  // the FFT setup can be done in the constructor, it doesn't need any info on the
  // matrix except the dimensions.
  int NX=numBlocks;// we want to transform in this direction (sequentially)
  int NY=blockSize;
#ifdef HAVE_FFTW  
  planf = fftw_create_plan(NX, FFTW_FORWARD, FFTW_ESTIMATE);
  planb = fftw_create_plan(NX, FFTW_BACKWARD, FFTW_ESTIMATE);
  data = new fftw_complex[NX*NY];
  work = new fftw_complex[NX*NY];
#else
      std::string errorMessage = "FFTW is not available, use another preconditioner!";
      throwError("LOCA::Epetra::xyztPrec::xyztPrec", errorMessage);
#endif  
  DEBUG("Leave BlockCirculantSolver constructor");
}

BlockCirculantSolver::
~BlockCirculantSolver()
{
#ifdef HAVE_FFTW
  fftw_destroy_plan(planf);
	fftw_destroy_plan(planb);
  delete [] data;
  delete [] work;
#endif  
}


int BlockCirculantSolver::
SetUseTranspose(bool UseTranspose)
{
  // Disable this option
  return false;
}

int BlockCirculantSolver::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Not implemented: Throw an error!
  std::cout << "ERROR: BlockCirculantSolver::Apply() - "
			<< "method is NOT implemented!!  " << std::endl; 
  throw "LOCA Error";
  return false;
}

int BlockCirculantSolver::
ApplyInverse(const Epetra_MultiVector& input,
	     Epetra_MultiVector& result) const
  {

#ifdef HAVE_FFTW

  int NX=numBlocks;// we want to transform in this direction (sequentially)
  int NY=blockSize;
  double argexp;
  double eigReal;
  double eigImag;

  if (input.NumVectors()>1)
    {
    return -1; // only works for one input/output vector right now
    }

  double *realPointer;

  double *realArray = new double[blocksOnDomain*blockSize];
  double *imagArray = new double[blocksOnDomain*blockSize];

  double *globalRealArray = new double[NX*NY];
  double *globalImagArray = new double[NX*NY];

  EPETRA_CHK_ERR(input(0)->ExtractView(&realPointer));
  EPETRA_CHK_ERR(timeComm->GatherAll(realPointer,globalRealArray,blockSize*blocksOnDomain));

  // fill data array for FFT
  for (int i=0;i<NX;i++)
    {
    for (int j=0;j<NY;j++)
      {
      data[i*NY+j].re = globalRealArray[i*NY+j];
      data[i*NY+j].im = 0.0;
      }
    }

  //fft rhs
  fftw(planf, NY, data, NY, 1, work, NY, 1);

  // fill epetra vectors for later use
  for (int ii=0;ii<blocksOnDomain;ii++)
    {
    int i = globalComm->FirstTimeStepOnDomain()+ii;
    for (int j=0;j<NY;j++)
      {
      (*rhsRe[ii])[j]=work[i*blockSize+j].re;

      (*rhsIm[ii])[j]=work[i*blockSize+j].im;
			}
    }

  //Solve system per block
  for (int ii=0;ii<blocksOnDomain;ii++)
    {
    int i = globalComm->FirstTimeStepOnDomain()+ii;
    argexp = i*2.0*3.14159265358979323846/numBlocks;
    eigReal = cos(argexp);//dt al in massBlock
    eigImag = -sin(argexp);
    complexProblem[ii]->UpdateValues(1.0,0.0,*jacobianBlock,eigReal,eigImag,*massBlock,        *(lhsRe[ii]),*(lhsIm[ii]),*(rhsRe[ii]),*(rhsIm[ii]));
    EPETRA_CHK_ERR(Solver[ii]->Solve());
    complexProblem[ii]->ExtractSolution(*(lhsRe[ii]),*(lhsIm[ii]));

	 //subtract nullspace
   //this->PostProcSolution(*(lhsRe[ii]));
   //this->PostProcSolution(*(lhsIm[ii]));

    // put solution in local array for gathering
    lhsRe[ii]->ExtractCopy(realArray+ii*blockSize);
    lhsIm[ii]->ExtractCopy(imagArray+ii*blockSize);
    }

  EPETRA_CHK_ERR(timeComm->GatherAll(realArray,globalRealArray,blockSize*blocksOnDomain));
  EPETRA_CHK_ERR(timeComm->GatherAll(imagArray,globalImagArray,blockSize*blocksOnDomain));

  //PREPARE IFFT 
  // fill data array for FFT
  for (int i=0;i<NX;i++)
    {
    for (int j=0;j<NY;j++)
      {
      data[i*NY+j].re = globalRealArray[i*NY+j];
	    data[i*NY+j].im = globalImagArray[i*NY+j];
      }
		}

  //fft solution
  // do back-transform
  fftw(planb, NY, data, NY, 1, work, NY, 1);

  // fill epetra vectors for later use
  for (int ii=0;ii<blocksOnDomain;ii++)
    {
    int i = globalComm->FirstTimeStepOnDomain()+ii;
    for (int j=0;j<blockSize;j++)
      {
      result[0][ii*blockSize+j] = (work[i*blockSize+j].re)/numBlocks;
      }
    }
	
	delete [] realArray;
  delete [] imagArray;
  delete [] globalRealArray;
  delete [] globalImagArray;
#endif
  return 0;
  }


double BlockCirculantSolver::
NormInf() const
{
  // Not implemented: Throw an error!
  std::cout << "ERROR: BlockCirculantSolver::NormInf() - "
            << "method is NOT implemented!!  " << std::endl; 
  throw "LOCA Error";
  return 0.0;
}

bool BlockCirculantSolver::
UseTranspose() const
{
  // Not implemented: Throw an error!
  std::cout << "ERROR: BlockCirculantSolver::UseTranspose() - "
            << "method is NOT implemented!!  " << std::endl; 
  throw "LOCA Error";
  return false;
}

bool BlockCirculantSolver::
HasNormInf() const
{
  // NormInf is not implemented
  return false;
}

const char* BlockCirculantSolver::
Label() const
{
   return label.c_str();
}

const Epetra_Comm& BlockCirculantSolver::
Comm() const
{
  return jacobian.Comm();
}

const Epetra_Map& BlockCirculantSolver::
OperatorDomainMap () const
{
  return jacobian.OperatorDomainMap();
}

const Epetra_Map& BlockCirculantSolver::
OperatorRangeMap () const
{
  return jacobian.OperatorRangeMap();
}

bool BlockCirculantSolver::
computePreconditioner(const Epetra_Vector& x,
		      Epetra_Operator& Prec,
		      Teuchos::ParameterList* p)
  {
  DEBUG("Enter BlockCirculantSolver::computePreconditioner");
  bool stat = true;

  DEBUG("construct matrix...");

  double argexp;
  double eigReal;
  double eigImag;

	for (int ii=0; ii<blocksOnDomain; ii++)
    {
    int i=globalComm->FirstTimeStepOnDomain()+ii;
    INFO("compute solver "<<i);
    argexp = i*2.0*3.14159265358979323846/numBlocks;
    eigReal = cos(argexp);//dt al in massBlock
    eigImag = -sin(argexp);
    complexProblem[ii]->UpdateValues(1.0,0.0,*jacobianBlock,eigReal,eigImag,*massBlock,*(lhsRe[ii]),*(lhsIm[ii]),*(rhsRe[ii]),*(rhsIm[ii]));

    DEBUG("symbolic...");
    Solver[ii]->SymbolicFactorization();
    DEBUG("numeric...");
    Solver[ii]->NumericFactorization();
    }

  DEBUG("Leave BlockCirculantSolver::computePreconditioner");
  return stat;
  }

void BlockCirculantSolver::
throwError(const std::string& functionName, const std::string& errorMsg) const
{
  std::cout << "BlockCirculantSolver" << functionName 
	        << " - " << errorMsg << std::endl;
  throw "LOCA Error";
}

Teuchos::RCP<const Epetra_Comm> BlockCirculantSolver::getTimeComm() const
  {
  return timeComm;
  }

int BlockCirculantSolver::PostProcSolution(Epetra_Vector& v) const
  {
  if (nullSpace_ != Teuchos::null)
    {
    for (int k=0;k<nullSpace_->NumVectors();k++)
      {
      double fac;

      EPETRA_CHK_ERR((*nullSpace_)(k)->Dot(v,&fac));
      EPETRA_CHK_ERR(v.Update(-fac,*((*nullSpace_)(k)),1.0));
      }
    }
  return 0;
  }


  void BlockCirculantSolver::buildTimeComm()
    {
    // create a communicator that communicates
    // between procs owning the same spatial part
    // but at different time-levels

    MPI_Comm old_comm = MPI_COMM_WORLD; //TODO: This is not very pretty, obviously!
    MPI_Comm new_comm;
    MPI_Group old_group, new_group;

    Epetra_Comm& sdComm = globalComm->SubDomainComm();

    int nproc_tot = globalComm->NumProc();
    int pid_tot = globalComm->MyPID();

    int nproc_xyz = sdComm.NumProc();
    int pid_xyz = sdComm.MyPID();
    int nproc_t   = nproc_tot/nproc_xyz;

    // will contain spatial rank of global proc i as element i
    int *xyz_ranks_ = new int[nproc_tot];
    int *xyz_ranks = new int[nproc_tot];

    for (int i=0;i<nproc_tot;i++) 
      {
      xyz_ranks_[i]=0;
      xyz_ranks[i]=0;
      }

    xyz_ranks_[pid_tot] = pid_xyz;
    globalComm->SumAll(xyz_ranks_,xyz_ranks, nproc_tot);

    int *row_ranks = new int[nproc_t];

    //  Create vector of ranks for processes owning same spatial part of grid:
    int pos=0;
    for (int k = 0; k < nproc_tot; k++)
      {
      if (xyz_ranks[k]==pid_xyz) row_ranks[pos++]=k;
      }

    /* ----------------------- */
    /* create the global group */
    /* ----------------------- */
    int ierr=MPI_Comm_group(old_comm, &old_group);

    /* --------------------- */
    /* Create the new group  */
    /* --------------------- */
    ierr=MPI_Group_incl(old_group, nproc_t, row_ranks, &new_group);
    if (ierr!=0)
      {
      std::string errorMessage = "MPI Error encountered!";
      throwError("LOCA::Epetra::xyztPrec::xyztPrec", errorMessage);
      }

    /* ---------------------------------------------------- */
    /* Create the new communicator; MPI assigns the context */
    /* ---------------------------------------------------- */

    ierr = MPI_Comm_create(old_comm, new_group, &new_comm);
    if (ierr!=0)
      {
      std::string errorMessage = "MPI Error encountered!";
      throwError("LOCA::Epetra::xyztPrec::xyztPrec", errorMessage);
      }

    // finally create the Epetra comm object
    timeComm = Teuchos::rcp(new Epetra_MpiComm(new_comm),false);

    delete [] xyz_ranks_;
    delete [] xyz_ranks;
    delete [] row_ranks;
    }

int BlockCirculantSolver::PutDirichlet(Epetra_CrsMatrix& A, int gid)
  {
  // find out which proc owns this row
  int lid, pid;

  EPETRA_CHK_ERR(A.RowMap().RemoteIDList(1,&gid,&pid,&lid));

  // find out how long that row is (how many nonzeros)
  int len;

  if (pid==A.Comm().MyPID())
    {
    EPETRA_CHK_ERR(A.NumMyRowEntries(lid,len));
    }

  EPETRA_CHK_ERR(A.Comm().Broadcast(&len,1,pid));

  int* indices=new int[len];
  double* values=new double[len];

  if (pid==A.Comm().MyPID())
    {
    int dummy_len;
    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(gid,len,dummy_len,values,indices));
    // set row to 0 and diagonal to 1
    for (int i=0;i<len;i++)
      {
      if (indices[i]==gid)
        {
        values[i]=1.0;
        }
      else
        {
        values[i]=0.0;
        }
      // put it back in
      EPETRA_CHK_ERR(A.ReplaceGlobalValues(gid,len,values,indices));
      }
    }

  // broadcast indices to everyone
  EPETRA_CHK_ERR(A.Comm().Broadcast(indices,len,pid));

  // we assume that the pattern of the matrix is symmetric and process all the rows in
  // indices, setting any coupling to gid to 0
  int *indices_i;
  double *values_i;
  int len_i;
  for (int i=0;i<len;i++)
    {
    int grid = indices[i];
    if (A.RowMap().MyGID(grid))
      {
      if (grid!=gid)
        {
        int lrid = A.LRID(grid);
        EPETRA_CHK_ERR(A.ExtractMyRowView(lrid,len_i,values_i,indices_i));
        for (int j=0;j<len_i;j++)
          {
          if (A.GCID(indices_i[j])==gid)
            {
            values_i[j]=0.0;
            }
          }
        }
      }
    }
  return 0;
  }

