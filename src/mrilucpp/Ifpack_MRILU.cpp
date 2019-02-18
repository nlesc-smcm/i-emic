/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#include "Ifpack_MRILU.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include <iomanip>
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "GlobalDefinitions.H"


// we still use one or two subroutines from the Fortran code,
// this might (and should) be changed in the future 
extern "C"
{
#ifdef HAVE_IFPACK_MRILU
	void mrilucpp_create(int* id, int* n, int* nnz, int* beg, int* jco, double* co);
	void mrilucpp_destroy(int* id);
	void mrilucpp_set_params(const int* id, int* blocksize, int* cutmck, int* scarow, 
							 int* xactelm,int* clsonce, 
							 double* nlsfctr, double* epsw, double* elmfctr, int* gusmod,
							 double* gusfctr,double* redfctr,
							 double* schtol, double* denslim, double* globfrac,
							 double* locfrac, double* sparslim, int* ilutype,
							 double* droptol, double* compfct, double* cpivtol,
							 double* lutol_, int* singlu, int* outlev);
	
	void mrilucpp_compute(const int* id);
	void mrilucpp_apply(const int* id, int* dim, const double *rhs, double   *sol);
#endif
}//extern


#ifdef DEBUGGING
#define DEBUG(s) std::cerr << comm->MyPID() << ": " << s << std::endl;
#else
#define DEBUG(s) 
#endif

#ifndef CHECK_ZERO
# define CHECK_ZERO(funcall) {int ierr = funcall;						\
		if (ierr) {std::cerr<<"Trilinos Error "<<ierr<<" returned from call "<<#funcall<<std::endl;}}
#endif


///////////////////////////////////////////////////////////////////////////////
// CLASS IFPACK_MRILU                                                          
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Constructor	                                                               
///////////////////////////////////////////////////////////////////////////////

// note: 'needs_setup' is currently irrelevant, the preconditioner is recomputed
//       completely every time.

Ifpack_MRILU::Ifpack_MRILU(Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::RCP<Epetra_Comm> comm_) : 
	mrilu_id(0),
	Matrix_(A),
	comm(comm_),
	is_initialized(false),is_computed(false)
{  
    std::string s1="MRILU(";
    std::string s2(A->Label());
    std::string s3=")";
	label = s1+s2+s3;
	DEBUG("Ifpack_MRILU constructor called");
	default_params();
}    //Constructor

Ifpack_MRILU::Ifpack_MRILU(Epetra_RowMatrix* A) :
	mrilu_id(0),
	is_initialized(false),is_computed(false)
{  
    std::string s1="MRILU(";
    std::string s2(A->Label());
    std::string s3=")";
	label = s1+s2+s3;
	Matrix_ = Teuchos::rcp(A,false);
	comm=Teuchos::rcp(&(A->Comm()),false);
	DEBUG("Ifpack_MRILU Ifpack constructor called");
	default_params();
	//std::cout << "PID " << comm->MyPID()<<", MRILU local matrix size = "<<A->NumMyRows()<<std::endl;
}    //Constructor

///////////////////////////////////////////////////////////////////////////////
// destructor                                                                  
///////////////////////////////////////////////////////////////////////////////

Ifpack_MRILU::~Ifpack_MRILU()
{
	DEBUG("Ifpack_MRILU destructor called");
	if (mrilu_id>0)
    {
		DEBUG("Destroy MRILU preconditioner...\n");
#ifdef HAVE_IFPACK_MRILU
		mrilucpp_destroy(&mrilu_id);
#endif
    }    
	is_initialized=false;
	is_computed=false;
}

//////////////

const Epetra_Comm& Ifpack_MRILU::Comm() const {return *comm;}

const Epetra_Map& Ifpack_MRILU::OperatorDomainMap() const
{return Matrix_->OperatorDomainMap();}

const Epetra_Map& Ifpack_MRILU::OperatorRangeMap() const
{return Matrix_->OperatorRangeMap();}


int Ifpack_MRILU::SetParameters(Teuchos::ParameterList& ifpParams)
{
	Teuchos::ParameterList& lsParams = ifpParams.sublist("MRILU");
	DEBUG("Parameters supplied by User: ");
	DEBUG(lsParams);

	blocksize =  lsParams.get("blocksize",blocksize);
	cutmck    =  lsParams.get("cutmck",cutmck);
	scarow    =  lsParams.get("scarow",scarow);
	xactelm   =  lsParams.get("xactelm",xactelm);
	clsonce   =  lsParams.get("clsonce",clsonce);
	nlsfctr   =  lsParams.get("nlsfctr",nlsfctr);
	epsw      =  lsParams.get("epsw",epsw);
	elmfctr   =  lsParams.get("elmfctr",elmfctr);
	gusmod    =  lsParams.get("gusmod",gusmod);
	gusfctr   =  lsParams.get("gusfctr",gusfctr);
	redfctr   =  lsParams.get("redfctr",redfctr);
	schtol    =  lsParams.get("schtol",schtol);
	denslim   =  lsParams.get("denslim", denslim);
	globfrac  =  lsParams.get("globfrac", globfrac);
	locfrac   =  lsParams.get("locfrac", locfrac);
	sparslim  =  lsParams.get("sparslim", sparslim);
	ilutype   =  lsParams.get("ilutype", ilutype);
	droptol   =  lsParams.get("droptol", droptol);
	compfct   =  lsParams.get("compfct", compfct);
	cpivtol   =  lsParams.get("cpivtol", cpivtol);
	lutol     =  lsParams.get("lutol", lutol);
	singlu    =  lsParams.get("singlu", singlu);
	outlev    =  lsParams.get("Output Level",outlev);

	DEBUG("Parameters used: ");
	DEBUG(lsParams);
	return 0;
}



// build preconditioner. This function is called by NOX/LOCA, 
// I am not sure what the arguments mean but we do not use them
// anyway
bool Ifpack_MRILU::computePreconditioner(const Epetra_Vector& x,
										 Epetra_Operator& Prec,
										 Teuchos::ParameterList* p)
{
	DEBUG("+++ Enter Ifpack_MRILU::computePreconditioner");
	if (p!=NULL) this->SetParameters(*p);
	if (!IsInitialized())
    {
		this->Initialize();
    }
	this->Compute();
	DEBUG("--- Leave Ifpack_MRILU::computePreconditioner");
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// apply inverse preconditioner step                                           
///////////////////////////////////////////////////////////////////////////////

int Ifpack_MRILU::
ApplyInverse(const Epetra_MultiVector& input,
			 Epetra_MultiVector& result) const
{
	if (!IsComputed())
    {
		Error("Ifpack_MRILU not yet computed!",__FILE__,__LINE__);
    }
	if (&input == &result)
		Error("aliased call to Ifpack_MRILU::ApplyInverse",__FILE__,__LINE__);
	
	DEBUG("+++ Enter Ifpack_MRILU::ApplyInverse");
	if (input.NumVectors()!=1)
    {
		this->Error("Multiple Rhs not supported by Ifpack_MRILU!",__FILE__,__LINE__);
    }
    
#ifdef HAVE_IFPACK_MRILU
	int n = input.MyLength();
	double *sol_array; //= new double[n];
	const double *rhs_array; //= new double[n];
	//CHECK_ZERO(input.ExtractCopy(rhs_array,n));
	
	rhs_array = input[0];
	sol_array = result[0];
	DEBUG("Apply MRILU preconditioner...");
/*
  if (outlev>3)
  {
  std::cout << "Apply " << label << std::endl;
  }
*/
	if (is_identity)
		result = input;
	else
	{
		mrilucpp_apply(&mrilu_id, &n, rhs_array,sol_array);
	}
	
	//for (int i=0;i<n;i++)    result[0][i]=sol_array[i];
	//delete [] sol_array;
	//delete [] rhs_array;
#else
	std::cout << "WARNING: MRILU is not available, using identity preconditioner."<<std::endl;
	result=input;
#endif
	DEBUG("done!");

	
	DEBUG("+++ Leave Ifpack_MRILU::ApplyInverse");
	return 0;
	
}



///////////////////////////////////////////////////////////////////////////////
// Apply preconditioner matrix (not available)                                 
///////////////////////////////////////////////////////////////////////////////

int Ifpack_MRILU::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
	this->Error("Ifpack_MRILU::Apply() - method is NOT implemented!!  ",__FILE__,__LINE__);
	return false;
}


///////////////////////////////////////////////////////////////////////////////
// inf-norm (n/a)                                                              
///////////////////////////////////////////////////////////////////////////////


double Ifpack_MRILU::NormInf() const
{
	this->Error("Ifpack_MRILU::NormInf() - method is NOT implemented!!  ",__FILE__,__LINE__);
	return -1.0;
}




//////////////////////////////////////////////////////////////////////////////
// remaining IFPACK interface                                               //
//////////////////////////////////////////////////////////////////////////////

//! Computes all it is necessary to initialize the preconditioner.
int Ifpack_MRILU::Initialize()
{
	DEBUG("+++ Enter Ifpack_MRILU::Initialize");

	// extract CRS array of the local Jacobian.
	int nnz   = Matrix_->NumMyNonzeros();
	int nrows = Matrix_->NumMyRows();
  
	// although we might get CRS arrays directly from Epetra,
	// this method is preferred because it is implementation independent
	int* beg   = new int[nrows+1];
	int* jco   = new int[nnz];
	double* co = new double[nnz];
  
	int len;
	beg[0]  = 0;
	int idx = 0;
	is_identity = true;

	for (int i = 0; i < nrows; i++)
    {
		CHECK_ZERO(Matrix_->ExtractMyRowCopy(i, nnz-beg[i], len, co+beg[i], jco+beg[i]));
		beg[i+1] = beg[i] + len;
		
		// Check whether this matrix is the identity
		if (is_identity)
		{
			for (int j = beg[i]; j != beg[i+1]; ++j)
			{
				if ( ((jco[idx] == i) && (std::abs(co[idx] - 1.0) > 1e-7)) ||
					 ((jco[idx] != i) && (std::abs(co[idx])       > 1e-7))   )
				{
					is_identity = false;
					break;
				}
				idx++;
			}
		}
    }

	DEBUG("Create preconditioner...");	

#ifdef HAVE_IFPACK_MRILU
	// the loca 'Reuse Policy' doesn't work for user
	// defined preconditioners, it seems. They just 
	// call this function whenever they want a new precond.
	mrilucpp_destroy(&mrilu_id);
	// note that the conversion to 1-based indexing is done in fortran:
	mrilucpp_create(&mrilu_id, &nrows, &nnz, beg, jco, co);
#endif    

	delete [] beg;
	delete [] co;
	delete [] jco;
	is_initialized = true;
	DEBUG("+++ Leave Ifpack_MRILU::Initialize");
	return 0;
}

//! Returns true if the  preconditioner has been successfully initialized, false otherwise.
bool Ifpack_MRILU::IsInitialized() const
{
    return is_initialized;
}

//! Computes all it is necessary to apply the preconditioner.
int Ifpack_MRILU::Compute()
{
    DEBUG("Enter Ifpack_MRILU::Compute");

#ifdef HAVE_IFPACK_MRILU

	//note: needs_setup is ignored up to now

	if (comm->NumProc()>1)
    {
		// this class should be run through an Ifpack_LocalFilter or the like:
		this->Error("The MRILU preconditioner is not intended for parallel use!",__FILE__,__LINE__);
    }

	if (!this->IsInitialized()) this->Initialize();

	DEBUG("set parameters...");
	mrilucpp_set_params( &mrilu_id, &blocksize, &cutmck ,  &scarow ,  &xactelm ,  
						 &clsonce,  &nlsfctr,  
						 &epsw   ,  &elmfctr,  &gusmod  ,  &gusfctr,  &redfctr,
						 &schtol ,  &denslim,  &globfrac,  &locfrac,  &sparslim,
						 &ilutype,  &droptol,  &compfct ,  &cpivtol,  &lutol,  
						 &singlu ,  &outlev);

	DEBUG("Compute factorization...");
	if (outlev>3)
    {
		std::cout << "Compute " << label << std::endl;
    }
	mrilucpp_compute(&mrilu_id);  
	DEBUG("done!");
  
	// after building the preconditioner, the internal csr matrix is destroyed
	is_computed=true;
	is_initialized=false; // always have to re-initialize before calling Compute again!

#else
	std::cout << "WARNING: MRILU is not available, using identity!"<<std::endl;
#endif    
    DEBUG("Leave Ifpack_MRILU::Compute");
    return 0;
}

//! Returns true if the  preconditioner has been successfully computed, false otherwise.
bool Ifpack_MRILU::IsComputed() const
{
    return is_computed;
}

//! Computes the condition number estimate, returns its value.
double Ifpack_MRILU::Condest(const Ifpack_CondestType CT,
							 const int MaxIters,
							 const double Tol,
							 Epetra_RowMatrix* Matrix)
{
    // this is not yet implemented, although it should be trivial to get
    // some estimate from the ILU factorization
    condest=-1.0;
    return condest;
}

//! Returns the computed condition number estimate, or -1.0 if not computed.
double Ifpack_MRILU::Condest() const
{
    return condest;
}

//! Returns a pointer to the matrix to be preconditioned.
const Epetra_RowMatrix& Ifpack_MRILU::Matrix() const 
{
    if (Matrix_==Teuchos::null)
	{
		Error("Matrix handle not available!",__FILE__,__LINE__);
	}
    return *Matrix_;
}

// TODO: none of the performance measuring routines below are implemented, yet

//! Returns the number of calls to Initialize().
int Ifpack_MRILU::NumInitialize() const
{
    return 0;
}

//! Returns the number of calls to Compute().
int Ifpack_MRILU::NumCompute() const
{
    return 0;
}

//! Returns the number of calls to ApplyInverse().
int Ifpack_MRILU::NumApplyInverse() const
{
    return 0;
}

//! Returns the time spent in Initialize().
double Ifpack_MRILU::InitializeTime() const
{
    return -1.0;
}

//! Returns the time spent in Compute().
double Ifpack_MRILU::ComputeTime() const
{
    return -1.0;
}

//! Returns the time spent in ApplyInverse().
double Ifpack_MRILU::ApplyInverseTime() const
{
    return -1.0;
}


//! Returns the number of flops in the initialization phase.
double Ifpack_MRILU::InitializeFlops() const
{
    return 0.0;
}

//! Returns the number of flops in the computation phase.
double Ifpack_MRILU::ComputeFlops() const
{
    return 0.0;
}

//! Returns the number of flops in the application of the preconditioner.
double Ifpack_MRILU::ApplyInverseFlops() const
{
    return 0.0;
}

//! Prints basic information on iostream. This function is used by operator<<.
std::ostream& Ifpack_MRILU::Print(std::ostream& os) const
{
    // TODO: print something interesting for the user to read and enjoy
    os << this->Label();
    return os;
}

void Ifpack_MRILU::default_params()
{
	blocksize=1;
	cutmck=false;
	scarow  =  true;
	xactelm = false;
	clsonce = false;
	nlsfctr  = 0.1;
	epsw     = 1e-2;
	elmfctr  = 0.1;
	gusmod  =  false;
	gusfctr  = 1.0;
	redfctr  = 0.8;
	schtol   = 0.0;
	denslim  = 0.01;
	globfrac = 0.1;
	locfrac  = 1.0;
	sparslim = 0.9;
	ilutype  = 9;
	droptol  = 1e-4;
	compfct  = 1.0;
	cpivtol  = 0.5;
	lutol   = 1e-10;
	singlu  = false;
	outlev  = 2;
}

void Ifpack_MRILU::Error(std::string msg, std::string file, int line) const
{
	std::cerr << "ERROR: "<<msg<<std::endl;
	std::cerr <<         "("<<file<<", line "<<line<<")\n";
#ifdef HAVE_MPI
	MPI_Abort(MPI_COMM_WORLD,-1);
#endif
	exit(-1);
}

