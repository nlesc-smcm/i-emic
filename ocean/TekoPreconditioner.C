#include "TekoPreconditioner.H"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactoryOperator.hpp"
#include "Teko_BlockedEpetraOperator.hpp"
#include "Teko_StridedEpetraOperator.hpp"

#include <EpetraExt_RowMatrixOut.h>

// constructor
TekoPreconditioner::TekoPreconditioner(Teuchos::RCP<const Epetra_RowMatrix> K, 
									   Teuchos::RCP<Teuchos::ParameterList> params)
	:
	numInitialize_(0),
	numCompute_(0),
	numApplyInverse_(0),
	flopsInitialize_(0.0),
	flopsCompute_(0.0),
	flopsApplyInverse_(0.0),
	timeInitialize_(0.0),
	timeCompute_(0.0),
	timeApplyInverse_(0.0),
	initialized_(false),
	computed_(false),
	matrix_(K),
	params_(params),
	comm_(Teuchos::rcp(&(K->Comm()),false)),
	normInf_(-1.0),
	useTranspose_(false),
	label_("TekoPreconditioner")
{
	std::cout << "TekoPreconditioner constructor called" << '\n';
}

// Ifpack_TekoPreconditioner interface

// Computes all it is necessary to initialize the preconditioner.
int TekoPreconditioner::Initialize()
{
	std::cout << "entering TekoPreconditioner::Initialize" << std::endl;

    Teuchos::RCP<Teko::InverseLibrary> invLib =
		Teko::InverseLibrary::buildFromParameterList(*params_);
	Teuchos::RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("GS-Outer");

	// Create the initial preconditioner, and build it from strided_A
	prec_ = Teuchos::rcp(new Teko::Epetra::InverseFactoryOperator(inverse));
	prec_->initInverse();
	
	initialized_ = true;
	numInitialize_++;

	std::cout << "leaving TekoPreconditioner::Initialize" << std::endl;
	return 0;
}

// Returns true if the  preconditioner has been successfully initialized, false otherwise.
bool TekoPreconditioner::IsInitialized() const {return initialized_;}

// Computes all it is necessary to apply the preconditioner.
int TekoPreconditioner::Compute()
{
	std::cout << "entering TekoPreconditioner::Compute" << std::endl;
    if (!IsInitialized())
	{
		// the user should normally call Initialize before Compute
		std::cout << "TekoPreconditioner not initialized. I'll do it for you."
				  << __FILE__ << __LINE__ << std::endl;
		this->Initialize();
	}

	int n = matrix_->NumMyRows() / 6;

	std::vector< std::vector<int> > vec(4);
    vec[0] = std::vector<int>(n*2);
    vec[1] = std::vector<int>(n);
	vec[2] = std::vector<int>(n);
	vec[3] = std::vector<int>(n*2);
    Epetra_Map const &map = matrix_->RowMatrixRowMap();
    for (int i = 0; i != n; ++i)
    {
        vec[0][i*2]    =  map.GID(i*6);   // u
		vec[0][i*2+1]  =  map.GID(i*6+1); // v
		vec[1][i]      =  map.GID(i*6+2); // w
		vec[2][i]      =  map.GID(i*6+3); // p
		vec[3][i*2]    =  map.GID(i*6+4); // T
		vec[3][i*2+1]  =  map.GID(i*6+5); // S
    }
	
	Teuchos::RCP<Teko::Epetra::BlockedEpetraOperator> op =
		Teuchos::rcp(new Teko::Epetra::BlockedEpetraOperator(vec,matrix_));
    //==================================================================
	// EXPERIMENTAL
	//==================================================================
	/*
	Teuchos::RCP<const Epetra_CrsMatrix> Block11 =
		Teuchos::rcp_const_cast<const Epetra_CrsMatrix>(op->GetBlock(1,1));
	
	int NumMyElements = Block11->Map().NumMyElements();
	int *MyGlobalElements = Block11->Map().MyGlobalElements();
	int index;
	double diag = 1.0;
    //asserten			
	for (int i = 0; i != NumMyElements; ++i)
	{
		Block11->ReplaceGlobalValues(MyGlobalElements[i], 1,
								  &diag, MyGlobalElements + i);
	}
	Block11->FillComplete();
	*/
    //==================================================================
	// /EXPERIMENTAL
	//==================================================================	
	
	std::string reorder = "[[0 1 2] 3]";
	Teuchos::RCP<const Teko::BlockReorderManager> brm =
		Teko::blockedReorderFromString(reorder);
	op->Reorder(*brm);
    prec_->buildInverseOperator(Teuchos::rcp_static_cast<const Epetra_Operator>(op.getConst()));
	
	computed_ = true;

	numCompute_++;
	std::cout << "leaving TekoPreconditioner::Compute" << std::endl;
	return 0;
}

// Sets all parameters for the preconditioner.
int TekoPreconditioner::SetParameters(Teuchos::ParameterList& List)
{
	std::cout << "NOT IMPLEMENTED" << std::endl;
	throw 1;
}

// Returns true if the  preconditioner has been successfully computed, false otherwise.
bool TekoPreconditioner::IsComputed() const {return computed_;}

// Applies the preconditioner to vector X, returns the result in Y.
int TekoPreconditioner::ApplyInverse(const Epetra_MultiVector& B,
									 Epetra_MultiVector& X) const
{
    numApplyInverse_++;
//    time_->ResetStartTime();
    prec_->ApplyInverse(B, X);
//    timeApplyInverse_+=time_->ElapsedTime();
    return 0;
}

// Returns a pointer to the matrix to be preconditioned.
const Epetra_RowMatrix& TekoPreconditioner::Matrix() const {return *matrix_;}

//TODO: the flops-counters currently do not include anything inside Belos

double TekoPreconditioner::InitializeFlops() const
{
    return 0;
}

double TekoPreconditioner::ComputeFlops() const
{
    return 0;
}

double TekoPreconditioner::ApplyInverseFlops() const
{
    return 0;
}


// Computes the condition number estimate, returns its value.
double TekoPreconditioner::Condest(const Ifpack_CondestType CT,
								   const int MaxIters,
								   const double Tol,
								   Epetra_RowMatrix* Matrix)
{
	return -1.0; // not implemented.
}

// Returns the computed condition number estimate, or -1.0 if not computed.
double TekoPreconditioner::Condest() const
{
    return -1.0;
}


// Returns the number of calls to Initialize().
int TekoPreconditioner::NumInitialize() const {return numInitialize_;}

// Returns the number of calls to Compute().
int TekoPreconditioner::NumCompute() const {return numCompute_;}

// Returns the number of calls to ApplyInverse().
int TekoPreconditioner::NumApplyInverse() const {return numApplyInverse_;}

// Returns the time spent in Initialize().
double TekoPreconditioner::InitializeTime() const {return timeInitialize_;}

// Returns the time spent in Compute().
double TekoPreconditioner::ComputeTime() const {return timeCompute_;}

// Returns the time spent in ApplyInverse().
double TekoPreconditioner::ApplyInverseTime() const {return timeApplyInverse_;}
