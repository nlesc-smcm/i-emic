#ifndef TEKO_PRECONDITIONER_H
#define TEKO_PRECONDITIONER_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Ifpack_Preconditioner.h"
#include "Ifpack_Container.h"
#include "Epetra_Time.h"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosOutputManager.hpp"

#include "Teko_InverseFactoryOperator.hpp"
#include "Teko_EpetraBlockPreconditioner.hpp"
#include "Teko_BlockedEpetraOperator.hpp"
#include "Teko_StridedEpetraOperator.hpp"

// forward declarations
class Epetra_Comm;
class Epetra_Map;
class Epetra_RowMatrix;
class Epetra_Import;

namespace Teuchos 
{
	class ParameterList;
}

//! main solver class in HYMLS

/*! This class 
  - sets parameters for the problem
  - creates an overlapping partitioning (HYMLS::OverlappingPartitioner)
  - computes a factorization of the subdomain matrices
  - creates a HYMLS::SchurComplement
  - creates a HYMLS::SchurTekoPreconditioner for the SC
  - creates a Belos iterative solver
  - solves the problem (in ApplyInverse())
  
  The solver has the following features:
  
  - recursive application: the SchurTekoPreconditioner object will check
  the parameter "Number of Levels" and see if it should function as 
  a direct or approximate solver. For instance, if you set "Number  
  of Levels" to 1, the solver does the domain decomposition and then
  solves the separator problem directly. Setting it to 2 gives the  
  standard two-level method where the preconditioner does orthogonal
  transformation and dropping and creates a direct solver for the   
  Vsum problem. For more levels, it will create another instance of 
  this "TekoPreconditioner" class.                                      
  
  - monitoring and deflating the most unstable modes: to activate   
  this, set the "Deflated Subspace Dimension" to a positive integer.
  The solver computes the eigenvalues closest to 0 of the precondi- 
  tioner. If an eigenvalue has croosed the imaginary axis, a message
  is issued. That's all for now, the deflation still has to be      
  implemented.
  
*/
class TekoPreconditioner : public Ifpack_Preconditioner

{
  
public:


	//!                                                     
	//! Constructor                                         
	//!                                                     
	//! the caller should typically just use                
	//!                                                     
	//! TekoPreconditioner(K,params);                           
	//!                                                     
	//! The constructor with pre-constructed partitioner is 
	//! used for recursive application of the method        
	//! (SchurTekoPreconditioner will construct a TekoPreconditioner
	//  object based on an OverlappingPartitioner created by
	//! SpawnNextLevel())                                   
	//!                                                     
	TekoPreconditioner(Teuchos::RCP<const Epetra_RowMatrix> K, 
					   Teuchos::RCP<Teuchos::ParameterList> params);

	//! destructor
	~TekoPreconditioner(){};

	//! Computes all it is necessary to initialize the preconditioner.
	int Initialize();

	//! Returns true if the  preconditioner has been successfully initialized, false otherwise.
	bool IsInitialized() const;

	int SetParameters(Teuchos::ParameterList& List);
	
	//! Computes all it is necessary to apply the preconditioner.
	int Compute();

	//! Returns true if the  preconditioner has been successfully computed, false otherwise.
	bool IsComputed() const;


	//! Computes the condition number estimate, returns its value.
	double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
				   const int MaxIters = 1550,
				   const double Tol = 1e-9,
				   Epetra_RowMatrix* Matrix = 0);

	//! Returns the computed condition number estimate, or -1.0 if not computed.
	double Condest() const;

	//! Applies the operator (not implemented)
	int Apply(const Epetra_MultiVector& X,
			  Epetra_MultiVector& Y) const {return -1;}

	//! Applies the preconditioner to vector X, returns the result in Y.
	int ApplyInverse(const Epetra_MultiVector& X,
					 Epetra_MultiVector& Y) const;

	//! Returns a pointer to the matrix to be preconditioned.
	const Epetra_RowMatrix& Matrix() const;

	//! Returns the number of calls to Initialize().
	int NumInitialize() const;

	//! Returns the number of calls to Compute().
	int NumCompute() const;

	//! Returns the number of calls to ApplyInverse().
	int NumApplyInverse() const;

	//! Returns the time spent in Initialize().
	double InitializeTime() const;

	//! Returns the time spent in Compute().
	double ComputeTime() const;

	//! Returns the time spent in ApplyInverse().
	double ApplyInverseTime() const;

	//! Returns the number of flops in the initialization phase.
	double InitializeFlops() const;

	//! Returns the number of flops in the computation phase.
	double ComputeFlops() const;

	//! Returns the number of flops in the application of the preconditioner.
	double ApplyInverseFlops() const;

	//! Prints basic information on iostream. This function is used by operator<<.
	ostream& Print(std::ostream& os) const {};

	int SetUseTranspose(bool UseTranspose)
		{
			useTranspose_=false; // not implemented.
			return -1;
		}
	//! not implemented.
	bool HasNormInf() const {return false;}

	//! infinity norm
	double NormInf() const {return normInf_;}

	//! label
	const char* Label() const {return label_.c_str();}
  
	//! use transpose?
	bool UseTranspose() const {return useTranspose_;}
  
	//! communicator
	const Epetra_Comm & Comm() const {return *comm_;}
  
    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return *rangeMap_;}

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const {return *rangeMap_;}
  
	//@}
   

protected:
	//! parameters
	Teuchos::RCP<const Teuchos::ParameterList> params_;
	//! communicator
	Teuchos::RCP<const Epetra_Comm> comm_;

	//! fake communicator
	Teuchos::RCP<const Epetra_Comm> serialComm_;
  
	//! matrix based on range map
	Teuchos::RCP<const Epetra_RowMatrix> matrix_;
  
	//! internal reordered matrix (first interior and then separator vars)
	Teuchos::RCP<Epetra_CrsMatrix> reorderedMatrix_;

	//! range/domain map of matrix_
	Teuchos::RCP<const Epetra_Map> rangeMap_;

	//! row map of this operator (first all interior and then all separator variables, no overlap)
	Teuchos::RCP<const Epetra_Map> rowMap_;
  
	//! importer from range to row map
	Teuchos::RCP<Epetra_Import> importer_;

	//~ Teuchos::RCP<Teko::Epetra::EpetraBlockPreconditioner> prec_;
	Teuchos::RCP<Teko::Epetra::InverseFactoryOperator> prec_;
	//@}
  
	//! use transposed operator?
	bool useTranspose_;
  
	//! infinity norm
	double normInf_;
  
	//! label
	std::string label_;
  
	//! timer
	mutable Teuchos::RCP<Epetra_Time> time_;
  
	//! has Initialize() been called?
	bool initialized_;

	//! has Compute() been called?
	bool computed_;
  
	//! true if SetBorder() has been called with non-null args
	bool haveBorder_;
  
	//! how often has Initialize() been called?
	int numInitialize_;

	//! how often has Compute() been called?
	int numCompute_;

	//! how often has ApplyInverse() been called?
	mutable int numApplyInverse_;
  
	//! flops during Initialize()
	double flopsInitialize_;

	//! flops during Compute()
	double flopsCompute_;

	//! flops during ApplyInverse()
	mutable double flopsApplyInverse_;

	//! time during Initialize()
	mutable double timeInitialize_;

	//! time during Compute()
	mutable double timeCompute_;

	//! time during ApplyInverse()
	mutable double timeApplyInverse_;
  
	//!@}
  
  
	//! geometric info, used for storing the ordering in MATLAB for visualization
	int dim_,dof_,nx_,ny_,nz_;
};

#endif /* TEKO_PRECONDITIONER_H */
