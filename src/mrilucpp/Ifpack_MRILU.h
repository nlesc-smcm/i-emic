/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
/*
  This is a light-weight interface so that MRILU (by Fred Wubs) can be
  called from Trilinos. The interface is such that parameters can be
  read from .xml files and the solvers look to Trilinos like
  sequential Ifpack_Preconditioners, so that they can be run in
  parallel using Additive Schwarz.
*/

/*
  To be able to use this with code using a PreconditionerFactory, Trilinos needs to be patched.
  Ifpack.h and Ifpack.cpp in Trilinos need to be adjusted, see patched examples in notes.
*/

#ifndef IFPACK_MRILU_H
#define IFPACK_MRILU_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Ifpack_Preconditioner.h"


class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Comm;
class Epetra_Operator;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;


//! wrapper for the MRILU preconditioner
class Ifpack_MRILU : public Ifpack_Preconditioner
{
      
public:
      
	//! Constructor
	Ifpack_MRILU(Teuchos::RCP<Epetra_CrsMatrix> A,Teuchos::RCP<Epetra_Comm> comm);

	//! this Constructor is needed for Ifpack_AdditiveScharz
	Ifpack_MRILU(Epetra_RowMatrix* A);
      
	//! Destructor
	virtual ~Ifpack_MRILU();
      
	//! Set transpose (not implemented => returns -1).
	int SetUseTranspose(bool UseTranspose) 
	{
        // return 0 (success) if transpose not requested,
        // otherwise -1 (transpose is not available)
        return -(int)UseTranspose;
	}
      
	//! Apply MRILU preconditioning operator (not implemented)
	int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
      
	//! Apply preconditioner operator inverse
	int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
      
	//! Computing infinity norm (not implemented)
	double NormInf() const;
      
	//! Label
	const char* Label() const {return label.c_str();}
      
	//! Transpose
	bool UseTranspose() const {return false;}
      
	//! Have norm-inf
	bool HasNormInf() const {return false;}
      
	/*! 
	 * \brief Returns a pointer to the Epetra_Comm communicator associated 
	 * with this operator.
	 */
	const Epetra_Comm& Comm() const;
      
	/*! 
	 * \brief Returns the Epetra_Map object associated with the domain of 
	 * this operator.
	 */
	const Epetra_Map& OperatorDomainMap() const;
      
	/*! 
	 * \brief Returns the Epetra_Map object associated with the range of 
	 * this operator.
	 */
	const Epetra_Map& OperatorRangeMap() const;
    
	//! Compute preconditioner \f$M\f$. 
      
	//! The list of valid parameters is formed by the original MRILU names like
	//! "epsw", "droptol", "cutmck", etc. The only exception is the output level (outlev)
	//! which is now called "Output Level". All parameters are either double or int, no logicals/bools
	//! are accepted (use 0/1 ints instead).
	bool computePreconditioner(const Epetra_Vector& x,
							   Epetra_Operator& Prec,
							   Teuchos::ParameterList* p = NULL);
      

    //! \name Ifpack_Preconditioner interface
    
    //@{

	//! set parameters
	int SetParameters(Teuchos::ParameterList& List);

	//! Computes all it is necessary to initialize the preconditioner.
	int Initialize();

	//! Returns true if the  preconditioner has been successfully initialized, false otherwise.
	bool IsInitialized() const;

	//! Computes all it is necessary to apply the preconditioner.
	int Compute();

	//! Returns true if the  preconditioner has been successfully computed, false otherwise.
	bool IsComputed() const;

	//! Computes the condition number estimate, returns its value.
	double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
				   const int MaxIters = 1550,
				   const double Tol = 1e-9,
				   Epetra_RowMatrix* Matrix_ = 0);

	//! Returns the computed condition number estimate, or -1.0 if not computed.
	double Condest() const;

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
        std::ostream& Print(std::ostream& os) const;

    
    //@}

      
protected:
                  
	//! Label for this operator.
	std::string label;
      
	//! identifies preconditioner instance in fortran
      
	//! during the first 'Initialize' call, we obtain a unique
	//! identifier from the fortran module m_mriluprec, which 
	//! stays the same ever after and is passed to all fortran
	//! calls.
	int mrilu_id;
    
	//! the matrix of the linear system we want to precondition
	Teuchos::RCP<Epetra_RowMatrix> Matrix_;

	//! communicator
	Teuchos::RCP<const Epetra_Comm> comm;
      
	//! matrix passed to MRILU?
	bool is_initialized;

	//! mrilu computed
	bool is_computed;
      
	//! condition number estimate
	double condest;

	//! if our matrix forms the identity
	bool is_identity;
	//! \name MRILU parameters
	//!@{

// TODO: documentation of parameters, 
//       reasonable default values
          
	int  blocksize;
	int cutmck;
	int scarow;
	int xactelm;
	int clsonce;
	double nlsfctr;
	double epsw;
	double elmfctr;
	int gusmod;
	double gusfctr;
	double redfctr;
	double schtol;

	double denslim;
	double globfrac;
	double locfrac;
	double sparslim;

	int ilutype;
	double droptol;
	double compfct;
	double cpivtol;
	double lutol;
	int singlu;

	//!@}

	int outlev; //! output level (0-5)
           
private:      
            
	//! set default values for MRILU params
	void default_params();      
  
	//! Error function. We don't use the one from Trilinos-THCM
	//! to avoid a dependency on its Filestreams/globdefs      
	void Error(std::string msg, std::string file, int line) const;

      
}; // end of class Ifpack_MRILU
    


#endif // MRILU_PREC_H
