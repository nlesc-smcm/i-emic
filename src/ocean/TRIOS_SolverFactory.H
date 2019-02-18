/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef TRIOS_SOLVERFACTORY_H
#define TRIOS_SOLVERFACTORY_H

#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_Preconditioner.h"

// typedef'd Teuchos pointers

class AztecOO;
class OceanModel;

namespace TRIOS {

  //! class with static functions for creating krylov solvers and
  //! preconditioners of all sorts.
  class SolverFactory
    {

    public:


      //! create a preconditinoer for a matrix
      //! verbose=5 doesn't change anything
      //! verbose=0 makes the preconditioner silent
      //! verbose=10 makes it chatter
      static Teuchos::RCP<Epetra_Operator> CreateAlgebraicPrecond(Epetra_CrsMatrix& A, Teuchos::ParameterList& plist,
                         int verbose=5);

      //! compute preconditinoer for a matrix
      static void ComputeAlgebraicPrecond(Teuchos::RCP<Epetra_Operator> P, Teuchos::ParameterList& plist);

      //! create a preconditinoer for a matrix
      //! verbose=5 doesn't change anything
      //! verbose=0 makes the solver silent
      //! verbose=10 makes it chatter
      static Teuchos::RCP<AztecOO> CreateKrylovSolver(Teuchos::ParameterList& plist,int verbose=5);

      //! convert parameterlist to Aztec options array
      static void ExtractAztecOptions(Teuchos::ParameterList& list, int* options, double* params);

      //! do eigen-analysis of operator P\A-I to analyze preconditioner properties
      static void AnalyzeSpectrum(Teuchos::ParameterList& plist, Teuchos::RCP<const Epetra_Operator> Prec);


    private:

      //! stores all levels of an ML hierarchy in files
      static void DumpMLHierarchy(Teuchos::RCP<Epetra_Operator> Prec);

    }; // end of class SolverFactory


//////////////////////////////////////////////////////////////////////


//! identity operator (used if "None" is selected as algebraic precond for
//! one of the subsystems)
class IdentityOperator : public Ifpack_Preconditioner
  {

 public:
//-----//


  //! constructor
  IdentityOperator(const Epetra_Map& RangeMap,const Epetra_Map& DomainMap,
        const Epetra_Comm& comm_)
    : comm(comm_), rangeMap(RangeMap), domainMap(DomainMap), useTranspose(false) {}

  //! destructor
  virtual ~IdentityOperator(){}

//!\name  Atribute set methods
//!@{

  //! If set true, transpose of this operator will be applied (n/a)
  int 	SetUseTranspose (bool UseTranspose) {useTranspose=UseTranspose; return 0;}

//!@}

//!\name Mathematical functions
//!@{

  //! apply operator Y=Op*X
  int 	Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
    Y=X;
    return 0;
    }

  //! apply inverse operator (n/a)
  int 	ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
    Y = X;
    return 0;
    }

  //! Returns the infinity norm of the global matrix.
  double NormInf () const {return 1.0;}


//!@}


//!\name Atribute access functions
//!@{

  //!Returns a character string describing the operator.
  const char * 	Label () const {return "Identity";}

  //! Returns the current UseTranspose setting.
  bool 	UseTranspose () const {return useTranspose;}

  //! Returns true if the this object can provide an approximate Inf-norm, false otherwise.
  bool 	HasNormInf () const {return true;}

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & 	Comm () const {return comm;}

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & 	OperatorDomainMap () const
    {
    return domainMap;
    }

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & 	OperatorRangeMap () const
    {
    return rangeMap;
    }

  //@}

  //! \name Ifpack_Preconditioner interface
  //@{
  int SetParameters(Teuchos::ParameterList& List) {return 0;}

  //! Computes all it is necessary to initialize the preconditioner.
  int Initialize() {return 0;}

  //! Returns true if the  preconditioner has been successfully initialized, false otherwise.
  bool IsInitialized() const {return true;}

  //! Computes all it is necessary to apply the preconditioner.
  int Compute() {return 0;}

  //! Returns true if the  preconditioner has been successfully computed, false otherwise.
  bool IsComputed() const {return true;}

  //! Computes the condition number estimate, returns its value.
  double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
                         Epetra_RowMatrix* Matrix = 0)
  {return -1.0;}

  //! Returns the computed condition number estimate, or -1.0 if not computed.
  double Condest() const {return -1.0;}


  //! Returns a pointer to the matrix to be preconditioned.
  const Epetra_RowMatrix& Matrix() const {throw("IdentityOperator::Matrix() not implemented!");}

  //! Returns the number of calls to Initialize().
  int NumInitialize() const {return -1;}

  //! Returns the number of calls to Compute().
  int NumCompute() const {return -1;}

  //! Returns the number of calls to ApplyInverse().
  int NumApplyInverse() const {return -1;}
  //! Returns the time spent in Initialize().
  double InitializeTime() const {return 0.0;}

  //! Returns the time spent in Compute().
  double ComputeTime() const {return 0.0;}

  //! Returns the time spent in ApplyInverse().
  double ApplyInverseTime() const {return 0.0;}

  //! Returns the number of flops in the initialization phase.
  double InitializeFlops() const {return 0;}

  //! Returns the number of flops in the computation phase.
  double ComputeFlops() const {return 0;}

  //! Returns the number of flops in the application of the preconditioner.
  double ApplyInverseFlops() const {return 0;}

  //! Prints basic information on iostream. This function is used by operator<<.
  std::ostream& Print(std::ostream& os) const {os << "(Identity Operator)"; return os;}

  //@}

 protected:
//--------//

	  //! communicator
	  const Epetra_Comm& comm;

	  //! range map (map for vectors y in Ax=y)
	  const Epetra_Map& rangeMap;

	  //! domain map (map for vectors x in Ax=y)
	  const Epetra_Map& domainMap;

	  bool useTranspose;

  }; // end of class IdentityOperator


//! this operator implements P\A-I, which is used for Eigen-Analysis of
//! preconditioned matrices in the solver factory:
class PAmI_Operator : public Epetra_Operator
  {

 public:
//-----//


  //! constructor
  PAmI_Operator(Teuchos::RCP<const Epetra_Operator> _A, Teuchos::RCP<const Epetra_Operator> _P)
    : A(_A), P(_P)
    {}

  //! destructor
  virtual ~PAmI_Operator(){}

//!\name  Atribute set methods
//!@{

  //! If set true, transpose of this operator will be applied (n/a)
  int 	SetUseTranspose (bool UseTranspose) {return -1;}

//!@}

//!\name Mathematical functions
//!@{

  //! apply operator Y=Op*X
  int 	Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
    Y=X;
    Epetra_MultiVector tmp=Y;
    EPETRA_CHK_ERR(A->Apply(X,tmp));
    EPETRA_CHK_ERR(P->ApplyInverse(tmp,Y));
    EPETRA_CHK_ERR(Y.Update(-1.0,X,1.0));
    return 0;
    }

  //! apply inverse operator (n/a)
  int 	ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
    return -1;
    }

  //! Returns the infinity norm of the global matrix.
  double NormInf () const {return -1.0;}


//!@}


//!\name Atribute access functions
//!@{

  //!Returns a character string describing the operator.
  const char * 	Label () const {return "P\\A-I";}

  //! Returns the current UseTranspose setting.
  bool 	UseTranspose () const {return false;}

  //! Returns true if the this object can provide an approximate Inf-norm, false otherwise.
  bool 	HasNormInf () const {return false;}

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & 	Comm () const {return A->Comm();}

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & 	OperatorDomainMap () const
    {
    return A->OperatorDomainMap();
    }

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & 	OperatorRangeMap () const
    {
    return A->OperatorRangeMap();
    }

  //@}

 protected:
//--------//

  Teuchos::RCP<const Epetra_Operator> A;
  Teuchos::RCP<const Epetra_Operator> P;

  }; // end of class PAmI_Operator

}//namespace TRIOS

#endif // SOLVERFACTORY_H
