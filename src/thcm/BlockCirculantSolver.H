#ifndef BLOCK_CIRCULANT_SOLVER_H
#define BLOCK_CIRCULANT_SOLVER_H

#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_Scaling.H"
#ifdef HAVE_MPI
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#ifdef HAVE_FFTW
#include "fftw.h"
#endif
#include "EpetraExt_RowMatrixOut.h"

class Epetra_LinearProblem;
class Komplex_LinearProblem;
class Amesos_BaseSolver;
class Ifpack_Preconditioner;

class BlockCirculantSolver : public Epetra_Operator,
							 public NOX::Epetra::Interface::Preconditioner
{
      
public:
      
	/*!                                                               
	  \brief Constructor                                              
	*/
	BlockCirculantSolver(const EpetraExt::BlockCrsMatrix& jacobian,
						 const EpetraExt::BlockVector& cloneVector,
						 const Teuchos::RCP<const Epetra_CrsMatrix> A,
						 const Teuchos::RCP<const Epetra_CrsMatrix> B,
						 const Teuchos::RCP<const EpetraExt::MultiComm> globalComm_,
						 double dt_, Teuchos::ParameterList& params);
      
	//! Destructor
	virtual ~BlockCirculantSolver();
      
	//! Set transpose.
	virtual int SetUseTranspose(bool UseTranspose);
      
	//! Apply XYZT preconditioner operator
	virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
      
	//! Apply XYZT preconditioner operator inverse
	virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

	//! setup preconditioner
	bool computePreconditioner(const Epetra_Vector& x,
							   Epetra_Operator& Prec, Teuchos::ParameterList* p = NULL);


      
	//! Computing infinity norm
	virtual double NormInf() const;
      
	//! Label
	virtual const char* Label() const;
      
	//! Transpose
	virtual bool UseTranspose() const;
      
	//! Have norm-inf
	virtual bool HasNormInf() const;
      
	/*! 
	 * \brief Returns a pointer to the Epetra_Comm communicator associated 
	 * with this operator.
	 */
	virtual const Epetra_Comm& Comm() const;
      
	/*! 
	 * \brief Returns the Epetra_Map object associated with the domain of 
	 * this operator.
	 */
	virtual const Epetra_Map& OperatorDomainMap() const;
      
	/*! 
	 * \brief Returns the Epetra_Map object associated with the range of 
	 * this operator.
	 */
	virtual const Epetra_Map& OperatorRangeMap() const;
            
	//! Exception handler for the XYZT preconditioner class.
	virtual void throwError(const std::string& functionName, 
							const std::string& errorMsg) const;

	//! creates a communicator in the t-direction
	Teuchos::RCP<const Epetra_Comm> getTimeComm() const;
            
	//! set scaling object for diagonal blocks
	void resetScaling(Teuchos::RCP<const NOX::Epetra::Scaling> scal)
        {
			scaling=scal;
        }
            
protected:
	//!globalMap standardMap
	Teuchos::RCP<const Epetra_Map> standardMap;
	Teuchos::RCP<const Epetra_BlockMap> globalMap;

	//! parameter list obtained from user
	Teuchos::ParameterList& params;
    
	//! Reference to global communicator for global XYZT operator.
	const Teuchos::RCP<const EpetraExt::MultiComm> globalComm;

	//! reference to actual jacobian (just used for cloning)
	const EpetraExt::BlockCrsMatrix& jacobian;
      

	//! communicator for only the spatial decomposition,
	//! i.e. all procs in this comm own the same time steps
	//! but on different parts of the computational domain.
	Teuchos::RCP<Epetra_MpiComm> localComm;

	//! communicator in the time-direction. All procs in this comm
	//! have the same part of the spatial domain but different time steps.
	Teuchos::RCP<const Epetra_Comm> timeComm;
      
	//! Pointer to time-averaged jacobian block A+1/dtB
	const Teuchos::RCP<const Epetra_CrsMatrix> jacobianBlock;
      
	//! Pointer to mass matrix block -1/dt B
	const Teuchos::RCP<const Epetra_CrsMatrix> massBlock;
      
	//! scaling object for diagonal blocks
	Teuchos::RCP<const NOX::Epetra::Scaling> scaling;
      
	//! time step-size
	double dt;
      
	//! Label for this operator.
	std::string label;
      

	//! number of blocks, blocksize, blocks per process
	int numBlocks;
	int blocksOnDomain;
	int blockSize;

	Teuchos::Array<Teuchos::RCP<Epetra_Vector> > rhsRe,rhsIm,lhsRe,lhsIm;
	Teuchos::Array<Teuchos::RCP<Amesos_BaseSolver> > Solver;

	//! real representation of complex linear problem
	Teuchos::Array<Teuchos::RCP<Komplex_LinearProblem> > complexProblem;

#ifdef HAVE_FFTW
	//! fft plan and temp vectors
	fftw_complex *data;
	fftw_complex *work;
	fftw_plan planf;
	fftw_plan planb;
#endif      
	Teuchos::RCP<const Epetra_MultiVector> nullSpace_;

	Teuchos::RCP<Epetra_MultiVector> FFTNullSpace;
	//Teuchos::RCP<Epetra_MultiVector> nullSpace_;

private:
    
	//! private setup function called by constructor
	void buildTimeComm(void);  
      
	//! postprocess solution from a single time step
	//! projects out the null-space component from a
	//! computed solution
	int PostProcSolution(Epetra_Vector& soln) const;

	//!get NullSpace 
	int getNullSpace();

	//! replace one row and column by a Dirichlet condition (0 everywhere, 1 on diagonal).
	//! This function assumes that the pattern of the matrix is symmetric!
	int PutDirichlet(Epetra_CrsMatrix& A, int gid);

                    
};

#endif // BLOCK_CIRCULANT_SOLVER_H
