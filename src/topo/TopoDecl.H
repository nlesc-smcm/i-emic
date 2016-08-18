//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef TOPODECL_H
#define TOPODECL_H

#include <vector>

#include <Epetra_MultiVector.h>

#include <Teuchos_RCP.hpp>
#include "BelosConfigDefs.hpp"
#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"
#include "BelosTypes.hpp"
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>

//! Forward declarations
template<typename Model>
struct Combined_Operator;

namespace Belos
{
	template <typename Model>
	class OperatorTraits <double, Epetra_MultiVector,
						  Combined_Operator<Model> >;
}

//! This class manages a Model that allows changing its topography. A
//! continuation with Topo<Model> model is then possible through the
//! continuation parameter delta_.
//!
//! We expect a Model, which should be a (reference counted) pointer
//! to an I-EMIC model. A Model should maintain its own Vector,
//! which we expect to be of pointer type as well and named VectorPtr.
//!
//! In addition to the requirements from the Continuation class, we
//! need a Model to implement member functions:
//!
//!   void          setLandMask()
//!   LandMask      getLandMask()
//! 
//! The templated type ParameterList should be a pointer 
//! to a parameter storing object having a .get() method.

template<typename Model, typename ParameterList>
class Topo
{

public:
	
	using VectorPtr = typename Model::element_type::VectorPtr;
	using Vector    = typename VectorPtr::element_type;
	using MatrixPtr = typename Model::element_type::MatrixPtr;
	using Matrix    = typename MatrixPtr::element_type;
	using PreconPtr = typename Model::element_type::PreconPtr;
	using Precon    = typename PreconPtr::element_type;
	using LandMask  = typename Model::element_type::LandMask;

private: 
	
	//! Model
	Model model_;

	//! Parameterlist
	ParameterList pars_;

	//! Number of masks we need to load
	int nMasks_;

	//! Vector containing all the filenames
	std::vector<std::string> landMaskFileNames_;

	//! Vector containing RCP's to the landmask arrays
	std::vector<LandMask> landMasks_;

	//! Continuation parameter
	//!  ranging continually from 0 to nMasks-1
	double delta_;

	//! The model's continuation parameter, which we might
	//! want to affect as well
	double modelPar_;

	//! State
	VectorPtr stateView_;
	VectorPtr stateCopy_;

	//! We keep our own stateDot
	VectorPtr stateDot_;

	//! Right hand side
	//!  is initialized with a copy from the model's rhs	
	VectorPtr rhs_;

	//! Derivative with respect to delta
	VectorPtr dFdPar_;

	//! Solution from solve()
	//!  is initialized with a view from the model's solution	
	VectorPtr solView_;

	//! Integer for selecting the right landmask.
	//!  k_ = floor(delta_)
	int k_;

	//! Homotopy factors
	double facA_, facB_;

 	//! Vectors of integers for selecting the right landmask.
	//!  a_[i+1] = a_[i] + 2 * ((k_+1) % 2); a_[0] = 0
	//!  b_[i+1] = b_[i] + 2 * ( k_ % 2 );   b_[0] = 1
	std::vector<int> a_, b_;

	//! Pointers to Jacobians
	MatrixPtr matA_, matB_;

	//! Pointers to preconditioner
	PreconPtr precA_, precB_;

	//! Combined operators
	Combined_Operator<Model> combMat_;
	Combined_Operator<Model> combPrec_;

	//! If every converged state deserves to be stored, enable this
	bool storeEverything_;

	//! Solver parameterlist
	ParameterList solverParams_, belosParamList_;

	Teuchos::RCP
	<Belos::LinearProblem
	 <double, Epetra_MultiVector, Combined_Operator<Model> > > problem_;	
	
	Teuchos::RCP
	<Belos::BlockGmresSolMgr
	<double, Epetra_MultiVector, Combined_Operator<Model> > > belosSolver_;

	//! Initialization flag
	bool solverInitialized_;

	//! Array keeping track of initialized preconditioners
	std::vector<bool> initPrecs_;

	//! Relaxation factors for predictors
	double rA_, rB_;

	//! Enable this to use predictors
	bool usePredictors_;

	//! Enable this to dynamically estimate relaxation factors
	bool calcRelax_;	

public:	
	//! default constructor
	Topo(){};
	
	//! constructor
	Topo(Model model, ParameterList pars);
	
	//! destructor
	~Topo();

	//! get continuation parameter
	double getPar();

	//! set continuation parameter
	//!  this also sets k_, facA_ and facB_
	void setPar(double par);

	//! get the state of the model
	//!   mode: 'C' copy
	//!         'V' view
	VectorPtr getState(char mode = 'C');

	//! get the right hand side of the model
	//!   mode: 'C' copy
	//!         'V' view
	VectorPtr getRHS(char mode = 'C');

	//! get the solution of the model after a call to solve()
	//!   mode: 'C' copy
	//!         'V' view
	VectorPtr getSolution(char mode = 'C');
	
	//! compute right hand side
	void computeRHS();

	//! compute Jacobian matrix
	void computeJacobian();

	//! compute derivative of RHS with respect to delta
	void computeDFDPar();

	//! build Preconditioner
	void buildPreconditioner();

	//! solve Jx=b
	void solve(VectorPtr b);

	//! apply Jacobian matrix J*v
	void applyMatrix(Vector const &v, Vector &out);

	//! apply Preconditioning inv(P)*v
	void applyPrecon(Vector const &v, Vector &out);

	//! perform some duties before important things happen
	void preProcess();

	//! perform some duties after important things have happened
	void postProcess();

	//! get a_ for testing
	std::vector<int> getA() { return a_; }

	//! get b_ for testing
	std::vector<int> getB() { return b_; }

	//! return number of masks
	int nMasks() { return nMasks_; }

	//! test different RHS
	void testRHS();
private:
  	//! load the mask filenames
	void loadMaskFileNames();

	//! load the masks
	void loadMasks();

	//! compute residual norm ||b - facA*A*x - facB*B*x||
	double computeResidual(VectorPtr b);

	//! initialize solver
	void initializeSolver();

};

//------------------------------------------------------------------
template <typename Model>
struct Combined_Operator
{
	using MatrixPtr = typename Model::element_type::MatrixPtr;
	using PreconPtr = typename Model::element_type::PreconPtr;
	MatrixPtr A;
	MatrixPtr B;
	PreconPtr PA;
	PreconPtr PB;
	double facA;
	double facB;
	bool isPrec;
	char precType;
};

//------------------------------------------------------------------
// Specialization of OperatorTraits
//------------------------------------------------------------------
namespace Belos
{
	template<typename Model>
	class OperatorTraits <double, Epetra_MultiVector,
						  Combined_Operator<Model> >
	{
	public:
		static void
		Apply (const Combined_Operator<Model> &Op,
			   const Epetra_MultiVector &x,
			   Epetra_MultiVector &y,
			   int trans = 0)
			{
				if (Op.isPrec)
				{
					if (Op.precType == 'H') // Use homotopy
					{
						Epetra_MultiVector tmp = y;	
						Op.PA->ApplyInverse(x, tmp);
						Op.PB->ApplyInverse(x, y);
						y.Update(Op.facA, tmp, Op.facB);
					}
					else if (Op.precType == 'M') 
					{
						if (Op.facB < 0.95)
							y.Update(1.0, x, 0.0);
						else
							Op.PB->ApplyInverse(x, y);
					}
					else if (Op.precType == 'N') // Apply scaling
					{
						if (Op.facA >= Op.facB)
						{
							Op.PA->ApplyInverse(x, y);
							y.Scale(1.0 / Op.facA);
						}
						else
						{
							Op.PB->ApplyInverse(x, y);
							y.Scale(1.0 / Op.facB);
						}
					}
				}
				else
				{
					if (Op.facB < 0.95)
					{
						Op.B->Apply(x, y);
						y.Update(Op.facA, x, Op.facB);
					}
					else
						Op.B->Apply(x, y);
				}
			}
		
		static bool
		HasApplyTranspose (const Combined_Operator<Model> &Op) {return false;}
	};
}

#endif
