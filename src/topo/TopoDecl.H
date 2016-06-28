//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef TOPODECL_H
#define TOPODECL_H

#include <vector>

#include "GMRESSolver.H"

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
	using Matrix    = typename VectorPtr::element_type;
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

	//! State
	//!  is initialized with a view from the model's state
	VectorPtr stateView_;

	//! Right hand side
	//!  is initialized with a copy from the model's rhs	
	VectorPtr rhs_;

	//! Solution from solve()
	//!  is initialized with a view from the model's solution	
	VectorPtr solView_;

	//! Integer for selecting the right landmask.
	//!  k_ = floor(delta_)
	int k_;

 	//! Vectors of integers for selecting the right landmask.
	//!  a_[i+1] = a_[i] + 2 * ((k_+1) % 2); a_[0] = 0
	//!  b_[i+1] = b_[i] + 2 * ( k_ % 2 );   b_[0] = 1
	std::vector<int> a_, b_;

	//! Pointers to Jacobians
	MatrixPtr matA_, matB_;

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

private:
  	//! load the mask filenames
	void loadMaskFileNames();

	//! load the masks
	void loadMasks();

};

#endif
