#ifndef IDRSolver_H
#define IDRSolver_H

#include <vector>

// Templated types are assumed to be pointers: we use -> in calls to
// their members.

// Model should be a class with members:
//    -applyMatrix(Vector v), performing matrix vector product t=Av
//    -applyPrecon(Vector v), applying the operation v = P^{-1} x
// Note that Model should be compatible with Vector

// Vector should be a class with members:
//    -update(double scalarA, Vector A, double scalarThis, performing
//      this = scalarA * A + scalarThis * this
//    -norm()
//    -copy construction


template<typename Model, typename Vector>
class IDRSolver
{
	Model  model_;
	Vector x_;

public:
	// constructor
	IDRSolver(Model A, Vector x, Vector r);

	// destructor
	~IDRSolver();
	
	Vector solve(Vector b);	
	
};

//====================================================================
//====================================================================
//====================================================================
template<typename Model, typename Vector>
IDRSolver<Model, Vector>::
IDRSolver(Model A, Vector x, Vector r)
	:
	A_(A),
	x_(x),
	model_(M)
{}

// ===================================================================
template<typename Model, typename Vector>
Vector IDRSolver<Model, Vector>::
solve(Vector b)
{
	// check for zero rhs
	// TODO

	// Number close to machine precision:
	double mp = 1e-13;

	// Initialize output paramater relres
	double relres = nan();

	Vector 
	
}

#endif
