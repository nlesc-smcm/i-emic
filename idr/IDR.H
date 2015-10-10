#ifndef IDRSolver_H
#define IDRSolver_H

#include <vector>
#include <memory>

// Templated types are assumed to be shared_pointers: we use -> in calls to
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

	double tol_;
	bool   smoothing;

public:
	// constructor
	IDRSolver(Model model, Vector x0);
	
	// destructor
	~IDRSolver();
	
	Vector solve(Vector b);	
	
};

//====================================================================
//====================================================================
//====================================================================
template<typename Model, typename Vector>
IDRSolver<Model, Vector>::
IDRSolver(Model model, Vector x0)
	:
	model_(model),  
	x_(x0),
{
	tol_ = 1e-6;
	smoothing_ = false;
}

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

	double normb = b->norm();
	double tolb  = tol_ * normb; // Relative tolerance

	// Compute residual
	Vector Ax = model_->applyMatrix(*x_);     // Perform matvec
	Vector r  = std::make_shared<Vector>(*b); // Initialize r with b
	r->update( 1.0,  *b, 0.0);                // Replace contents with rhs
	r->update(-1.0, *Ax, 1.0);                // r = b - A*x

	Vector x_s, r_s;
	if (smoothing_)
	{
		x_s = std::make_shared<Vector>(*x_);
		r_s = std::make_shared<Vector>(*r);
	}

	double normr = r->norm();
	std::vector<double> resvec();
	resvec.push_back(norm);
	double trueres = 0.0;

	// Initial guess is a good enough solution
	// TODO

	int iter = 0;
	while (normr > tolb && iter < maxit)
	{
		// kijk in fortran code
	}
	
}

#endif
