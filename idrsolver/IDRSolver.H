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


template<typename ModelPointer, typename VectorPointer>
class IDRSolver
{
	// We require the pointers to be of shared_ptr/RCP type
	using  Model  = typename ModelPointer::element_type;	
	using  Vector = typename VectorPointer::element_type;
	
	ModelPointer  model_;
	VectorPointer x_;

	double tol_;
	bool   smoothing_;
	int    s_;
	int    maxit_;
	// smoothing vectors
	Vector xs_;
	Vector rs_;

	// shadow space?
	std::vector<Vector> P_;

public:
	// constructor
	IDRSolver(ModelPointer model, VectorPointer x0);
	
	// destructor
	~IDRSolver();
	
	VectorPointer solve(VectorPointer b);

	void test();
	
};

//====================================================================
//====================================================================
//====================================================================
template<typename ModelPointer, typename VectorPointer>
IDRSolver<ModelPointer, VectorPointer>::
IDRSolver(ModelPointer model, VectorPointer x0)
	:
	model_(model),  
	x_(x0)
{
	tol_ = 1e-6;
	smoothing_ = false;
	s_ = 4;
	maxit_ = 500;
	double alpha;
	for (int j = 0; j < s_; ++j)
	{
		Vector p(*x_);
		p.random();
		P_.push_back(p);
		for (int k = 0; k < j; ++k)
		{
			alpha = (P_[k]).dot(P_[j]);
			P_[j].update(-alpha, P_[k], 1.0);
		}
		P_[j].scale(1.0/P_[j].norm());
	}

}

//====================================================================
template<typename ModelPointer, typename VectorPointer>
IDRSolver<ModelPointer, VectorPointer>::
~IDRSolver()
{}

// ===================================================================
template<typename ModelPointer, typename VectorPointer>
VectorPointer IDRSolver<ModelPointer, VectorPointer>::
solve(VectorPointer b)
{
	// check for zero rhs
	// TODO

	// Number close to machine precision:
	double mp = 1e-13;

	// Initialize output paramater relres
	double relres = nan("");

	double normb = b->norm();
	double tolb  = tol_ * normb; // Relative tolerance

	// Compute residual
	VectorPointer Ax = model_->applyMatrix(*x_);     // Perform matvec
	VectorPointer r  = std::make_shared<VectorPointer>(*b); // Initialize r with b
	r->update( 1.0,  *b, 0.0);                // Replace contents with rhs
	r->update(-1.0, *Ax, 1.0);                // r = b - A*x

	// Constructing smoothing vectors xs_ and rs_
	if (smoothing_)
	{
		xs_(*x_);
		rs_(*r);
	}
	
	double normr = r->norm();
	std::vector<double> resvec();
	resvec.push_back(normr);
	double trueres = 0.0;

	// Initial guess is a good enough solution
	// TODO

	int iter = 0;
	std::vector<double> f;
	while (normr > tolb && iter < maxit_)
	{
		// Create new right hand side for small system:
		for (int i = 0; i < s_; ++i)
			f.push_back(r->dot(P_[i]));		
	}
	
}

//====================================================================
template<typename ModelPointer, typename VectorPointer>
void IDRSolver<ModelPointer, VectorPointer>::
test()
{
	std::cout << "========== Testing IDR solver ===============" << std::endl;
	int psize = P_.size();
	std::cout << "size of P: " << P_.size() << std::endl;
	for (int i = 0; i < psize; ++i)
	{
		for (int j = 0; j < psize; ++j)
		{
			std::cout << "i=" << i << " j=" << j << " P(:,i)^T P(:,j)="
					  << P_[i].dot(P_[j]) << std::endl;
		}
	}
	std::cout << "========== Testing IDR solver finished  =====" << std::endl;
}

#endif
