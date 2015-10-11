#ifndef IDRSolver_H
#define IDRSolver_H

#include <vector>
#include <memory>

#include "IDRSolverDecl.H"

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
	tol_       = 1e-6;
	smoothing_ = false;
	s_         = 4;
	maxit_     = 500;

	createP();
}

//====================================================================
template<typename ModelPointer, typename VectorPointer>
IDRSolver<ModelPointer, VectorPointer>::
~IDRSolver()
{}

//====================================================================
template<typename ModelPointer, typename VectorPointer>
void IDRSolver<ModelPointer, VectorPointer>::
createP()
{
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
	
	int   iter = 0;
	size_t dim = s_;
	std::vector<double> f;
	std::vector<std::vector<double> > M(dim, std::vector<double>(dim, 0.0));
	while (normr > tolb && iter < maxit_)
	{
		// Create new right hand side for small system:
		for (int i = 0; i < s_; ++i)
			f.push_back(r->dot(P_[i]));

		for (int k = 0; k < s_; ++k)
		{
			// Solve small lower triangular system and make v orthogonal to P:
			double gamma;
			for (int i = k; i < s_; ++i)
			{
				gamma = f[i];
				for (int j = k; j < i-1; ++j)
				{
					
				}				
			}
		}
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
