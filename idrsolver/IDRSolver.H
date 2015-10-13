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
	//-->TODO

	// Number close to machine precision:
	double mp = 1e-13;

	// Initialize output paramater relres
	double relres = nan("");

	double normb = b->norm();
	double tolb  = tol_ * normb; // Relative tolerance

	// Compute residual
	VectorPointer Ax = model_->applyMatrix(*x_);     // Perform matvec
	Vector r(*b);                                    // Initialize r with b
	r.update( 1.0,  *b, 0.0);                        // Replace contents with rhs
	r.update(-1.0, *Ax, 1.0);                        // r = b - A*x

	// Constructing smoothing vectors xs_ and rs_
	if (smoothing_)
	{
		xs_ = Vector(*x_);
		rs_ = Vector(r);
	}
	
	double normr = r.norm();
	std::vector<double> resvec;
	resvec.push_back(normr);
	double trueres = 0.0;

	// Initial guess is a good enough solution
	//-->TODO

	//-->RECYCLING!!
	bool inispace = false; // true if an initial search space is available
	
	int   iter = 0; 
	int     ii = 0;   // inner iteration counter
	int     jj = 0;   // G-space counter
	
	size_t dim = s_;
	double om  = 1.0;
	
	std::vector<double> f;
	std::vector<double> gamma(dim, 0.0);
	
	std::vector<VectorPointer>        G(dim, VectorPointer());
	std::vector<VectorPointer>        U(dim, VectorPointer());
	std::vector<std::vector<double> > M(dim, std::vector<double>(dim, 0.0));
	
	while (normr > tolb && iter < maxit_)
	{
		// Create new right hand side for small system:
		for (int i = 0; i < s_; ++i)
			f.push_back(r.dot(P_[i]));

		for (int k = 0; k < s_; ++k)
		{
			// Update inner iteration counter
			ii = ii + 1;
				
			// Create new vectors v,t
			Vector v(r);

			if (jj > 0)
			{				
				// Solve small lower triangular system and make v orthogonal to P:
				for (int i = k; i < s_; ++i)
				{
					gamma[i] = f[i];
					for (int j = k; j < i-1; ++j)
					{
						gamma[i] = gamma[i] - M[i][j] * gamma[j];
					}
					gamma[i] = gamma[i] / M[i][i];
					v.update(-gamma[i], *(G[i]), 1.0);
				}

				// Preconditioning
				//-->needs checks, now we are assuming we do right prec
				VectorPointer t = model_->applyPrecon(v);
				t->scale(om);

				// Compute new U(:,k)
				for (int i = k; i < s_; ++i)
					t->update(gamma[i], *(G[i]), 1.0);
				U[k] = t;

				// Compute Hessenberg matrix
				//-->TODO
			}
			else if (!inispace) //only when an initial search space unavailable
			{
				U[k] = model_->applyPrecon(v);
			}

			// Compute new G(:,k), G(:,k) is in space G_j
			G[k] = model_->applyMatrix(*(U[k]));

			// Bi-Orthogonalise the new basis vectors:
			//-->TODO
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
