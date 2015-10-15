#ifndef IDRSolver_H
#define IDRSolver_H

#include <vector>
#include <memory>

#include "IDRSolverDecl.H"

//====================================================================
template<typename ModelPointer, typename VectorPointer>
IDRSolver<ModelPointer, VectorPointer>::
IDRSolver(ModelPointer model,
		  VectorPointer x0,
		  VectorPointer b)
	:
	model_(model),  // model
	x_(x0),         // initial guess
	b_(b)           // RHS
{
	tol_         = 1e-6;
	angle_       = 0.7;
	smoothing_   = false;
	s_           = 4;
	maxit_       = 500;
	mp_          = 1e-13; // number close to machine precision
	replacement_ = false;
	trueres_     = false;

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
int IDRSolver<ModelPointer, VectorPointer>::
solve()
{
	// check for zero rhs
	std::cout << "IDR: norm rhs = " << b_->norm() << std::endl;
	//-->TODO

	// Number close to machine precision:
	double mp = 1e-13;

	// Initialize output paramater relres
	double relres = nan("");

	double normb = b_->norm();
	double tolb  = tol_ * normb; // Relative tolerance

	// Compute residual r = b - Ax
	VectorPointer r = model_->applyMatrix(*x_); // Ax
	std::cout << "IDR: Ax = " << r->norm() << std::endl;
	r->update( 1.0,  *b_, -1.0);                 // b - Ax
	
	// Constructing smoothing vectors xs_ and rs_
	if (smoothing_)
	{
		xs_ = Vector(*x_);
		rs_ = Vector(*r);
	}
	
	double normr = r->norm();
	std::cout << "IDR: initial residual = " << normr << std::endl;
	resvec_.push_back(normr);
	double trueres = 0.0;

	// Initial guess is a good enough solution
	//-->TODO

	//-->RECYCLING!!
	bool inispace = false; // true with initial search space

	int   flag = 0;   // a flag!
	int   iter = 0; 
	int     ii = 0;   // inner iteration counter
	int     jj = 0;   // G-space counter

	int replacements_ = 0;
	
	size_t dim = s_;
	double om  = 1.0;
	
	std::vector<double> f;
	std::vector<double> gamma(dim, 0.0);
	
	std::vector<VectorPointer>        G(dim, VectorPointer());
	std::vector<VectorPointer>        U(dim, VectorPointer());
	std::vector<std::vector<double> > M
		(dim, std::vector<double>(dim, 0.0));
	// Declare vector v
	Vector v;
	
	while (normr > tolb && iter < maxit_)
	{
		std::cout << normr << " " << tolb << " " << iter << std::endl;
		// Create new right hand side for small system:
		for (int i = 0; i < s_; ++i)
			f.push_back(r->dot(P_[i]));

		for (int k = 0; k < s_; ++k)
		{
			// Update inner iteration counter
			ii = ii + 1;
				
			// Compute new v (copy construction)
			v = Vector(*r);

			if (jj > 0)
			{				
				// Solve small lower triangular system and make
				// v orthogonal to P:
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
				//-->needs checks, now we are assuming
				//   we do right preconditioning
				VectorPointer t = model_->applyPrecon(v);
				t->scale(om);

				// Compute new U(:,k)
				for (int i = k; i < s_; ++i)
					t->update(gamma[i], *(G[i]), 1.0);
				U[k] = t;

				// Compute Hessenberg matrix
				//-->TODO
			}
			else if (!inispace) 
			{
				U[k] = model_->applyPrecon(v);
			}

			// Compute new G(:,k), G(:,k) is in space G_j
			G[k] = model_->applyMatrix(*(U[k]));

			// Bi-Orthogonalise the new basis vectors:
			double alpha;
			for (int i = 0; i < k-1; ++i)
			{
				alpha = P_[i].dot(*(G[k])) / M[i][i];
				G[k]->update(-alpha, *(G[i]), 1.0);
				U[k]->update(-alpha, *(U[i]), 1.0);
			}

			// New column of M = P'*G  (first k-1 entries are zero)
			for (int i = k; i < s_; ++i)
			{
				M[i][k] = G[k]->dot(P_[i]);
			}
			if (M[k][k] == 0)
			{
				flag = 3;
				std::cout << "WARNING: M[k][k] == 0"  << std::endl;
				return -1;
			}

			// Make r orthogonal to p_i, i = 1..k, update solution and residual
			double beta = f[k] / M[k][k];
			r->update( -beta, *(G[k]), 1.0); // r = r - beta*G(:,k);
			x_->update(beta, *(U[k]), 1.0);  // x = x + beta*U(:,k);

			// Check whether we need to replace residual
			normr = r->norm();			
			if (replacement_ && normr > tolb / mp_) trueres_ = true;

			// Smoothing
			//-->TODO

			// New f = P'*r (first k  components are zero)
			//-->Here we use s-1 since we are zero-based!
			if (k < s_-1)
			{
				for (int i = k+1; i < s_; ++i)
				{
					f[i] = f[i] - beta*M[i][k];
				}
			}

			// Check for convergence
			normr = r->norm();
			resvec_.push_back(normr);
			iter = iter + 1;
			if (normr < tolb || iter == maxit_)
				break;
		}

		// Check for convergence
		if (normr < tolb || iter == maxit_)
			break;

		// Update G-space counter
		jj = jj + 1;

		// Compute first residual in G_j+1
		// Note: r is already perpendicular to P so v = r
		// Preconditioning: v = inv(P)*r
		v = *(model_->applyPrecon(*r));

		// Matrix-vector multiplication: t = A*v
		VectorPointer t = model_->applyMatrix(v);

		// Computation of a new omega
		om = calc_omega(*t, *r, angle_);
		
		// Update solution and residual:
		r->update(-om, *t, 1.0);  // r = r - om*t
		x_->update(om, v, 1.0);   // x = x + om*v 

		normr = r->norm();
		
		// Residual replacement?
		if (replacement_ && normr > tolb / mp_) trueres_ = 1;
		if (trueres_ && normr < normb )
		{
			// r = b - A*x;
			r = model_->applyMatrix(*x_); // Ax
			r->update( 1.0,  *b_, -1.0);   // b - Ax
			trueres_ = false;
			replacements_ = replacements_ + 1;
		}

		// Smoothing
		//-->TODO

		resvec_.push_back(normr);
		iter = iter + 1;		
	}

	// Smoothing
	//-->TODO

	
	
}

//====================================================================
template<typename ModelPointer, typename VectorPointer>
void IDRSolver<ModelPointer, VectorPointer>::
printResVec()
{
	int ctr = 0;
	for (auto &i: resvec_)
		std::cout << ctr << " " << i << std::endl;
}

//====================================================================
template<typename ModelPointer, typename VectorPointer>
double IDRSolver<ModelPointer, VectorPointer>::
calc_omega(Vector const &t, Vector const &s, double angle)
{
	double ns = s.norm();
	double nt = t.norm();
	double ts = t.dot(s);

	double rho = std::abs(ts / (nt * ns));

	double om = ts / (nt * nt);

	if (rho < angle)
		om = om * angle / rho;

	return om;
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
