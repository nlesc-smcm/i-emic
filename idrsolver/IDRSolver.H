#ifndef IDRSolver_H
#define IDRSolver_H

#include <vector>
#include <memory>

#include "IDRSolverDecl.H"

//====================================================================
template<typename Model, typename VectorPointer>
IDRSolver<Model, VectorPointer>::
IDRSolver(Model &model)
	:
	model_         (model),  // model
	haveInitSol_   (false),
	haveRHS_       (false),
	smoothing_     (false),
	replacement_   (false),
	trueres_       (false),
	verbosity_     (0),
	inispace_      (false),
	mp_            (1e-13),  // number close to machine precision
	tol_           (1e-8),
	angle_         (0.7),
	s_             (4),
	maxit_         (1000),
	iter_          (0)
{}

//====================================================================
template<typename Model, typename VectorPointer>
IDRSolver<Model, VectorPointer>::
IDRSolver(Model &model,
		  VectorPointer x0,
		  VectorPointer b)
	:
	model_         (model),  // model
	x_             (x0),     // initial guess
	b_             (b),      // RHS
	haveInitSol_   (true),
	haveRHS_       (true),
	smoothing_     (false),
	replacement_   (false),
	trueres_       (false),
	verbosity_     (0),
	inispace_      (false),
	mp_            (1e-13),  // number close to machine precision
	tol_           (1e-6),
	tolb_          (1e-6),
	angle_         (0.7),
	s_             (4),
	maxit_         (1000),
	iter_          (0)
{}

//====================================================================
template<typename Model, typename VectorPointer>
IDRSolver<Model, VectorPointer>::
~IDRSolver() {}

//====================================================================
template<typename Model, typename VectorPointer>
template<typename ParListPtr>
void IDRSolver<Model, VectorPointer>::
setParameters(ParListPtr pars)
{
	if (verbosity_ > 5)
		std::cout << "IDR: updating parameters" << std::endl;
	
	s_           = pars->get("IDR s", 4);
	angle_       = pars->get("IDR angle", 0.7);
	tol_         = pars->get("IDR tolerance", 1e-8);
	maxit_       = pars->get("IDR iterations", 500);
	inispace_    = pars->get("IDR save search space", true);
	smoothing_   = pars->get("IDR use smoothing", false);
	replacement_ = pars->get("IDR replace residuals", false);
	verbosity_   = pars->get("IDR verbosity", 0);

}

//====================================================================
template<typename Model, typename VectorPointer>
void IDRSolver<Model, VectorPointer>::
createP()
{
	P_.clear();
	P_ = std::vector<Vector>(s_, Vector());
	double alpha;
	for (int j = 0; j < s_; ++j)
	{
		Vector p(*x_);
		p.random();
		P_[j] = p;
		for (int k = 0; k < j; ++k)
		{
			alpha = (P_[k]).dot(P_[j]);
			P_[j].update(-alpha, P_[k], 1.0);
		}
		P_[j].scale(1.0/P_[j].norm());

		if (verbosity_ > 9) P_[j].print();
	}
	
}

// ===================================================================
template<typename Model, typename VectorPointer>
int IDRSolver<Model, VectorPointer>::
solve()
{
	if (P_.empty())
		createP();

	if (!haveInitSol_ || !haveRHS_)
	{
		std::cout << "Problem not setup correctly!"
				  << haveInitSol_ << " "
				  << haveRHS_ << std::endl;
		return 1;
	}			
	
	// check for zero rhs
	//-->TODO
	
	double normb = b_->norm();
	tolb_  = tol_ * normb; // Relative tolerance

	// Compute residual r = b - Ax
	Vector r(*x_);
	model_.applyMatrix(*x_, r);
	r.update( 1.0,  *b_, -1.0);               // b - Ax
	
	// Constructing smoothing vectors xs_ and rs_
	if (smoothing_)
	{
		xs_ = Vector(*x_);
		rs_ = Vector(r);
	}
	
	normr_ = r.norm();
	
	resvec_.clear();
	resvec_.push_back(normr_);
	
	// Initial guess is a good enough solution
	//-->TODO

	int   flag = 0;   // a flag!
	int     ii = 0;   // inner iteration counter
	int     jj = 0;   // G-space counter

	int replacements_ = 0;
	
	size_t dim = s_;
	double om  = 1.0;
	
	std::vector<double> f(dim, 0.0);
	std::vector<double> gamma(dim, 0.0);
	
	std::vector<Vector> G(dim, Vector(*x_));
	
	// Determine whether the search space needs initialization
	bool createInitU =
		(U_init_.size() == 0) || (!inispace_)  ||
		U_init_.empty()      || (iter_ < s_);

	std::vector<Vector> U;
	if (createInitU)
	{
		U       = std::vector<Vector>(dim, Vector());
		U_init_ = std::vector<Vector>(dim, Vector(*x_));
	}	
	else
		U = U_init_;
	
	std::vector<std::vector<double> > M
		(dim, std::vector<double>(dim, 0.0));
	
	Vector v; 
	Vector t(r);
	
	iter_ = 0;
	if (verbosity_ > 4) printIterStatus();
	while (normr_ > tolb_ && iter_ < maxit_)
	{
		
		// Create new right hand side for small system:
		for (int i = 0; i < s_; ++i)
			f[i] = r.dot(P_[i]);

		for (int k = 0; k < s_; ++k)
		{			
			
			// Update inner iteration counter
			ii = ii + 1;
			
			// Compute new v (copy construction)
			v = r;
			
			if (jj > 0)
			{				
				// Solve small lower triangular system and make
				// v orthogonal to P:
				for (int i = k; i < s_; ++i)
				{
					gamma[i] = f[i];
					for (int j = k; j < i; ++j)
					{
						gamma[i] = gamma[i] - M[i][j] * gamma[j];
					}
					gamma[i] = gamma[i] / M[i][i];
					v.update(-gamma[i], G[i], 1.0);
				}
				
				// Preconditioning
				//--> needs checks, now we are assuming
				//    we do right preconditioning
				model_.applyPrecon(v, t);
				t.scale(om);

				// Compute new U(:,k)
				for (int i = k; i < s_; ++i)
					t.update(gamma[i], U[i], 1.0);
				U[k] = t;

				// Compute Hessenberg matrix
				//-->TODO
			}
			else if (createInitU) 
			{
				// create initial search space
				model_.applyPrecon(v, U_init_[k]);
				U[k] = U_init_[k];
			}

			// Compute new G(:,k), G(:,k) is in space G_j
			model_.applyMatrix(U[k], G[k]);

			// Bi-Orthogonalise the new basis vectors:
			double alpha;
			for (int i = 0; i < k; ++i)
			{
				alpha = P_[i].dot(G[k]) / M[i][i];
				G[k].update(-alpha, G[i], 1.0);
				U[k].update(-alpha, U[i], 1.0);
			}

			// New column of M = P'*G  (first k-1 entries are zero)
			for (int i = k; i < s_; ++i)
			{
				M[i][k] = G[k].dot(P_[i]);
			}
			if (M[k][k] == 0)
			{
				flag = 3;
				std::cout << "WARNING: M[k][k] == 0"  << std::endl;
				return -1;
			}

			// Make r orthogonal to p_i, i = 1..k, update solution and residual
			double beta = f[k] / M[k][k];
			r.update( -beta, G[k], 1.0);  // r = r - beta*G(:,k);
			x_->update(beta, U[k], 1.0);  // x = x + beta*U(:,k);

			// Check whether we need to replace residual
			normr_ = r.norm();
			if (replacement_ && normr_ > tolb_ / mp_) trueres_ = true;

			// Smoothing
			//-->TODO

			// New f = P'r (first k  components are zero)
			//-->Here we use s-1 since we are zero-based!
			if (k < s_-1)
			{
				for (int i = k+1; i < s_; ++i)
				{
					f[i] = f[i] - beta*M[i][k];
				}
			}

			// Check for convergence
			normr_ = r.norm();
			resvec_.push_back(normr_);
			iter_ = iter_ + 1;
			if (verbosity_ > 4) printIterStatus();
			if (normr_ < tolb_ || iter_ == maxit_)
				break;
		}
		
		// Check for convergence
		if (normr_ < tolb_ || iter_ == maxit_)
			break;

		// Update G-space counter
		jj = jj + 1;

		// Compute first residual in G_j+1
		// Note: r is already perpendicular to P so v = r
		// Preconditioning: v = inv(P)r
		model_.applyPrecon(r, v);
		
		// Matrix-vector multiplication: t = A*v
		model_.applyMatrix(v, t);

		// Computation of a new omega
		om = calc_omega(t, r, angle_);
		
		// Update solution and residual:
		r.update(-om,  t, 1.0);   // r = r - om*t
		x_->update(om, v, 1.0);   // x = x + om*v 

		normr_ = r.norm();
		
		// Residual replacement?
		if (replacement_ && normr_ > tolb_ / mp_) trueres_ = true;

		if (trueres_ && normr_ < normb )
		{
			std::cout << "IDR: replacing residual..." << std::endl;
			// r = b - A*x;
			model_.applyMatrix(*x_, r);   // Ax
			r.update( 1.0,  *b_, -1.0);   // b - Ax
			trueres_ = false;
			replacements_ = replacements_ + 1;
		}

		// Smoothing
		//-->TODO

		resvec_.push_back(normr_);
		iter_ = iter_ + 1;		
		if (verbosity_ > 4) printIterStatus();
	}

	// Smoothing
	//-->TODO
	
	return flag;
}

//====================================================================
template<typename Model, typename VectorPointer>
void IDRSolver<Model, VectorPointer>::
printIterStatus()
{
	std::cout << "iteration: "   << iter_
			  << " residual: "   << normr_
			  << " rel. tol: "   << tolb_
			  << std::endl;
}
	
//====================================================================
template<typename Model, typename VectorPointer>
void IDRSolver<Model, VectorPointer>::
printResVec()
{
	int ctr = 0;
	for (auto &i: resvec_)
		std::cout << ctr << " " << i << std::endl;
}

//====================================================================
template<typename Model, typename VectorPointer>
double IDRSolver<Model, VectorPointer>::
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
template<typename Model, typename VectorPointer>
void IDRSolver<Model, VectorPointer>::
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

//====================================================================
template<typename Model, typename VectorPointer>
double IDRSolver<Model, VectorPointer>::
explicitResNorm()
{
	Vector r(*x_);
	model_.applyMatrix(*x_, r); // Ax
	r.update( 1.0,  *b_, -1.0); // b - Ax
	return r.norm() / b_->norm();
}

//====================================================================
template<typename Model, typename VectorPointer>
void IDRSolver<Model, VectorPointer>::
writeVector(std::vector<double> &vector, const std::string &filename)
{
	if (!vector.empty())
	{
		std::ofstream atmos_ofstream;
		atmos_ofstream.open(filename);
		for (auto &i : vector)
			atmos_ofstream << std::setprecision(12) << i << '\n';
		atmos_ofstream.close();
	}
	else
		WARNING("vector is empty", __FILE__, __LINE__);
}

#endif
