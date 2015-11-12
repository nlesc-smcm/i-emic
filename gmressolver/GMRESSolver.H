#ifndef GMRESSOLVER_H
#define GMRESSOLVER_H

#include "GMRESSolverDecl.H"
#include <vector>

//====================================================================
template<typename Model, typename VectorPointer>
GMRESSolver<Model, VectorPointer>::
GMRESSolver(Model &model)
	:
	model_         (model),  // model
	haveInitSol_   (false),
	haveRHS_       (false),
	tol_           (1e-8),
	resid_         (1.0),
	maxit_         (1000),
	m_             (1000),
	iter_          (0),
	leftPrec_      (false),
	verbosity_     (0)
{}

//====================================================================
template<typename Model, typename VectorPointer>
GMRESSolver<Model, VectorPointer>::
GMRESSolver(Model &model,
			VectorPointer x0,
			VectorPointer b)
	:
	model_         (model),  // model
	x_             (x0),     // initial guess
	b_             (b),      // RHS
	haveInitSol_   (true),
	haveRHS_       (true),
	tol_           (1e-8),
	resid_         (1.0),
	maxit_         (1000),
	m_             (1000),
	iter_          (0),
	leftPrec_      (false),
	verbosity_     (0)
{}

//====================================================================
template<typename Model, typename VectorPointer>
GMRESSolver<Model, VectorPointer>::
~GMRESSolver() {}


//*****************************************************************************
template<typename Model, typename VectorPointer>
template<typename ParListPtr>
void GMRESSolver<Model, VectorPointer>::
setParameters(ParListPtr pars)
{
	if (verbosity_ > 5)
		std::cout << "GMRES: updating parameters" << std::endl;
	
	tol_         = pars->get("GMRES tolerance", 1e-8);
	maxit_       = pars->get("GMRES iterations", 500);
	m_           = pars->get("GMRES restart", maxit_);
	verbosity_   = pars->get("GMRES verbosity", 0);
	leftPrec_    = pars->get("GMRES left prec", false);
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
int GMRESSolver<Model, VectorPointer>::
solve()
{
	if (!haveInitSol_ || !haveRHS_)
	{
		std::cout << "Problem not setup correctly!"
				  << haveInitSol_ << " "
				  << haveRHS_ << std::endl;
		return 1;
	}
	int i, k;
	iter_ = 0;

	// use zero initial guess
	//x_->zero();
	
	STLVector s (m_+1, 0.0);
	STLVector cs(m_+1, 0.0);
	STLVector sn(m_+1, 0.0);

	Matrix H(m_+1, STLVector(m_, 0.0));
	
	Vector tmp(*x_);
	Vector r  (*x_);
	Vector w  (*x_);

	double normb = b_->norm();

	if (leftPrec_)
	{
		model_.applyMatrix(*x_, tmp); // Ax
		tmp.update(1.0, *b_, -1.0);   // b - Ax
		model_.applyPrecon(tmp, r);   // r = M.solve(b - A * x);
	}
	else
	{
		model_.applyMatrix(*x_, r); // Ax
		r.update(1.0, *b_, -1.0);   // b - Ax
	}
	
	double beta = r.norm();
	
	if (normb == 0.0)
		normb = 1;

	resid_ = beta / normb;
	if (resid_ <= tol_)
	{
		iter_ = 0;
		return 0;
	}

	std::vector<Vector> V(m_+1, Vector());
	
	while (iter_ <= maxit_)
	{
		r.scale(1.0 / beta);
		V[0] = r;
		s.assign(m_+1, 0.0);
		s[0] = beta;
		
		for (i = 0; i < m_ && iter_ <= maxit_; i++, iter_++)
		{
			if (verbosity_ > 4) printIterStatus();
			
			if (leftPrec_)
			{
				model_.applyMatrix(V[i], tmp); 
				model_.applyPrecon(tmp, w);     // M.solve(A * v[i])
			}
			else
			{
				model_.applyPrecon(V[i], tmp);
				model_.applyMatrix(tmp, w);    // w =  A * M^{-1} v[i]
			}
			
			for (k = 0; k <= i; k++)
			{
				H[k][i] = w.dot(V[k]);            // H(k, i) = dot(w, v[k]);
				w.update(-H[k][i], V[k], 1.0);    // w -= H(k, i) * v[k];
			}
			H[i+1][i] = w.norm();
			w.scale(1.0 / H[i+1][i]); //  w / H(i+1, i)
			V[i+1]    = w;  
			
			for (k = 0; k < i; k++)
				ApplyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);
      
			GeneratePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
			
			resid_ = std::abs(s[i+1]) / normb;
			if (resid_ < tol_)
			{
				Update(i, H, s, V);
				return 0;
			}
		}
		Update(m_ - 1, H, s, V);

		if (leftPrec_)
		{
			model_.applyMatrix(*x_, tmp); // Ax
			tmp.update(1.0, *b_, -1.0);   // b - Ax
			model_.applyPrecon(tmp, r);   // r = M.solve(b - A * x);
		}
		else
		{
			model_.applyMatrix(*x_, r); // Ax
			r.update(1.0, *b_, -1.0);   // b - Ax
		}
			
		beta = r.norm();
		
		resid_ = beta / normb;
		if (resid_ < tol_)
		{
			return 0;
		}
	}

	return 1;
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
	if (dy == 0.0)
	{
		cs = 1.0;
		sn = 0.0;
	}
	else if (abs(dy) > abs(dx))
	{
		double temp = dx / dy;
		sn = 1.0 / sqrt( 1.0 + temp*temp );
		cs = temp * sn;
	}
	else
	{
		double temp = dy / dx;
		cs = 1.0 / sqrt( 1.0 + temp*temp );
		sn = temp * cs;
	}
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
	double temp  =  cs * dx + sn * dy;
	dy = -sn * dx + cs * dy;
	dx = temp;
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
Update(int k, Matrix &H, STLVector &s, std::vector<Vector> &V)
{
	STLVector y(s);
	// Backsolve:  
	for (int i = k; i >= 0; i--)
	{
		y[i] /= H[i][i];
		for (int j = i - 1; j >= 0; j--)
			y[j] -= H[j][i] * y[i];
	}
	
	if (leftPrec_)
	{
		for (int j = 0; j <= k; j++)
			x_->update(y[j], V[j], 1.0); //x += v[j] * y(j);
	}
	else
	{
		Vector tmp1(*x_);
		Vector tmp2(*x_);
		tmp1.zero(); // tmp1 = 0
		for (int j = 0; j <= k; j++)
			tmp1.update(y[j], V[j], 1.0); // tmp1 += v[j] * y(j);
		model_.applyPrecon(tmp1, tmp2);   // tmp2  = inv(M)*V*y
		x_->update(1.0, tmp2, 1.0);			
	}
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
double GMRESSolver<Model, VectorPointer>::
residual()
{
	return resid_;
}

//====================================================================
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
printIterStatus()
{
	std::cout << "iteration: "   << iter_
			  << " residual: "   << resid_
			  << " rel. tol: "   << tol_
			  << std::endl;
}

#endif
