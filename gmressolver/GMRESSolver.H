#ifndef GMRESSOLVER_H
#define GMRESSOLVER_H

#include "GMRESSolverDecl.H"
#include "GMRESMacros.H"
#include <vector>


//====================================================================
// constructor 1
template<typename Model, typename VectorPointer>
GMRESSolver<Model, VectorPointer>::
GMRESSolver(Model &model)
	:
	model_           (model),  // model
	haveInitSol_     (false),
	haveRHS_         (false),
	tol_             (1e-8),
	resid_           (1.0),
	maxit_           (1000),
	m_               (1000),
	iter_            (0),
	prec_            (true),
	leftPrec_        (false),
	verbosity_       (0)
{}

//====================================================================
// constructor 2
template<typename Model, typename VectorPointer>
GMRESSolver<Model, VectorPointer>::
GMRESSolver(Model &model,
			VectorPointer x0,
			VectorPointer b)
	:
	model_           (model),  // model
	x_               (x0),     // initial guess
	b_               (b),      // RHS
	haveInitSol_     (true),
	haveRHS_         (true),
	tol_             (1e-8),
	resid_           (1.0),
	maxit_           (1000),
	m_               (1000),
	iter_            (0),
	prec_            (true),
	leftPrec_        (false),
	verbosity_       (0)
{}

//====================================================================
// destructor
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
		PRINT("GMRES: updating parameters");
	
	tol_            = pars->get("GMRES tolerance"       , tol_);
	maxit_          = pars->get("GMRES iterations"      , maxit_);
	m_              = pars->get("GMRES restart"         , m_);
	verbosity_      = pars->get("GMRES verbosity"       , verbosity_);
	prec_           = pars->get("GMRES preconditioning" , prec_);
	leftPrec_       = pars->get("GMRES left prec"       , leftPrec_);
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
	
	STLVector s (m_+1, 0.0);
	STLVector cs(m_+1, 0.0);
	STLVector sn(m_+1, 0.0);

	Matrix H(m_+1, STLVector(m_, 0.0));
	
	Vector tmp(*x_);
	Vector r  (*x_);
	Vector w  (*x_);

	double normb = b_->norm();

	if (prec_ && leftPrec_)
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
			if ((verbosity_ > 2 && iter_ % 10 == 0) || verbosity_ > 7) printIterStatus();
			
			// Compute w
			if (leftPrec_ && prec_)
			{
 				model_.applyMatrix(V[i], tmp); 
				model_.applyPrecon(tmp, w);     // M.solve(A * v[i])
			}
			else if (prec_) // Right preconditioning (default)
			{
				model_.applyPrecon(V[i], tmp);
				model_.applyMatrix(tmp, w);    // w =  A * M^{-1} v[i]
			}
			else
			{
				model_.applyMatrix(V[i], w);   // w =  A * v[i]
			}
			
			for (k = 0; k <= i; k++)
			{
				H[k][i] = w.dot(V[k]);            // H(k, i) = dot(w, v[k]);
				w.update(-H[k][i], V[k], 1.0);    // w -= H(k, i) * v[k];
			}			
			
			H[i+1][i] = w.norm();
			w.scale(1.0 / H[i+1][i]); //  w / H(i+1, i)
			V[i+1]      = w;
			
			for (k = 0; k < i; k++)
				ApplyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);
      
			GeneratePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
			
			resid_ = std::abs(s[i+1]) / normb;
			if (resid_ < tol_)
			{
				if (verbosity_ > 0)
				{
					PRINT("GMRES implicit residual passed...");
					PRINT("         iterations = " << iter_);
					PRINT("  implicit residual = " << resid_);
				}
				break;
			}
		}		
		
		if (i == m_)
			Update(0, m_ - 1, H, s, V);
		else
			Update(0, i, H, s, V);
		
		// Calculate explicit residual
		if (prec_ && leftPrec_)
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
		beta    = r.norm();
		resid_  = beta / normb;
		
		if (verbosity_ > 0)
			PRINT("  explicit residual = " << resid_);

		if (resid_ < tol_)
		{
			if (verbosity_ > 0)
				PRINT("GMRES explicit residual passed...");
			return 0;
		}
		else if (verbosity_ > 0)
			PRINT("   restart ");
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
Update(int first, int last, Matrix &H, STLVector &s, std::vector<Vector> &V)
{
	STLVector y(s);
	// Backsolve:  
	for (int i = last; i >= first; i--)
	{
		y[i] /= H[i][i];
		for (int j = i - 1; j >= first; j--)
			y[j] -= H[j][i] * y[i];
	}
	
	if (!prec_ || leftPrec_)
	{
		for (int j = first; j <= last; j++)
			x_->update(y[j], V[j], 1.0); //x += v[j] * y(j);
	}
	else // Right preconditioning (default)
	{
		Vector tmp1(*x_);
		Vector tmp2(*x_);
		tmp1.zero(); // tmp1 = 0
		for (int j = first; j <= last; j++)
			tmp1.update(y[j], V[j], 1.0); // tmp1 += v[j] * y(j);
		model_.applyPrecon(tmp1, tmp2);   // tmp2  = inv(M)*V*y
		x_->update(1.0, tmp2, 1.0);		  // x += inv(M)*V*y
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
	PRINT("iteration: "   << iter_  << " residual: "   << resid_
		  << " rel. tol: "   << tol_);
}

#endif
