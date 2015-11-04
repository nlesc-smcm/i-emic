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
	maxit_         (1000),
	m_             (1000),
	iter_          (0),
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
	tol_           (1e-6),
	maxit_         (1000),
	m_             (1000),
	iter_          (0),
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
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
int GMRESSolver<Model, VectorPointer>::
solve()
{
	int i, j = 1, k;
	
	STLVector s (m_+1, 0.0);
	STLVector cs(m_+1, 0.0);
	STLVector sn(m_+1, 0.0);

	Matrix H(m_+1, STLVector(m_, 0.0));
	
	Vector tmp(*x_);
	Vector r  (*x_);
	Vector w  (*x_);

	model_.applyPrecon(*b_, tmp);
	double normb = tmp.norm();	  // norm(M.solve(b))

	model_.applyMatrix(*x_, tmp); // Ax
	tmp.update(1.0, *b_, -1.0);   // b - Ax
	model_.applyPrecon(tmp, r);   // r = M.solve(b - A * x);
	double beta = r.norm();

	if (normb == 0.0)
		normb = 1;

	double resid = beta / normb;
	if (resid <= tol_)
	{
		tol_  = resid;
		iter_ = 0;
		return 0;
	}

	std::vector<Vector> V(m+1, Vector());
	
	while (j <= maxit_)
	{
		V[0] = r.scale(1.0 / beta);
		s.assign(m_+1, 0.0);
		s[0] = beta;
		
		for (i = 0; i < m_ && j <= maxit_; i++, j++)
		{
			model_.applyMatrix(V[i], tmp); 
			model_.applyPrecon(tmp, w);     // M.solve(A * v[i])
			
			for (k = 0; k <= i; k++)
			{
				H[k][i] = w.dot(V[k]);            // H(k, i) = dot(w, v[k]);
				w.update(-H[k][i], V[k], 1.0);    // w -= H(k, i) * v[k];
			}
			H[i+1][i] = w.norm();
			V[i+1]    = w.scale(1.0 / H[i+1][i]); //  w / H(i+1, i)
			
			for (k = 0; k < i; k++)
				ApplyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);
      
			GeneratePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

			resid = std::abs(s[i+1]) / normb;
			if (resid < tol)
			{
				Update(i, H, s, V);
				tol   = resid;
				iter_ = j;
				return 0;
			}
		}
		Update(m - 1, H, s, V);
		model_.applyMatrix(*x_, tmp); // Ax
		tmp.update(1.0, *b_, -1.0);   // b - Ax
		model_.applyPrecon(tmp, r);   // r = M.solve(b - A * x);
		beta = r.norm();
		
		resid = beta / normb;
		if (resid < tol)
		{
			tol   = resid;
			iter_ = j;
			return 0;
		}
	}
	
	tol = resid;
	return 1;
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver::GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
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
void GMRESSolver::ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
	double temp  =  cs * dx + sn * dy;
	dy = -sn * dx + cs * dy;
	dx = temp;
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver::Update(int k, Matrix &h, STLVector &s, std::vector<Vector> v)
{
	STLVector y(s);
	// Backsolve:  
	for (int i = k; i >= 0; i--)
	{
		y[i] /= h[i][i];
		for (int j = i - 1; j >= 0; j--)
			y[j] -= h[j,i] * y[i];
	}
	
	for (int j = 0; j <= k; j++)
		x_->update(y[j], V[j], 1.0); //x += v[j] * y(j);
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
double GMRESSolver<Model, VectorPointer>::
residual()
{
	return 666;
}


#endif
