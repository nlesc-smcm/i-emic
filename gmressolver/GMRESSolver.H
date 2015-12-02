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
	model_           (model),
	haveInitSol_     (false),
	haveRHS_         (false),
	replaceResidual_ (false),
	minimizeScheme_  ('B'),
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
GMRESSolver(Model &model, VectorPointer x0, VectorPointer b)
	:
	GMRESSolver(model)
{
	x_ = x0;
	b_ = b;
}

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
	PRINT("GMRES: updating parameters", verbosity_);
	
	tol_             = pars->get("GMRES tolerance"       , tol_);
	maxit_           = pars->get("GMRES iterations"      , maxit_);
	m_               = pars->get("GMRES restart"         , m_);
	verbosity_       = pars->get("GMRES verbosity"       , verbosity_);
	prec_            = pars->get("GMRES preconditioning" , prec_);
	leftPrec_        = pars->get("GMRES left prec"       , leftPrec_);
	replaceResidual_ = pars->get("GMRES replace residual", replaceResidual_);
	minimizeScheme_  = pars->get("GMRES minimizer scheme", minimizeScheme_);
}

// Lapack least squares solver:
extern "C" void dgels_(char* TRANS, int *M, int *N, int *NRHS, double *A,
					   int *LDA, double *B, int *LDB, double *WORK, int *LWORK,
					   int *INFO);

extern "C" void dgesv_(int *N, int *NRHS, double *A,
					   int *LDA, int *IPIV, double *B,
					   int *LDB, int *INFO);

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

	Vector tmp    (*x_);
	Vector r      (*x_);
	Vector w      (*x_);

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
	double proj;
	while (iter_ <= maxit_)
	{
		int first = 0;
		xCopy_ = *x_; // initialize local copy
		beta   = r.norm();
		r.scale(1.0 / beta);
		V[0] = r;
		s.assign(m_+1, 0.0);
		s[0] = beta;

		// We do not exit from this loop when iter_ > maxit,
		// whether this is a bug or a feature I'm not sure...
		for (i = 0; i < m_; i++, iter_++)
		{
			if ((verbosity_ > 2 && !(iter_ % 10)) || verbosity_ > 7)
				printIterStatus();
			
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
			
			if (minimizeScheme_ == 'B')
			{
				for (k = 0; k < i; k++)
					ApplyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);
      
				GeneratePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
				ApplyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
				ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
				
			
				resid_ = std::abs(s[i+1]) / normb;
			}
			
			// if necessary replace residual with true residual
			if (replaceResidual_)
			{
				UpdateCopy(0, i, H, s, V);
				first = i+1;
				
				if (prec_ && leftPrec_)
				{
					model_.applyMatrix(*x_, tmp); // Ax
					tmp.update(1.0, *b_, -1.0);   // b - Ax
					model_.applyPrecon(tmp, r);   // r = M.solve(b - A * x);
				}
				else
				{
					model_.applyMatrix(xCopy_, r);   // A*x
					r.update(1.0, *b_, -1.0);        // b - A*x
				}
				
				beta   = r.norm();                   // ||b - A*x||
				PRINT("      explicit residual = " << beta / normb, verbosity_);
				r.update(.999, V[i+1], .001);
				for (k = 0; k <= i; k++) // orthogonalize r w.r.t. V[0:i]
				{
					proj = r.dot(V[k]);          // pr = dot(r, v[k]);
					r.update(-proj, V[k], 1.0);  // r -= pr * v[k];
				}			
				beta = r.norm();
				r.scale(1.0 / beta);

				V[i+1] = r;
			}
			
			if (resid_ < tol_)
			{
				PRINT("GMRES implicit residual passed...", verbosity_);
				PRINT("         iterations = " << iter_, verbosity_);
				PRINT("  implicit residual = " << resid_, verbosity_);
				break;
			}
		}		
		
		if (replaceResidual_)
			*x_ = xCopy_;
		else if (i == m_)
			Update(0, m_-1, H, s, V);
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
		
		PRINT("  explicit residual = " << resid_, verbosity_);

		if (resid_ < tol_)
		{
			PRINT("GMRES explicit residual passed...", verbosity_);
			return 0;
		}
		else
			PRINT("   restart ", verbosity_);
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
	
	if (minimizeScheme_ == 'B') // Backsolve:
		backsolve(last, H, y);	
	else if (minimizeScheme_ == 'Q')   // Use lapack QR solve
		LLSsolve(last, H, y);

	if (!prec_ || leftPrec_)
	{
		for (int j = first; j <= last; j++)
			x_->update(y[j], V[j], 1.0);  // x += v[j] * y(j);
	}
	else // Right preconditioning (default)
	{
		Vector tmp1(*x_);
		Vector tmp2(*x_);
		tmp1.zero(); // tmp1 = 0
		for (int j = first; j <= last+1; j++)
			tmp1.update(y[j], V[j], 1.0); // tmp1 += v[j] * y(j);
		model_.applyPrecon(tmp1, tmp2);   // tmp2  = inv(M)*V*y
		x_->update(1.0, tmp2, 1.0);		  // x += inv(M)*V*y
	}
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
UpdateCopy(int first, int last, Matrix &H, STLVector &s, std::vector<Vector> &V)
{
	xCopy_ = *x_; 	// Copy solution into xCopy
	
	STLVector y(s);

	if (minimizeScheme_ == 'B') // Backsolve:
		backsolve(last, H, y);	
	else if (minimizeScheme_ == 'Q')   // Use lapack QR solve
		LLSsolve(last, H, y);
	
	if (!prec_ || leftPrec_)
	{
		for (int j = first; j <= last; j++)
			xCopy_.update(y[j], V[j], 1.0); //x += v[j] * y(j);
	}
	else // Right preconditioning (default)
	{
		Vector tmp1(*x_);
		Vector tmp2(*x_);
		tmp1.zero(); // tmp1 = 0
		//PRINT("  updating entries " << first << " to " << last, verbosity_ - 5);
		for (int j = first; j <= last; j++)
			tmp1.update(y[j], V[j], 1.0); // tmp1 += v[j] * y(j);
		model_.applyPrecon(tmp1, tmp2);   // tmp2  = inv(M)*V*y
		xCopy_.update(1.0, tmp2, 1.0);		  // x += inv(M)*V*y
	}
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
backsolve(int m, Matrix &H, STLVector &s)
{
	for (int i = m; i >= 0; i--)
	{
		s[i] /= H[i][i];
		for (int j = i - 1; j >= 0; j--)
			s[j] -= H[j][i] * s[i];
	}
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
LLSsolve(int m, Matrix &H, STLVector &s)
{
	STLVector exportH((m_+2)*(m_+1), 0.0);

	int idx = 0;
	for (int col = 0; col <= m; ++col)
		for (int row = 0; row <= m+1; ++row)
		{
			exportH[idx] = H[row][col];
			idx++;
		}

	char trans = 'N';
	int  M     =  m+1;
	int  N     =  M;
	int  nrhs  =  1;
	int  lda   =  M;
	int  ldb   =  M;
	int  lwork =  4*N;  // some size needed by lapack
	int  ipiv[N];

	double  work[lwork];
	
	int  info;

	// lapack LLS solver
	dgels_(&trans, &M, &N, &nrhs, &exportH[0], &lda,
		   &s[0], &ldb, work, &lwork, &info);

//	dgesv_(&M, &nrhs, &exportH[0], &lda, ipiv, &s[0], &ldb, &info);		   
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
		  << " rel. tol: "   << tol_, verbosity_);
}

#endif
