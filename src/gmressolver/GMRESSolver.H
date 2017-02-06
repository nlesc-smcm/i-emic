#ifndef GMRESSOLVER_H
#define GMRESSOLVER_H

#include "GMRESSolverDecl.H"
#include "GMRESMacros.H"
#include "GlobalDefinitions.H"
#include <vector>
#include <math.h>

//====================================================================
// constructor 1
template<typename Model, typename VectorPointer>
GMRESSolver<Model, VectorPointer>::
GMRESSolver(Model &model)
	:
	model_           (model),
	haveInitSol_     (false),
	haveRHS_         (false),
	minimizeScheme_  ('B'),
	flexible_        (true),
	computeExplResid_(false),
	tol_             (1e-4),
	resid_           (1.0),
	maxit_           (500),
	m_               (400),
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
	
	tol_              = pars->get("GMRES tolerance"        , tol_);
	maxit_            = pars->get("GMRES iterations"       , maxit_);
	m_                = pars->get("GMRES restart"          , m_);
	verbosity_        = pars->get("GMRES verbosity"        , verbosity_);
	prec_             = pars->get("GMRES preconditioning"  , prec_);
	leftPrec_         = pars->get("GMRES left prec"        , leftPrec_);
	minimizeScheme_   = pars->get("GMRES minimizer scheme" , minimizeScheme_);
	flexible_         = pars->get("GMRES flexible"         , flexible_);
	computeExplResid_ = pars->get("GMRES explicit residual", computeExplResid_);
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
		model_.applyPrecon(tmp, r);   // r = inv(M) * (b - A * x);
	}
	else
	{
		model_.applyMatrix(*x_, r); // Ax
		r.update(1.0, *b_, -1.0);   // b - Ax
	}
	
	double beta = r.norm();
	
	if (normb == 0.0)
		normb = 1;

	resid_     = beta / normb;
	explResid_ = resid_;
	if (resid_ <= tol_)
	{
		iter_ = 0;
		return 0;
	}
	
	std::vector<Vector> Z(m_+1, Vector()); // for FGMRES
	std::vector<Vector> V(m_+1, Vector()); 
	int spaceSize;	// keeping track of the size of Z,V
	
	while (iter_ <= maxit_)
	{
		xCopy_ = *x_; // initialize local copy
		beta   = r.norm();
		r.scale(1.0 / beta);
		
		V[0] = r;
		
		s.assign(m_+1, 0.0);
		s[0] = beta;

		for (i = 0; i < m_ && iter_ <= maxit_; i++, iter_++)
		{
			if ((verbosity_ > 2 && !(iter_ % 10)) || verbosity_ > 7)
				printIterStatus();
			
			// Compute w
			TIMER_START("GMRES: compute w...");
			if (prec_ && leftPrec_)               // Left preconditioning
			{
 				model_.applyMatrix(V[i], tmp); 
				model_.applyPrecon(tmp, w);       // inv(M) * (A * v[i])
			}
			else if (prec_ && flexible_)          // Right preconditioned FGMRES
			{
				model_.applyPrecon(V[i], tmp);    // z[i] = M^{-1} * v[i]
				Z[i] = tmp;
				model_.applyMatrix(tmp, w);       // w    = A * z[i]
			}
			else if (prec_)                       // Right preconditioning (default)
			{
				model_.applyPrecon(V[i], tmp);				
				model_.applyMatrix(tmp, w);       // w =  A * M^{-1} * v[i]
			}
			else
			{
				model_.applyMatrix(V[i], w);      // w =  A * v[i]
			}
			TIMER_STOP("GMRES: compute w...");

			// Orthogonalize w
			TIMER_START("GMRES: orthogonalization...");
			for (k = 0; k <= i; k++)
			{
				H[k][i] = w.dot(V[k]);            // H(k, i) = dot(w, v[k]);
				w.update(-H[k][i], V[k], 1.0);    // w -= H(k, i) * v[k];
			}
			TIMER_STOP("GMRES: orthogonalization...");

			// Normalize and assign to space
			H[i+1][i] = w.norm();
			w.scale(1.0 / H[i+1][i]);             //  w / H(i+1, i)			
			V[i+1]    =  w;
			spaceSize = i;
			
			if (minimizeScheme_ == 'B')
			{
				for (k = 0; k < i; k++)
					ApplyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);
				
				GeneratePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
				ApplyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
				ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
				
				resid_ = std::abs(s[i+1]) / normb;
			}
			else 
				resid_ =  compute_r(i, H, s) / normb;

			if (computeExplResid_)
			{
				if (flexible_)
					explResid_ = compute_explicit_residual(i, H, s, Z) / normb;
				else
					explResid_ = compute_explicit_residual(i, H, s, V) / normb;

				resid_ = std::max(resid_, explResid_);
			}
			
			if (resid_ < tol_)
			{
				PRINT("GMRES residual passed...", verbosity_);
				PRINT("           iterations = " << iter_, verbosity_);
				PRINT("             residual = " << resid_, verbosity_);
				break;
			}
		}		

		// Update solution 
		if (flexible_)
			Update(spaceSize, H, s, Z); // xm = x0 + Z*ym
		else
			Update(spaceSize, H, s, V); // xm = x0 + inv(M)*V*ym
		
		// Calculate explicit residual
		if (prec_ && leftPrec_)
		{
			model_.applyMatrix(*x_, tmp); // Ax
			tmp.update(1.0, *b_, -1.0);   // b - Ax
			model_.applyPrecon(tmp, r);   // r = inv(M) * (b - A * x);
		}
		else
		{
			model_.applyMatrix(*x_, r); // Ax
			r.update(1.0, *b_, -1.0);   // b - Ax
		}
		beta   = r.norm();
		explResid_ = beta / normb;
		
		PRINT("    true residual = " << explResid_, verbosity_);

		if (resid_ < tol_)
		{
			PRINT("GMRES explicit residual passed...", verbosity_);
			return 0;
		}
		else
			PRINT("    :( ...  restart ", verbosity_);
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
Update(int last, Matrix &H, STLVector &s, std::vector<Vector> &V)
{
	compute_y(last, H, s);
	
	if (!prec_ || leftPrec_ || flexible_)
	{
		for (int j = 0; j <= last; j++)
			x_->update(y_[j], V[j], 1.0);  // x += v[j] * y(j);
	}
	else    // Right preconditioning (default)
	{
		Vector tmp1(*x_);
		Vector tmp2(*x_);
		tmp1.zero(); // tmp1 = 0
		for (int j = 0; j <= last; j++)
			tmp1.update(y_[j], V[j], 1.0); // tmp1 += v[j] * y(j);
		model_.applyPrecon(tmp1, tmp2);    // tmp2  = inv(M)*V*y
		x_->update(1.0, tmp2, 1.0);		   // x    += inv(M)*V*y
	}
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
compute_y(int last, Matrix &H, STLVector &s)
{
	y_ = s;

	if (minimizeScheme_ == 'B') // Backsolve:
		backSolve(last, H, y_);	
	else if (minimizeScheme_ == 'Q')   // Use lapack QR solve
		LLSSolve(last, H, y_);
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
double GMRESSolver<Model, VectorPointer>::
compute_r(int last, Matrix &H, STLVector &s)
{
	compute_y(last, H, s);
	
	int M = last+2;
	int N = last+1;

	STLVector tmp(M,0.0);

	double nrm = 0;
	
	for (int row = 0; row < M; ++row)
	{
		for (int col = 0; col < N; ++col)
			tmp[row] += H[row][col]*y_[col];
		tmp[row] -= s[row];
	}
	
	for (int row = 0; row < M; ++row)
		nrm += tmp[row] * tmp[row];
	
	nrm = sqrt(nrm);
	return nrm;
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
double GMRESSolver<Model, VectorPointer>::
compute_explicit_residual(int last, Matrix &H, STLVector &s, std::vector<Vector> &V)
{
	// update local copy of solution
	UpdateCopy(last, H, s, V);

	// compute residual
	Vector r(*x_);
	model_.applyMatrix(xCopy_, r);   // A*x
	r.update(1.0, *b_, -1.0);        // b - A*x
	double beta   = r.norm();        // ||b - A*x||
	return beta;
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
UpdateCopy(int last, Matrix &H, STLVector &s, std::vector<Vector> &V)
{
	xCopy_ = *x_; 	// Copy solution into xCopy		

	compute_y(last, H, s);
		
	if (!prec_ || leftPrec_ || flexible_)
	{
		for (int j = 0; j <= last; j++)
			xCopy_.update(y_[j], V[j], 1.0); //x += v[j] * y(j);
	}
	else // Right preconditioning (default)
	{
		Vector tmp1(*x_);
		Vector tmp2(*x_);
		tmp1.zero(); // tmp1 = 0
		//PRINT("  updating entries " << first << " to " << last, verbosity_ - 5);
		for (int j = 0; j <= last; j++)
			tmp1.update(y_[j], V[j], 1.0); // tmp1 += v[j] * y(j);
		model_.applyPrecon(tmp1, tmp2);   // tmp2  = inv(M)*V*y
		xCopy_.update(1.0, tmp2, 1.0);		  // x += inv(M)*V*y
	}
}

//*****************************************************************************
template<typename Model, typename VectorPointer>
void GMRESSolver<Model, VectorPointer>::
backSolve(int m, Matrix &H, STLVector &s)
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
LLSSolve(int last, Matrix &H, STLVector &s)
{
	int  M     =  last+2;  // FOM -> last+1, GMRES last+2 ??
	int  N     =  last+1;
	
	STLVector exportH(M * N, 0.0);

	int idx = 0;
	for (int col = 0; col < N; ++col)
		for (int row = 0; row < M; ++row)
		{
			exportH[idx] = H[row][col];
			idx++;
		}

	if (verbosity_ > 10)
	{
		PRINT("  testing exportH -----", verbosity_);
		for (auto &it: exportH)
			PRINT("   " << it, verbosity_);
		PRINT("  ", verbosity_);
		PRINT("  " << last << "  " << exportH.size(),  verbosity_);
		PRINT("  testing exportH -----", verbosity_);
	}

	// Variables required by lapack:
	char    trans = 'N';
	int     nrhs  =  1;
	int     lda   =  M;
	int     ldb   =  M;
	int     lwork =  4*N;  // some size needed by lapack
	int     info;
	double  work[lwork];
	
	// lapack LLS solver
	dgels_(&trans, &M, &N, &nrhs, &exportH[0], &lda,
		   &s[0], &ldb, work, &lwork, &info);
	
	PRINT("  LLS: DGELS_ info = " << info, verbosity_ - 10);
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
	PRINT("iteration: "   << iter_  << " impl res: "   << resid_
		  << " expl res: "   << explResid_, verbosity_);
}

#endif
