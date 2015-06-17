#include "Atmosphere.H"
#include "AtmosphereDefinitions.H"
#include <math.h>
#include <iostream>

Atmosphere::Atmosphere()
{
	
	// Filling the parameters
	rhoa_   = 1.25       ; //! atmospheric density 
	hdima_  = 8400.      ; //! atmospheric scale height 
	cpa_    = 1000.      ; //! heat capacity 
	d0_     = 3.1e+06    ; //! constant eddy diffusivity 
	arad_   = 216.0      ; //! radiative flux param A
	brad_   = 1.5        ; //! radiative flux param B
	sun0_   = 1360.      ; //! solar constant 
	c0_     = 0.43       ; //! atmospheric absorption coefficient
	ce_     = 1.3e-03    ; //! exchange coefficient 
	ch_     = 0.94 * ce_ ; //! exchange coefficient 
	uw_     = 8.5        ; //! mean atmospheric surface wind speed
	t0_     = 15.0       ; //! reference temperature
	udim_   = 0.1e+00    ; //! typical horizontal velocity of the ocean
	r0dim_  = 6.37e+06   ; //! radius of the earth

	// Filling the coefficients
	muoa_ =  rhoa_ * ch_ * cpa_ * uw_;
	amua_ = (arad_ + brad_ * t0_) / muoa_;
	bmua_ =  brad_ / muoa_;
	Ai_   =  rhoa_ * hdima_ * cpa_ * udim_ / (r0dim_ * muoa_);
	Ad_   =  rhoa_ * hdima_ * cpa_ * d0_ / (muoa_ * r0dim_ * r0dim_);
	As_   =  sun0_ * (1 - c0_) / (4 * muoa_);		

	// Set problem size and domain limits (temporary!)
	n_ = 16;
	m_ = 16;
	l_ = 1;

	np_  = NP_;   // all neighbouring points including the center
	nun_ = NUN_;  // only temperature

	// Construct dependency grid:
	Al_ = std::make_shared<DependencyGrid>(n_, m_, l_, np_, nun_);

	xmin_ = 286 * PI_ / 180;
	xmax_ = 350 * PI_ / 180;
	ymin_ = 10  * PI_ / 180;
	ymax_ = 74  * PI_ / 180;

	// Set the grid increments
	dx_ = (xmax_ - xmin_) / n_;
	dy_ = (ymax_ - ymin_) / m_;

	// Fill x
	for (int i = 0; i != n_+1; ++i)
	{
		xu_.push_back(xmin_ + i * dx_);
		xc_.push_back(xmin_ + (i - 0.5) * dx_);
	}
	
	// Fill y and latitude-based arrays
	for (int j = 0; j != m_+1; ++j)
	{
		yv_.push_back( ymin_ + j * dy_ );
		yc_.push_back( ymin_ + (j - 0.5) * dy_ );
		
		albe_.push_back(0.3);
		datc_.push_back(0.9 + 1.5 * exp(-12 * yc_[j] * yc_[j] / PI_));
		datv_.push_back(0.9 + 1.5 * exp(-12 * yv_[j] * yv_[j] / PI_));
		suna_.push_back(As_*(1 - .482 * (3 * pow(sin(yc_[j]), 2) - 1.) / 2.) *
						(1 - albe_[j]));
	}
}

Atmosphere::~Atmosphere()
{}


void Atmosphere::fillDependencyGrid()
{
	Atom tc (n_, m_, l_, np_);
	Atom tc2(n_, m_, l_, np_);
	Atom txx(n_, m_, l_, np_);
	Atom tyy(n_, m_, l_, np_);

	discretize(1, tc);
	discretize(2, tc2);
	discretize(3, txx);
	discretize(4, tyy);
	
	// Al(:,:,:,:,TT,TT) = - Ad * (txx + tyy) - tc + bmua*tc2
	txx.update(-Ad_, -Ad_, tyy, -1, tc, bmua_, tc2);
	//--> range array should have size 8 here!
	Al_->set({1,n_,1,m_,1,l_,1,np_}, 1, 1, txx);
	
	test(); // ?
}

void computeRHS()
{
	
}

void Atmosphere::discretize(int type, Atom &atom)
{
	switch (type)
	{
		double val2, val4, val5, val6, val8;
	case 1: // tc
		atom.set({1,n_,1,m_,1,l_}, 2,  -1.0);
		atom.set({1,n_,1,m_,1,l_}, 14,  1.0);
		break;
	case 2: // tc2
		atom.set({1,n_,1,m_,1,l_}, 5, 1.0);
		break;
	case 3: // txx
		double cosdx2i;
		for (int i = 1; i != n_+1; ++i)
			for (int j = 1; j != m_+1; ++j)
			{
				cosdx2i = 1.0 / pow(cos(yc_[j]), 2);
				val2 = datc_[j] * cosdx2i;
				val8 = val2;
				val5 = -2 * val2;
				
				for (int k = 1; k != l_+1; ++k)
				{
					atom.set(i, j, k, 2, val2);
					atom.set(i, j, k, 8, val8);
					atom.set(i, j, k, 5, val5);
				}
			}
		break;
	case 4: // tyy
		double dy2i = 1.0 / pow(dy_, 2);
		for (int i = 1; i != n_+1; ++i)
			for (int j = 1; j != m_+1; ++j)
			{
				val4 = dy2i * datv_[j-1] * cos(yv_[j-1]) / cos(yc_[j]);
				val6 = dy2i * datv_[j]   * cos(yv_[j])   / cos(yc_[j]);
				val5 = -(val4 + val6);
				
				for (int k = 1; k != l_+1; ++k)
				{
					atom.set(i, j, k, 4, val4);
					atom.set(i, j, k, 6, val6);
					atom.set(i, j, k, 5, val5);
				}
			}
		break;
	}
}

void Atmosphere::test()
{}


//=============================================================================
// DependencyGrid implementation
//=============================================================================
DependencyGrid::DependencyGrid(int n, int m, int l, int np, int nun)
	:
	n_(n),
	m_(m),
	l_(l),
	np_(np),
	nun_(nun),
	grid_(n, m, l, np, nun, nun)
{
}

//-----------------------------------------------------------------------------
DependencyGrid::~DependencyGrid()
{
}

//-----------------------------------------------------------------------------
double DependencyGrid::get(int i, int j, int k, int loc, int A, int B)
{
	// converting to 0-based
	return grid_(i-1, j-1, k-1, loc-1, A-1, B-1);
}

//-----------------------------------------------------------------------------
void DependencyGrid::set(int i, int j, int k, int loc, int A, int B, double value)
{
	// converting to 0-based
	grid_(i-1, j-1, k-1, loc-1, A-1, B-1) = value;
}

//-----------------------------------------------------------------------------
void DependencyGrid::set(int const (&range)[8], int A, int B, Atom &atom)
{
	for (int i = range[0]; i != range[1]+1; ++i)
		for (int j = range[2]; j != range[3]+1; ++j)
			for (int k = range[4]; k != range[5]+1; ++k)
				for (int loc = range[6]; loc != range[7]+1; ++loc)
				{
					// converting to 0-based
					grid_(i-1, j-1, k-1, loc-1, A-1, B-1) =
						atom.get(i,j,k,loc);
 				}
}

//=============================================================================
// Atom implementation
//=============================================================================
Atom::Atom(int n, int m, int l, int np)
	:
	n_(n),
	m_(m),
	l_(l),
	np_(np),
	atom_(n, m, l, np)
{
}

//-----------------------------------------------------------------------------
Atom::~Atom()
{
}

//-----------------------------------------------------------------------------
double Atom::get(int i, int j, int k, int loc)
{
	// converting to 0-based
	return atom_(i-1, j-1, k-1, loc-1);
}

//-----------------------------------------------------------------------------
void Atom::set(int i, int j, int k, int loc, double value)
{
	// converting to 0-based
	atom_(i-1, j-1, k-1, loc-1) = value;
}

//-----------------------------------------------------------------------------
void Atom::set(int const (&range)[6], int loc, double value)
{
	for (int i = range[0]; i != range[1]+1; ++i)
		for (int j = range[2]; j != range[3]+1; ++j)
			for (int k = range[4]; k != range[5]+1; ++k)
			{
				atom_(i-1, j-1, k-1, loc-1) = value;
			}
}

//-----------------------------------------------------------------------------
// this = scalThis*this+scalA*A+scalB*B+scalC*C
void Atom::update(double scalarThis,
				  double scalarA, Atom &A,
				  double scalarB, Atom &B,
				  double scalarC, Atom &C)
{
	for (int i = 1; i != n_+1; ++i)
		for (int j = 1; j != m_+1; ++j)
			for (int k = 1; k != l_+1; ++k)
				for (int loc = 1; loc != np_+1; ++loc)
				{
					// converting to 0-based
					atom_(i-1, j-1, k-1, loc-1) =
						scalarThis * atom_(i-1, j-1, k-1, loc-1) +
						scalarA*A.get(i, j, k, loc) +
						scalarB*B.get(i, j, k, loc) +
						scalarC*C.get(i, j, k, loc);
				}	
}
