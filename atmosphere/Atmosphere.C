#include "Atmosphere.H"
#include <math.h>
#include <iostream>
#define PI 3.14159265358979323846

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

	np_  = 27; // all neighbouring points including the center
	nun_ = 1;  // only temperature

	// Construct dependency grid:
	Al_ = std::make_shared<DependencyGrid>(n_, m_, l_, np_, nun_);

	xmin_ = 286 * PI / 180;
	xmax_ = 350 * PI / 180;
	ymin_ = 10  * PI / 180;
	ymax_ = 74  * PI / 180;

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
		datc_.push_back(0.9 + 1.5 * exp(-12 * yc_[j] * yc_[j] / PI));
		datv_.push_back(0.9 + 1.5 * exp(-12 * yv_[j] * yv_[j] / PI));
		suna_.push_back(As_*(1 - .482 * (3 * pow(sin(yc_[j]), 2) - 1.) / 2.) *
						(1 - albe_[j]));
	}
}

Atmosphere::~Atmosphere()
{}


void Atmosphere::fillDependencyGrid()
{
	//! ------------------------------------------------------------------
	//! atmosphere layer
	//! ------------------------------------------------------------------
	DependencyGrid tc (n_, m_, l_, np_, 1);
	DependencyGrid tc2(n_, m_, l_, np_, 1);
	DependencyGrid txx(n_, m_, l_, np_, 1);
	DependencyGrid tyy(n_, m_, l_, np_, 1);

	
	discretize(1, tc);
	/*
	discretize(2, tc2);
	discretize(3, txx);
	discretize(4, tyy);
	*/
}

void computeRHS()
{
	
}

void Atmosphere::discretize(int type, DependencyGrid &grid)
{
	switch (type)
	{
	case 1:
		for (int i = 1; i != n_+1; ++i)
			for (int j = 1; j != m_+1; ++j)
				for (int k = 1; k != l_+1; ++k)
				{
					grid.set(i,j,k,5, 1,1, -1.0);
					grid.set(i,j,k,14,1,1,  1.0);
				}
	}
}

void Atmosphere::test()
{
	std::cout << "Atmosphere test:" << std::endl;
	std::cout << "n: " << n_ << ", m: " << m_ << std::endl;
	std::cout << "xmin: " << xmin_ << ", xmax: " << xmax_ << std::endl;
	std::cout << "ymin: " << ymin_ << ", ymax: " << ymax_ << std::endl;
	std::cout << "length xu: " << xu_.size() << std::endl;
	std::cout << "length xc: " << xc_.size() << std::endl;
	std::cout << "length yv: " << yv_.size() << std::endl;
	std::cout << "length yc: " << yc_.size() << std::endl;
	std::cout << "xu:" << std::endl;
	for (int i = 0; i != xu_.size(); ++i)
		std::cout << xu_[i] << " ";
	std::cout << std::endl;
	std::cout << "Al[1][1][1][1][1][1]="
			  << Al_->get(1,1,1,1,1,1) << std::endl;
	std::cout << "Al[1][1][1][2][1][1]="
			  << Al_->get(1,1,1,2,1,1) << std::endl;
	std::cout << "Al[1][1][1][3][1][1]="
			  << Al_->get(1,1,1,3,1,1) << std::endl;
	std::cout << "Al[1][2][1][1][1][1]="
			  << Al_->get(1,2,1,1,1,1) << std::endl;
}

//==================================================================
// DependencyGrid implementation
//==================================================================
DependencyGrid::DependencyGrid(int n, int m, int l, int np, int nun)
	:
	n_(n),
	m_(m),
	l_(l),
	np_(np),
	nun_(nun)	
{
	// Build dependency grid
	grid_ = new double*****[n_+1];
	for (int i = 1; i != n_+1; ++i)
	{
		grid_[i] = new double****[m_+1];
		for (int j = 1; j != m_+1; ++j)
		{
			grid_[i][j] = new double***[l_+1];
			for (int k = 1; k != l_+1; ++k)
			{
				grid_[i][j][k] = new double**[np_+1];
				for (int loc = 1; loc != np_+1; ++loc)
				{
					grid_[i][j][k][loc] = new double*[nun_+1];
					for (int A = 1; A != nun_+1; ++A)
					{
						grid_[i][j][k][loc][A] = new double[nun_+1];
						for (int B = 1; B != nun_+1; ++B)
						{
							grid_[i][j][k][loc][A][B] = 0.0;
						}
					}
				}
			}
		}
	}
}

DependencyGrid::~DependencyGrid()
{
	// Destroy dependency grid
 	for (int i = 1; i != n_+1; ++i)
	{
		for (int j = 1; j != m_+1; ++j)
		{
			for (int k = 1; k != l_+1; ++k)
			{
				for (int loc = 1; loc != np_+1; ++loc)
				{
					for (int A = 1; A != nun_+1; ++A)
					{
						delete grid_[i][j][k][loc][A];
					}
					delete grid_[i][j][k][loc];
				}
				delete grid_[i][j][k];
			}
			delete grid_[i][j];
		}
		delete grid_[i];
	}
	delete grid_;
}

double DependencyGrid::get(int i, int j, int k, int loc, int A, int B)
{
	return grid_[i][j][k][loc][A][B];
}

void DependencyGrid::set(int i, int j, int k, int loc, int A, int B, double value)
{
	grid_[i][j][k][loc][A][B] = value;
}
