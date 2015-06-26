#include "Atmosphere.H"
#include "AtmosphereDefinitions.H"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vector.H"

//-----------------------------------------------------------------------------
Atmosphere::Atmosphere()
{
	// Continuation parameters
	ampl_    = 0.0        ; //! amplitude of forcing
	amplEnd_ = 1.0        ; //!
	
	// Filling the parameters --> xml
	rhoa_    = 1.25       ; //! atmospheric density 
	hdima_   = 8400.      ; //! atmospheric scale height 
	cpa_     = 1000.      ; //! heat capacity 
	d0_      = 3.1e+06    ; //! constant eddy diffusivity 
	arad_    = 216.0      ; //! radiative flux param A
	brad_    = 1.5        ; //! radiative flux param B
	sun0_    = 1360.      ; //! solar constant 
	c0_      = 0.43       ; //! atmospheric absorption coefficient
	ce_      = 1.3e-03    ; //! exchange coefficient 
	ch_      = 0.94 * ce_ ; //! exchange coefficient 
	uw_      = 8.5        ; //! mean atmospheric surface wind speed
	t0_      = 15.0       ; //! reference temperature
	udim_    = 0.1e+00    ; //! typical horizontal velocity of the ocean
	r0dim_   = 6.37e+06   ; //! radius of the earth

	
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
	nun_ = NUN_;  // only temperature TT_
	
	// Initialize RHS with zeros
	rhs_ = std::make_shared<std::vector<double> >(n_ * m_ * l_, 0.0);

	// Initialize solution
	sol_ = std::make_shared<std::vector<double> >(n_ * m_ * l_, 0.0);

	// Initialize state
	state_ = std::make_shared<std::vector<double> >(n_ * m_ * l_, 0.0);

	// Initialize ocean temperature
	oceanTemp_ = std::vector<double>(n_ * m_, 0.0);
	
	// Initialize forcing with zeros
	frc_ = std::vector<double>(n_ * m_ * l_, 0.0);

	// Initialize dense matrix:
	denseA_ = std::vector<double>(pow(n_ * m_ * l_, 2), 0.0);
	
	// Create pivot array for use in lapack
	ipiv_ = new int[n_*m_*l_+1];
	
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

	// put some values in oceanTemp and the state
	double value;
	int row;
	for (int i = 1; i <= n_; ++i)
		for (int j = 1; j <= m_; ++j)
		{
			value = 10*cos(PI_*(yc_[j]-ymin_)/(ymax_-ymin_));
			row   = find_row(i,j,l_,TT_)-1;
			oceanTemp_[value] = value;
			(*state_)[row] = value;
		}
}

//-----------------------------------------------------------------------------
Atmosphere::~Atmosphere()
{
	delete[] ipiv_;
}

//-----------------------------------------------------------------------------
void Atmosphere::computeJacobian()
{
	Atom tc (n_, m_, l_, np_);
	Atom tc2(n_, m_, l_, np_);
	Atom txx(n_, m_, l_, np_);
	Atom tyy(n_, m_, l_, np_);

	discretize(1, tc);
	discretize(2, tc2);
	discretize(3, txx);
	discretize(4, tyy);
	
	// Al(:,:,:,:,TT,TT) = Ad * (txx + tyy) + tc - bmua*tc2
	txx.update(Ad_, Ad_, tyy, 1, tc, -bmua_, tc2);
	//--> range array should have size 8 here!
	Al_->set({1,n_,1,m_,1,l_,1,np_}, 1, 1, txx);
	boundaries();
	assemble();
}

//-----------------------------------------------------------------------------
void Atmosphere::computeRHS()
{
	// If necessary compute a new Jacobian
   //	if (recomputeJacobian_)
	computeJacobian();

	// Compute the forcing
	forcing();
	
	// Compute the right hand side rhs_	
	double value;
	int row;
	for (int i = 1; i <= n_; ++i)
		for (int j = 1; j <= m_; ++j)
		{
			row   = find_row(i, j, l_, TT_);
			value = matvec(row) + frc_[row-1];
			(*rhs_)[row-1] = value;
		}
}

//-----------------------------------------------------------------------------
double Atmosphere::matvec(int row)
{
	// Returns inner product of a row in the matrix with the state.
	// > ugly stuff with 1 to 0 based...
	int first = beg_[row-1];
	int last  = beg_[row] - 1;
	double result = 0.0;
	for (int j = first; j <= last; ++j)
		result += ico_[j-1] * (*state_)[jco_[j-1]-1];
	
	return result;
}

//-----------------------------------------------------------------------------
void Atmosphere::forcing()
{
	std::cout << "forcing amplitude: " << ampl_ << std::endl;
	double value;
	int row;
	for (int j = 1; j <= m_; ++j)
		for (int i = 1; i <= n_; ++i)
		{
			row = find_row(i, j, l_, TT_);
			value = oceanTemp_[row-1] + ampl_*(suna_[j] - amua_);
			frc_[row-1] = value;
		}
}

//-----------------------------------------------------------------------------
void Atmosphere::discretize(int type, Atom &atom)
{
	switch (type)
	{
		double val2, val4, val5, val6, val8;
	case 1: // tc
		atom.set({1,n_,1,m_,1,l_}, 5,  -1.0);
		// coupling of ocean is more convenient through forcing
		// atom.set({1,n_,1,m_,1,l_}, 14,  1.0);
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


//-----------------------------------------------------------------------------
void Atmosphere::assemble()
{
	// clear old CRS matrix
	beg_.clear();
	ico_.clear();
	jco_.clear();
	
	int i2,j2,k2; // will contain neighbouring grid pointes given by shift()
	int row;
	int elm_ctr = 1;
	double value;
	for (int k = 1; k <= l_; ++k)
		for (int j = 1; j <= m_; ++j)
			for (int i = 1; i <= n_; ++i)
				for (int A = 1; A <= nun_; ++A)
				{
					// Filling new row:
					//  find the row corresponding to A at (i,j,k):
					row = find_row(i, j, k, A);
					//  put element counter in beg:					
					beg_.push_back(elm_ctr);
					for (int loc = 1; loc <= np_; ++loc)
					{
						// find index of neighbouring point loc
						shift(i,j,k,i2,j2,k2,loc);						
						for (int B = 1; B <= nun_; ++B)
						{
							value = Al_->get(i,j,k,loc,A,B);
							if (std::abs(value) > 0)
							{
								// store value								
								ico_.push_back(value);
								jco_.push_back(find_row(i2,j2,k2,B));
								// increment the element counter
								++elm_ctr;
							}
						}
					}
				}
	
	// final element of beg
	beg_.push_back(elm_ctr);
	
	// create dense A and its LU for solving with dgetrs()
	buildDenseA();
}

//-----------------------------------------------------------------------------
// Declaring a few LAPACK functions needed around this point
//-----------------------------------------------------------------------------
// Create LU
extern "C" void dgetrf_(int* M, int *N, double *A,
						int *LDA, int *IPIV, int *INFO);

// Solve LU 
extern "C" void dgetrs_(char* TRANS, int *N, int *NRHS, double *A,
						int *LDA, int *IPIV, double *B,
						int *LDB, int *INFO);

// Make LU and solve
extern "C" void dgesv_(int *N, int *NRHS, double *A,
					   int *LDA, int *IPIV, double *B,
					   int *LDB, int *INFO);

//-----------------------------------------------------------------------------
void Atmosphere::solve(std::shared_ptr<Vector> rhs)
{
	char trans = 'N';
	
	int dim    = n_*m_*l_;
	int nrhs   = 1;
	int lda    = dim;
	int ldb    = dim;
	int info;

	if (rhs == nullptr)
	{
		(*sol_) = std::vector<double>((*rhs_));
		std::cout << "solve... using own rhs" << std::endl;
	}
	else
	{
		(*sol_) = std::vector<double>(*(rhs->getStdVector()));
		std::cout << "solve... using ext rhs" << std::endl;
	}

	dgetrs_(&trans, &dim, &nrhs, &denseA_[0], &lda, ipiv_,
			&(*sol_)[0], &ldb, &info);
	
	std::cout << "solve info: " << info << std::endl;
}

//-----------------------------------------------------------------------------
void Atmosphere::buildDenseA()
{
	//------------------------------------------------
	//--> weird stuff: in the future stick to sparse plx
	// Create dense matrix:
	denseA_.assign(pow(n_*m_*l_, 2), 0.0);
	std::vector<double> ivals(jco_.size(), 0.0);
	int rw  = 1; // row
	int idx = 1;
	while (rw  <= m_*n_*l_)
	{
		for (int k = beg_[rw-1]; k != beg_[rw]; ++k)
		{
			ivals[idx-1] = rw ;
			++idx;
		}
		++rw ;
	}
	int i,j;
	double val;
	for (int cntr = 0; cntr != ico_.size(); ++cntr)
	{
		val = ico_[cntr];
		i   = ivals[cntr];
		j   = jco_[cntr];
		idx = i + (j-1)*m_*n_*l_ - 1;
		denseA_[idx] = val;
	}

	// create LU factorisation, store it in denseA_
	int dim  = n_*m_*l_;
	int lda  = dim;
	int info;

	std::cout << " Creating LU " << std::endl;
	dgetrf_(&dim, &dim, &denseA_[0], &lda, ipiv_, &info);
}


//-----------------------------------------------------------------------------
void Atmosphere::boundaries()
{
	//! size of stencil/neighbourhood:
	//! +----------++-------++----------+
	//! | 12 15 18 || 3 6 9 || 21 24 27 |
	//! | 11 14 17 || 2 5 8 || 20 23 26 |
	//! | 10 13 16 || 1 4 7 || 19 22 25 |
	//! |  below   || center||  above   |
	//! +----------++-------++----------+
	
	int west, east, north, south, bottom, top;
	for (int i = 1; i <= n_; ++i)
		for (int j = 1; j <= m_; ++j)
			for (int k = 1; k <= l_; ++k)
			{
				west   = i-1;
				east   = i+1;
				north  = j+1;
				south  = j-1;
				bottom = k-1;
				top    = k+1;

				// western boundary
				if (west == 0)
				{
					Al_->set(i,j,k,5,TT_,TT_, 
							 Al_->get(i,j,k,5,TT_,TT_) +
							 Al_->get(i,j,k,2,TT_,TT_));
					Al_->set(i,j,k,2,TT_,TT_, 0.0);
				}
				
				// eastern boundary
				if (east == n_+1)
				{
					Al_->set(i,j,k,5,TT_,TT_, 
							 Al_->get(i,j,k,5,TT_,TT_) +
							 Al_->get(i,j,k,8,TT_,TT_));
					Al_->set(i,j,k,8,TT_,TT_, 0.0);
				}
				
				// northern boundary
				if (north == m_+1)
				{
					Al_->set(i,j,k,5,TT_,TT_, 
							 Al_->get(i,j,k,5,TT_,TT_) +
							 Al_->get(i,j,k,6,TT_,TT_));
					Al_->set(i,j,k,6,TT_,TT_, 0.0);
				}

				// southern boundary
				if (south == 0)
				{
					Al_->set(i,j,k,5,TT_,TT_, 
							 Al_->get(i,j,k,5,TT_,TT_) +
							 Al_->get(i,j,k,4,TT_,TT_));
					Al_->set(i,j,k,4,TT_,TT_, 0.0);
				}
			}
}

//-----------------------------------------------------------------------------
void Atmosphere::writeAll()
{
	std::cout << "Writing everything to output files..." << std::endl;
	//--> this calls for some abstraction
	// Write solution
	if (!sol_->empty())
	{
		std::ofstream atmos_sol;
		atmos_sol.open("atmos_sol.txt");
		for (auto &i : *sol_)
			atmos_sol << std::setprecision(12) << i << '\n';
		atmos_sol.close();
	}
	else
		std::cout << " solution vector is empty" << std::endl;

	// Write rhs
	if (!rhs_->empty())
	{
		std::ofstream atmos_rhs;
		atmos_rhs.open("atmos_rhs.txt");
		for (auto &i : *rhs_)
			atmos_rhs << std::setprecision(12) << i << '\n';
		atmos_rhs.close();
	}
	else
		std::cout << " rhs vector is empty" << std::endl;

	// Write ico
	if (!ico_.empty())
	{
		std::ofstream atmos_ico;
		atmos_ico.open("atmos_ico.txt");
		for (auto &i : ico_)
			atmos_ico << std::setprecision(12) << i << '\n';
		atmos_ico.close();
	}
	else
		std::cout << " ico vector is empty" << std::endl;

	// Write jco
	if (!jco_.empty())
	{
		std::ofstream atmos_jco;
		atmos_jco.open("atmos_jco.txt");
		for (auto &i : jco_)
			atmos_jco << std::setprecision(12) << i << '\n';
		atmos_jco.close();
	}
	else
		std::cout << " jco vector is empty" << std::endl;

	// Write beg
	if (!beg_.empty())
	{
		std::ofstream atmos_beg;
		atmos_beg.open("atmos_beg.txt");
		for (auto &i : beg_)
			atmos_beg << std::setprecision(12) << i << '\n';
		atmos_beg.close();
	}
	else
		std::cout << " beg vector is empty" << std::endl;

	// Write frc
	if (!frc_.empty())
	{
		std::ofstream atmos_frc;
		atmos_frc.open("atmos_frc.txt");
		for (auto &i : frc_)
			atmos_frc << std::setprecision(12) << i << '\n';
		atmos_frc.close();
	}
	else
		std::cout << " frc vector is empty" << std::endl;

	// Write oceanTemp
	if (!oceanTemp_.empty())
	{
		std::ofstream atmos_oceanTemp;
		atmos_oceanTemp.open("atmos_oceanTemp.txt");
		for (auto &i : oceanTemp_)
			atmos_oceanTemp << std::setprecision(12) << i << '\n';
		atmos_oceanTemp.close();
	}
	else
		std::cout << " oceanTemp vector is empty" << std::endl;

	// Write state
	if (!state_->empty())
	{
		std::ofstream atmos_state;
		atmos_state.open("atmos_state.txt");
		for (auto &i : *state_)
			atmos_state << std::setprecision(12) << i << '\n';
		atmos_state.close();
	}
	else
		std::cout << " state vector is empty" << std::endl;
}

//-----------------------------------------------------------------------------
int Atmosphere::find_row(int i, int j, int k, int XX)
{
	// 1-based	
	return nun_ * ((k-1)*n_*m_ + n_*(j-1) + (i-1)) + XX;
}

//-----------------------------------------------------------------------------
void Atmosphere::shift(int i, int j, int k,
					   int &i2, int &j2, int &k2, int loc)
{
	if (loc < 10)
	{
		k2 = k;
        //   +-------+    +---------+
		//   | 3 6 9 |	  | 1  1  1 |
		//   | 2 5 8 | -> | 0  0  0 |
		//   | 1 4 7 |	  |-1 -1 -1 |
		//   +-------+	  +---------+
		j2 = j + ((loc + 2) % 3) - 1;
		//   +-------+    +---------+
		//   | 3 6 9 |	  |-1  0  1 |
		//   | 2 5 8 | -> |-1  0  1 |
		//   | 1 4 7 |	  |-1  0  1 |
		//   +-------+	  +---------+
		i2 = i + (loc - 1) / 3 - 1;
	}
	else if (loc < 19)
	{
		k2 = k - 1;
        //   +----------+     +---------+
		//   | 12 15 18 |	  | 1  1  1 |
		//   | 11 14 17 | ->  | 0  0  0 |
		//   | 10 13 16 |	  |-1 -1 -1 |
	    //   +----------+	  +---------+
		j2 = j + ((loc + 2) % 3) - 1;
		//   +----------+     +---------+
		//   | 12 15 18 |	  |-1  0  1 |
		//   | 11 14 17 | ->  |-1  0  1 |
		//   | 10 13 16 |	  |-1  0  1 |
	    //   +----------+	  +---------+
		i2 = i + (loc - 10) / 3 - 1;
	}
	else
	{
		k2 = k + 1;
		//   +----------+     +---------+
		//   | 21 24 27 |	  | 1  1  1 |
		//   | 20 23 26 | ->  | 0  0  0 |
		//   | 19 22 25 | 	  |-1 -1 -1 |
		//   +----------+	  +---------+
		j2 = j + ((loc + 2) % 3) - 1;
		//   +----------+     +---------+
		//   | 21 24 27 |	  |-1  0  1 |
		//   | 20 23 26 | ->  |-1  0  1 |
		//   | 19 22 25 | 	  |-1  0  1 |
		//   +----------+	  +---------+
		i2 = i + (loc - 19) / 3 - 1;
	}
}

//-----------------------------------------------------------------------------
void Atmosphere::test()
{
	std::cout << "Atmosphere: tests" << std::endl;
	
	std::cout << "  dense A size: " << denseA_.size() << std::endl;
	std::cout << "  outputting dense A in atmos_denseA.txt: " << std::endl;

	std::ofstream atmos_denseA;
	atmos_denseA.open("atmos_denseA.txt");
	for (auto &i : denseA_)
		atmos_denseA << std::setprecision(12) << i << '\n';
	atmos_denseA.close();
	
	std::cout << "            testing shift..." << std::endl;

	int i = 2;
	int j = 3;
	int k = 1;

	int i2, j2, k2;
	
	int loc = 21;
	
	std::cout << "  +----------++-------++----------+ " << '\n'
			  << "  | 12 15 18 || 3 6 9 || 21 24 27 | " << '\n'
			  << "  | 11 14 17 || 2 5 8 || 20 23 26 | " << '\n'
			  << "  | 10 13 16 || 1 4 7 || 19 22 25 | " << '\n'
			  << "  |  below   || center||  above   | " << '\n'
			  << "  +----------++-------++----------+ " << std::endl;
	
	std::cout << "(i,j,k,loc) = "
			  << "(" << i << ", " << j << ", "
			  << k << ", " << loc << ")" << std::endl;
	
	shift(i,j,k,i2,j2,k2,loc);
	std::cout << "(i,j,k,loc) = "
			  << "(" << i << ", " << j << ", "
			  << k << ", " << loc << ")" << std::endl;

	std::cout << "shift gives: " 
			  << "(" << i2 << ", " << j2 << ", " << k2 << ")" << std::endl;
	
	std::cout << "  test matvec(row), result is in atmos_test.txt" << std::endl;
	
	computeJacobian();
	
	std::vector<double> test;
	int row;
	for (int j = 1; j <= m_; ++j)
	{
		for (int i = 1; i <= n_; ++i)
		{
			row = find_row(i,j,l_,TT_);
			test.push_back(matvec(row));
		}
	}
	std::ofstream atmos_test;
	atmos_test.open("atmos_test.txt");
	for (auto &i : test)
		atmos_test << std::setprecision(12) << i << '\n';
	atmos_test.close();
}

//-----------------------------------------------------------------------------
std::shared_ptr<Vector> Atmosphere::getVector
(char mode, std::shared_ptr<std::vector<double> > vec)
{
	if (mode == 'C') // copy
	{
		std::shared_ptr<std::vector<double> > copy =
			std::make_shared<std::vector<double> >(*vec); 
		std::shared_ptr<Vector> ptr = std::make_shared<Vector>(copy);
		return ptr;
	}
	else if (mode == 'V') // view
	{
		std::shared_ptr<Vector> ptr = std::make_shared<Vector>(vec);
		return ptr;
	}
	else
	{
		std::cout << "WARNING (Atmosphere::getVector): Invalid mode"
				  << __FILE__ <<  __LINE__ << std::endl;
		return nullptr;
	}	
}

//-----------------------------------------------------------------------------
std::shared_ptr<Vector> Atmosphere::getSolution(char mode)
{
	return getVector(mode, sol_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<Vector> Atmosphere::getState(char mode)
{
	return getVector(mode, state_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<Vector> Atmosphere::getRHS(char mode)
{
	return getVector(mode, rhs_);
}

//-----------------------------------------------------------------------------

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
{}

//-----------------------------------------------------------------------------
DependencyGrid::~DependencyGrid()
{}

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
{}

//-----------------------------------------------------------------------------
Atom::~Atom()
{}

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
