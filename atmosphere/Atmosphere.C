#include "Atmosphere.H"
#include "AtmosphereDefinitions.H"
#include "GlobalDefinitions.H"
#include "SuperVector.H"

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm> // std::fill in assemble

#include <hdf5.h>

//---------------------------------------------------------------------------
// Constructor, specify horizontal grid dimensions
Atmosphere::Atmosphere(ParameterList params)
	:
	params_          (params),
	
	n_               (params->get("Global Grid-Size n", 16)),
	m_               (params->get("Global Grid-Size m", 16)),
	l_               (params->get("Global Grid-Size l", 1)),
	periodic_        (params->get("Periodic", false)),

	solvingScheme_   (params->get("Solving scheme", 'B')),
	preconditioner_  (params->get("Preconditioner", 'J')),
	gmresSolver_     (*this),
	gmresInitialized_(false),
	
	rhoa_            (params->get("atmospheric density",1.25)),
	hdima_           (params->get("heat capacity",8400.)),
	cpa_             (params->get("heat capacity",1000.)),
	d0_              (params->get("constant eddy diffusivity",3.1e+06)),
	arad_            (params->get("radiative flux param A",216.0)),
	brad_            (params->get("radiative flux param B",1.5)),
	sun0_            (params->get("solar constant",1360.)),
	c0_              (params->get("atmospheric absorption coefficient",0.43)),
	ce_              (params->get("exchange coefficient ce",1.3e-03)),
	ch_              (params->get("exchange coefficient ch",0.94 * ce_)),
	uw_              (params->get("mean atmospheric surface wind speed",8.5)),
	t0_              (params->get("reference temperature",15.0)),
	udim_            (params->get("horizontal velocity of the ocean",0.1e+00)),
	r0dim_           (params->get("radius of the earth",6.37e+06)),
	loadState_       (params->get("Load state", false)),
	saveState_       (params->get("Save state", false)),
	inputFile_       (params->get("Input file", "atmos.h5")),
	outputFile_      (params->get("Output file", "atmos.h5"))
{
	INFO("Atmosphere: constructor...");

	// Dimension of atmosphere
	dim_ = m_ * n_ * l_;

	// Number of super and sub-diagonal bands in banded matrix
	ksub_ = std::max(n_, m_);
	ksup_ = std::max(n_, m_);

	// Leading dimension of banded matrix
	ldimA_  = 2 * ksub_ + 1 + ksup_;
	
	// Continuation parameters
	ampl_    = 0.0        ; //! amplitude of forcing
	amplEnd_ = 1.0        ; //!
	
	// Filling the coefficients
	muoa_ =  rhoa_ * ch_ * cpa_ * uw_;
	amua_ = (arad_ + brad_ * t0_) / muoa_;
	bmua_ =  brad_ / muoa_;
	Ai_   =  rhoa_ * hdima_ * cpa_ * udim_ / (r0dim_ * muoa_);
	Ad_   =  rhoa_ * hdima_ * cpa_ * d0_ / (muoa_ * r0dim_ * r0dim_);
	As_   =  sun0_ * (1 - c0_) / (4 * muoa_);		
	
	np_  = ATMOS_NP_;   // all neighbouring points including the center
	nun_ = ATMOS_NUN_;  // only temperature ATMOS_TT_
	
	// Initialize state, rhs and solution of linear solve with zeros
	rhs_   = std::make_shared<std::vector<double> >(dim_, 0.0);
	sol_   = std::make_shared<std::vector<double> >(dim_, 0.0);
	state_ = std::make_shared<std::vector<double> >(dim_, 0.0);

	// Initialize ocean temperature
	oceanTemp_ = std::vector<double>(n_ * m_, 0.0);
	
	// Initialize forcing with zeros
	frc_ = std::vector<double>(n_ * m_ * l_, 0.0);

	// Initialize dense matrix:
	if (solvingScheme_ == 'D')
		denseA_ = std::vector<double>(n_ * m_ * l_ * n_ * m_ * l_, 0.0);

	// Initialize banded storage
	bandedA_     = std::vector<double>(ldimA_  * dim_, 0.0);
	buildLU_     = true;
	
	// Create pivot array for use in lapack
	ipiv_ = new int[n_*m_*l_+1];

	// Construct dependency grid:
	Al_ = std::make_shared<DependencyGrid>(n_, m_, l_, np_, nun_);

	xmin_ = params->get("Global Bound xmin", 286.0) * PI_ / 180.0;
	xmax_ = params->get("Global Bound xmax", 350.0) * PI_ / 180.0;
	ymin_ = params->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
	ymax_ = params->get("Global Bound ymax", 74.0)  * PI_ / 180.0;
	
	// Set the grid increments
	dx_ = (xmax_ - xmin_) / n_;
	dy_ = (ymax_ - ymin_) / m_;
	
	// Fill x
	xu_.reserve(n_+1);
	xc_.reserve(n_+1);
	for (int i = 0; i != n_+1; ++i)
	{
		xu_.push_back(xmin_ + i * dx_);
		xc_.push_back(xmin_ + (i - 0.5) * dx_);
	}
	
	// Fill y and latitude-based arrays
	yv_.reserve(m_+1);
	yc_.reserve(m_+1);
	albe_.reserve(m_+1);
	datc_.reserve(m_+1);
	datv_.reserve(m_+1);
	suna_.reserve(m_+1);
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

	// If specified we load a pre-existing state and parameter (x,l)
	if (loadState_)
		loadStateFromFile(inputFile_);

	INFO("Atmosphere: constructor... done");
}

//-----------------------------------------------------------------------------
// Destructor
Atmosphere::~Atmosphere()
{
	delete[] ipiv_;
}

//-----------------------------------------------------------------------------
void Atmosphere::initializeGMRES()
{
	gmresSolver_.setParameters(params_);
	gmresInitialized_ = true;
}

//-----------------------------------------------------------------------------
void Atmosphere::GMRESSolve(std::shared_ptr<SuperVector> rhs)
{} // NOT IMPLEMENTED

//-----------------------------------------------------------------------------
void Atmosphere::GMRESSolve(std::shared_ptr<SuperVector> rhs,
							std::shared_ptr<SuperVector> out)
{
	if (!gmresInitialized_)
		initializeGMRES();

	int verbosity = params_->get("GMRES verbosity", 0);
	
	TIMER_START("Atmosphere: solve...");

	gmresSolver_.setSolution (out);
	gmresSolver_.setRHS      (rhs);
	
	if (verbosity > 4) INFO("Atmosphere: GMRES solve...");
	gmresSolver_.solve();
	if (verbosity > 4) INFO("Atmosphere: GMRES solve... done");

	int iters    = gmresSolver_.getNumIters();
	double nrm   = gmresSolver_.residual();

	if (verbosity > 4) INFO("Atmosphere GMRES, i = " << iters << " residual = " << nrm);
	
	TRACK_ITERATIONS("Atmosphere GMRES iterations...", iters);
	TIMER_STOP("Atmosphere: solve...");
}


//-----------------------------------------------------------------------------
void Atmosphere::idealizedOcean()
{
	// put idealized values in oceanTemp
	double value;
	int row;
	for (int i = 1; i <= n_; ++i)
		for (int j = 1; j <= m_; ++j)
		{
			value = ampl_*cos(PI_*(yc_[j]-ymin_)/(ymax_-ymin_));
			row   = find_row(i,j,l_,ATMOS_TT_)-1;
			oceanTemp_[row] = value;
		}
}

//-----------------------------------------------------------------------------
void Atmosphere::idealizedState()
{
	// put idealized values in atmosphere
	double value;
	int row;
	for (int i = 1; i <= n_; ++i)
		for (int j = 1; j <= m_; ++j)
		{
			value = ampl_ * cos(PI_*(yc_[j]-ymin_)/(ymax_-ymin_));
			row   = find_row(i,j,l_,ATMOS_TT_)-1;
			(*state_)[row] = value;
		}
}

//-----------------------------------------------------------------------------
void Atmosphere::zeroState()
{
	// Set state to zero
	int dim = n_ * m_ * l_;
	*state_ = std::vector<double>(dim, 0.0);
}

//-----------------------------------------------------------------------------
void Atmosphere::zeroOcean()
{
	// Set ocean to zero
	oceanTemp_ = std::vector<double>(n_ * m_, 0.0);
}

//-----------------------------------------------------------------------------
void Atmosphere::setOceanTemperature(std::vector<double> const &sst)
{
	// Set ocean temperature (copy)
	oceanTemp_ = sst;
}

//-----------------------------------------------------------------------------
void Atmosphere::computeJacobian()
{
	TIMER_START("Atmosphere: compute Jacobian...");
	
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
	buildLU_ = true;

	//jac_ = getJacobian();
	
	TIMER_STOP("Atmosphere: compute Jacobian...");
}

//-----------------------------------------------------------------------------
std::shared_ptr<Atmosphere::CRSMat> Atmosphere::getJacobian()
{
	std::shared_ptr<CRSMat> jacMap;
	jacMap = std::make_shared<CRSMat>();
	(*jacMap)["ico"] = ico_;	
	(*jacMap)["jco"] = jco_;
	(*jacMap)["beg"] = beg_;
	return jacMap;		
}

//-----------------------------------------------------------------------------
void Atmosphere::computeRHS()
{
	TIMER_START("Atmosphere: compute RHS...");

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
			row   = find_row(i, j, l_, ATMOS_TT_);
			value = matvec(row) + frc_[row-1];
			(*rhs_)[row-1] = value;
		}
	
	TIMER_STOP("Atmosphere: compute RHS...");
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
	double value;
	int row;
	for (int j = 1; j <= m_; ++j)
		for (int i = 1; i <= n_; ++i)
		{
			row = find_row(i, j, l_, ATMOS_TT_);
			value = oceanTemp_[row-1] + ampl_ * (suna_[j] - amua_);
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
		atom.set({1,n_,1,m_,1,l_}, 5, -1.0);
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
	// Create CRS matrix storage and/or padded banded storage 

	// clear old CRS matrix
	beg_.clear();
	ico_.clear();
	jco_.clear();

	// Clear banded storage (just to be sure)
	std::fill(bandedA_.begin(), bandedA_.end(), 0.0);
	
	int i2,j2,k2; // will contain neighbouring grid pointes given by shift()
	int row;
	int rowb; // for banded storage
	int col;
	int colb; // for banded storage
	int idx;
	int kdiag = ksub_ + ksup_ + 1; // for banded storage
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
								// CRS --------------------------------------
								ico_.push_back(value);
								col = find_row(i2,j2,k2,B);
								jco_.push_back(col);

                                // increment the element counter
								++elm_ctr;

								// BND --------------------------------------
								// get row index for banded storage
								rowb = row - col + kdiag;
	
								// put matrix values in column major fashion
								// in the array
                                //  > go from 1 to 0-based
								rowb = rowb;
								colb = col;
								idx  = rowb + (colb - 1) * ldimA_ - 1;
								bandedA_[idx] = value;								
							}
						}
					}
				}
	
	// final element of beg
	beg_.push_back(elm_ctr);

	// create dense A and its LU for solving with dgetrs()
	if (solvingScheme_ == 'D')
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

// Solve banded system stored in bandedA_
extern "C" void dgbsv_(int *N, int *KL, int *KU, int *NRHS, double *AB,
					   int *LDAB, int *IPIV, double *B,
					   int *LDB, int *INFO);

// Solve banded system stored in bandedA_ EXPERT MODE
extern "C" void dgbsvx_(char *FACT, char *TRANS,
						int *N, int *KL, int *KU, int *NRHS, double *AB,
						int *LDAB, double *AFB, int *LDAFB, int *IPIV,
						char *EQUED, double *R, double *C,
						double *B, int *LDB,
						double *X, int *LDX,
						double *RCOND, double *FERR, double *BERR,
						double *WORK, int *IWORK, int *INFO);

// Create LU factorization of banded system
extern "C" void dgbtrf_(int *M, int *N, int *KL, int *KU, double *AB,
						int *LDAB, int *IPIV, int *INFO);

// Solve system using LU factorization given by dgbtrf
extern "C" void dgbtrs_(char *TRANS, int *N, int *KL, int *KU, int *NRHS,
						double *AB, int *LDAB, int *IPIV, double *B, int *LDB,
						int *INFO);

//-----------------------------------------------------------------------------
void Atmosphere::solve(std::shared_ptr<SuperVector> rhs)
{
	TIMER_START("Atmosphere: solve...");
	
	char trans  = 'N';
	int dim     = n_*m_*l_;
	int nrhs    = 1;
	int lda     = dim;
	int ldb     = dim;
	int info;

	if (rhs == nullptr)
		(*sol_) = std::vector<double>((*rhs_));
	else
		(*sol_) = std::vector<double>(*(rhs->getAtmosVector()));
	
	if (solvingScheme_ == 'D')
	{
		dgetrs_(&trans, &dim, &nrhs, &denseA_[0], &lda, ipiv_,
				&(*sol_)[0], &ldb, &info);
	}
	else if (solvingScheme_ == 'B')
	{
		// we use a copy to make sure bandedA_ does not get corrupted
		std::vector<double> bandedAcopy(bandedA_);
		dgbsv_(&dim, &ksub_, &ksup_, &nrhs, &bandedAcopy[0],
			   &ldimA_, ipiv_, &(*sol_)[0], &ldb, &info);
	}
	else if (solvingScheme_ == 'X')
	{// TODO		
	}
	else
		ERROR("Invalid solving scheme...", __FILE__, __LINE__);

	double residual = computeResidual(rhs);
	INFO("Atmosphere: solve... residual=" << residual);
	TIMER_STOP("Atmosphere: solve...");
}
// FACTORIZE!!!!!!!!
//-----------------------------------------------------------------------------
void Atmosphere::solve(SuperVector const &rhs, SuperVector &out)
{
	TIMER_START("Atmosphere: solve...");
	
	char trans  = 'N';
	
	int dim     = n_*m_*l_;
	int nrhs    = 1;
	int lda     = dim;
	int ldb     = dim;
		
	int info;
	
	if (solvingScheme_ == 'D')
	{
		out.assign(rhs.getAtmosVector());
		std::shared_ptr<std::vector<double> > sol = out.getAtmosVector();	
		dgetrs_(&trans, &dim, &nrhs, &denseA_[0], &lda, ipiv_,
				&(*sol)[0], &ldb, &info);
	}
	else if (solvingScheme_ == 'B')
	{ 
		if (buildLU_)
		{
			dgbtrf_(&dim, &dim, &ksub_, &ksup_, &bandedA_[0],
					&ldimA_, ipiv_, &info);
			buildLU_ = false; // until next request
		}
		
		out.assign(rhs.getAtmosVector());
		std::shared_ptr<std::vector<double> > sol = out.getAtmosVector();

		dgbtrs_(&trans, &dim, &ksub_, &ksup_, &nrhs,
				&bandedA_[0], &ldimA_, ipiv_,  &(*sol)[0], &ldb, &info);
	}
	else if (solvingScheme_ == 'G')
	{
		std::shared_ptr<SuperVector> b =
			std::make_shared<SuperVector>(rhs.getAtmosVector());
		std::shared_ptr<SuperVector> v =
			std::make_shared<SuperVector>(out.getAtmosVector());
		GMRESSolve(b, v);
		out.assign(v->getAtmosVector());
	}
	
	TIMER_STOP("Atmosphere: solve...");
}

//----------------------------------------------------------------------------
double Atmosphere::computeResidual(std::shared_ptr<SuperVector> rhs)
{
	std::shared_ptr<SuperVector> x = getSolution('C'); // obtain solution
	std::shared_ptr<SuperVector> b =
		std::make_shared<SuperVector>(rhs->getAtmosVector());
	x->linearTransformation(getJacobian());            // calculate J*x
	x->update(-1, *b, 1);                              // calculate Jx - b
	return x->norm();                                  // return ||b-Jx||
}

//----------------------------------------------------------------------------
double Atmosphere::computeResidual(SuperVector const &rhs,
								   SuperVector const &x)
{
	SuperVector r;
	r.assign(x.getAtmosVector());
	
	SuperVector b;
	b.assign(rhs.getAtmosVector());
	
	applyMatrix(x, r);
	r.update(1, b, -1);
	
	return r.norm();     // return ||b-Jx||
}

//-----------------------------------------------------------------------------
void Atmosphere::buildDenseA()
{
	//------------------------------------------------
	// This is horrible, only use when desparate (solvingScheme == 'D')
	if (solvingScheme_ != 'D')
		WARNING("this is not supposed to happen", __FILE__, __LINE__);

	// create dense matrix
	// reset denseA array
	denseA_.assign(n_ * m_ * l_ * n_ * m_ * l_, 0.0);
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
	for (size_t cntr = 0; cntr != ico_.size(); ++cntr)
	{
		val = ico_[cntr];
		i   = ivals[cntr];
		j   = jco_[cntr];
		idx = i + (j-1)*m_*n_*l_ - 1;
		denseA_[idx] = val;
	}

	write(denseA_, "atmos_denseA.txt");
	// create LU factorisation, store it in denseA_
	int dim  = n_*m_*l_;
	int lda  = dim;
	int info;
	dgetrf_(&dim, &dim, &denseA_[0], &lda, ipiv_, &info);
}


//-----------------------------------------------------------------------------
//! size of stencil/neighbourhood:
//! +----------++-------++----------+
//! | 12 15 18 || 3 6 9 || 21 24 27 |
//! | 11 14 17 || 2 5 8 || 20 23 26 |
//! | 10 13 16 || 1 4 7 || 19 22 25 |
//! |  below   || center||  above   |
//! +----------++-------++----------+
void Atmosphere::boundaries()
{
	
	int west, east, north, south;
	for (int i = 1; i <= n_; ++i)
		for (int j = 1; j <= m_; ++j)
			for (int k = 1; k <= l_; ++k)
			{
				west   = i-1;
				east   = i+1;
				north  = j+1;
				south  = j-1;

				// western boundary
				if (west == 0 && !periodic_)
				{
					Al_->set(i,j,k,5,ATMOS_TT_,ATMOS_TT_, 
							 Al_->get(i,j,k,5,ATMOS_TT_,ATMOS_TT_) +
							 Al_->get(i,j,k,2,ATMOS_TT_,ATMOS_TT_));
					Al_->set(i,j,k,2,ATMOS_TT_,ATMOS_TT_, 0.0);
				}
				
				// eastern boundary
				if (east == n_+1 && !periodic_)
				{
					Al_->set(i,j,k,5,ATMOS_TT_,ATMOS_TT_, 
							 Al_->get(i,j,k,5,ATMOS_TT_,ATMOS_TT_) +
							 Al_->get(i,j,k,8,ATMOS_TT_,ATMOS_TT_));
					Al_->set(i,j,k,8,ATMOS_TT_,ATMOS_TT_, 0.0);
				}
				
				// northern boundary
				if (north == m_+1)
				{
					Al_->set(i,j,k,5,ATMOS_TT_,ATMOS_TT_, 
							 Al_->get(i,j,k,5,ATMOS_TT_,ATMOS_TT_) +
							 Al_->get(i,j,k,6,ATMOS_TT_,ATMOS_TT_));
					Al_->set(i,j,k,6,ATMOS_TT_,ATMOS_TT_, 0.0);
				}

				// southern boundary
				if (south == 0)
				{
					Al_->set(i,j,k,5,ATMOS_TT_,ATMOS_TT_, 
							 Al_->get(i,j,k,5,ATMOS_TT_,ATMOS_TT_) +
							 Al_->get(i,j,k,4,ATMOS_TT_,ATMOS_TT_));
					Al_->set(i,j,k,4,ATMOS_TT_,ATMOS_TT_, 0.0);
				}
			}
}

//-----------------------------------------------------------------------------
void Atmosphere::write(std::vector<double> &vector, const std::string &filename)
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

//-----------------------------------------------------------------------------
void Atmosphere::writeAll()
{
	std::cout << "Writing everything to output files..." << std::endl;

	write(*sol_, "atmos_sol.txt"); 
	write(*rhs_, "atmos_rhs.txt"); 
	write(ico_, "atmos_ico.txt"); 
	write(jco_, "atmos_jco.txt"); 
	write(beg_, "atmos_beg.txt"); 
	write(bandedA_, "atmos_bandedA.txt"); 	
	write(frc_, "atmos_frc.txt");           
	write(oceanTemp_, "atmos_oceanTemp.txt");
	write(*state_, "atmos_state.txt");       
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
			row = find_row(i,j,l_,ATMOS_TT_);
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
std::shared_ptr<SuperVector> Atmosphere::getVector
(char mode, std::shared_ptr<std::vector<double> > vec)
{
	if (mode == 'C')      // copy
	{
		std::shared_ptr<std::vector<double> > copy =
			std::make_shared<std::vector<double> >(*vec); 
		std::shared_ptr<SuperVector> ptr = std::make_shared<SuperVector>(copy);
		return ptr;
	}
	else if (mode == 'V') // view
	{
		std::shared_ptr<SuperVector> ptr = std::make_shared<SuperVector>(vec);
		return ptr;
	}
	else
	{
		WARNING("invalid mode", __FILE__, __LINE__);
		return std::shared_ptr<SuperVector>();
	}	
}

//-----------------------------------------------------------------------------
std::shared_ptr<SuperVector> Atmosphere::getSolution(char mode)
{
	return getVector(mode, sol_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<SuperVector> Atmosphere::getState(char mode)
{
	return getVector(mode, state_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<SuperVector> Atmosphere::getRHS(char mode)
{
	return getVector(mode, rhs_);
}

//-----------------------------------------------------------------------------
std::shared_ptr<SuperVector> Atmosphere::applyMatrix(SuperVector const &v)
{
	std::shared_ptr<std::vector<double> > atmosVector = v.getAtmosVector();
	std::shared_ptr<std::vector<double> > result =
		std::make_shared<std::vector<double> > (atmosVector->size(),0.0);
	
	int first;
	int last;	
	
	// Perform matrix vector product
	// 1->0 based... horrible... 
	for (size_t row = 1; row <= atmosVector->size(); ++row)
	{
		first   = beg_[row-1];
		last    = beg_[row] - 1;
		for (int col = first; col <= last; ++col)
			(*result)[row-1] += ico_[col-1] * (*atmosVector)[jco_[col-1]-1];
	}
	return getVector('V', result);
}

//-----------------------------------------------------------------------------
void Atmosphere::applyMatrix(SuperVector const &v, SuperVector &out)
{
 	std::shared_ptr<std::vector<double> > atmosVector = v.getAtmosVector();
	std::shared_ptr<std::vector<double> > result = out.getAtmosVector();
	
	int first;
	int last;	
	
	// Perform matrix vector product
	// 1->0 based... horrible... 
	for (size_t row = 1; row <= atmosVector->size(); ++row)
	{
		first = beg_[row-1];
		last  = beg_[row] - 1;
		
		(*result)[row-1] = 0;
		for (int col = first; col <= last; ++col)
			(*result)[row-1] += ico_[col-1] * (*atmosVector)[jco_[col-1]-1];
	}
}

//-----------------------------------------------------------------------------
std::shared_ptr<SuperVector> Atmosphere::applyPrecon(SuperVector const &v)
{
	// Do nothing
	std::shared_ptr<SuperVector> result = std::make_shared<SuperVector>(v);
	return result;
}

//-----------------------------------------------------------------------------
void Atmosphere::applyPrecon(SuperVector const &v, SuperVector &out)
{
	if (preconditioner_ == 'J')
	{
		std::shared_ptr<std::vector<double> > atmosVector = v.getAtmosVector();
		std::shared_ptr<std::vector<double> > result = out.getAtmosVector();
		int first;
		int last;
		double scaling = 1;
		
		// Perform matrix vector product
		// 1->0 based... horrible... 
		for (size_t row = 1; row <= atmosVector->size(); ++row)
		{
			first = beg_[row-1];
			last  = beg_[row] - 1;
			
			(*result)[row-1] = 0;
			for (int col = first; col <= last; ++col)
				if (row == jco_[col-1])
					(*result)[row-1] = scaling
						* (*atmosVector)[jco_[col-1]-1]
						/ ico_[col-1] ;
		}
	}
	else if (preconditioner_ == 'D')
	{
		solve(v, out);
	}
	else
		out.assign(v.getAtmosVector());
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<double> > Atmosphere::getOceanBlock()
{
	// The contribution of the ocean in the atmosphere is a
	// diagonal of ones.
	return std::make_shared<std::vector<double> >(m_*n_, 1.0);
}

//-----------------------------------------------------------------------------
std::shared_ptr<std::vector<int> > Atmosphere::getAtmosRows()
{
	std::shared_ptr<std::vector<int> > rows =
		std::make_shared<std::vector<int> >();
	for (int j = 1; j <= m_; ++j)
		for (int i = 1; i <= n_; ++i)
		{
			// 1-based!!
			rows->push_back(find_row(i,j,l_,ATMOS_TT_));
		}	
	return rows;	
}

//-----------------------------------------------------------------------------
void Atmosphere::setLandMask(std::shared_ptr<std::vector<int> > landm)
{
	if ((int) landm->size() != (n_+2)*(m_+2))
		WARNING("landm->size() not ok:",  __FILE__, __LINE__);

	landm_ = landm;

	INFO("Printing land surface mask...");
	int ctr = 0;
	std::ostringstream string;
	for (auto &l: *landm)
	{
		ctr++;
		string << l;
		if (ctr % (n_+2) == 0)
		{
			INFO(string.str().c_str());
			string.str("");
			string.clear();
		}
	}
}

//-----------------------------------------------------------------------------
void Atmosphere::preProcess()
{
	// nothing to do here (yet)
}

//-----------------------------------------------------------------------------
void Atmosphere::postProcess()
{
	if (saveState_)
		saveStateToFile("atmos.h5");
	
	write(*state_, "atmos_state.txt");       
}

//-----------------------------------------------------------------------------
void Atmosphere::saveStateToFile(std::string const &filename)
{
	INFO("Writing to " << filename);
	
    hid_t       file_id, group_id, dataspace_id, dataset_id;
 	hsize_t     dim_state = state_->size();
	hsize_t     dim_par = 1;
	int         nstat = 9;
	herr_t      status[nstat];
	int         si = 0;
	
	// Create a new file 
	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create a group 
	group_id = H5Gcreate2(file_id, "/State", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create the data space and dataset
	dataspace_id = H5Screate_simple(1, &dim_state, NULL);
	dataset_id   = H5Dcreate2(file_id, "/State/Values", H5T_IEEE_F64LE, dataspace_id,
							  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write to the dataset
	status[si++] = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
							&(*state_)[0]);	

	// Close dataset, dataspace and group
	status[si++] = H5Dclose(dataset_id);
	status[si++] = H5Sclose(dataspace_id );
	status[si++] = H5Gclose(group_id);

	// Create a new group 
	group_id = H5Gcreate2(file_id, "/Continuation parameter", H5P_DEFAULT,
						  H5P_DEFAULT, H5P_DEFAULT);

	// Create data space and dataset
	dataspace_id = H5Screate_simple(1, &dim_par, NULL);
	dataset_id   = H5Dcreate2(file_id, "/Continuation parameter/Value", H5T_IEEE_F64LE,
							  dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	// Write
	status[si++] = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
							&ampl_);	
	
	// Close everything
	status[si++] = H5Dclose(dataset_id);
	status[si++] = H5Sclose(dataspace_id );
	status[si++] = H5Gclose(group_id);
	status[si++] = H5Fclose(file_id); 

	// Check for errors
	for (int i = 0; i != nstat; ++i)
		if (status[i] != 0)
			WARNING("Status[" << i << "] not ok", __FILE__, __LINE__);

}

//-----------------------------------------------------------------------------
void Atmosphere::loadStateFromFile(std::string const &filename)
{
	TIMER_START("Atmosphere::loadStateFromFile()...");
	INFO("Loading from " << filename);

	hid_t    file_id, dataset_id;
	int      nstat = 5;
	herr_t   status[nstat];
	int      si = 0;

	// Check whether file exists
	std::ifstream file(filename);
	if (!file)
	{
		WARNING("Can't open " << filename
				<< " continue with trivial state", __FILE__, __LINE__);
		return;
	}
	else file.close();
	
	// Open existing file
	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	
	// Open state dataset
	dataset_id = H5Dopen2(file_id, "/State/Values", H5P_DEFAULT);
	
	// Read state from dataset and close it
	status[si++] = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
						   &(*state_)[0]);
	status[si++] = H5Dclose(dataset_id);

	// Open continuation parameter dataset
	dataset_id = H5Dopen2(file_id, "/Continuation parameter/Value", H5P_DEFAULT);

	// Read from dataset and close it
	status[si++] = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
						   &ampl_);

	// Close dataset and file
	status[si++] = H5Dclose(dataset_id);
	status[si++] = H5Fclose(file_id); 
	
	// Check for errors
	for (int i = 0; i != nstat; ++i)
		if (status[i] != 0)
			WARNING("Status[" << i << "] not ok", __FILE__, __LINE__);

	TIMER_STOP("Atmosphere::loadStateFromFile()...");
}

//=============================================================================
// / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / //
//                                                                           //
//                                                                           //
// DependencyGrid implementation                                             //
//                                                                           //
//                                                                           //
// / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / //
//=============================================================================
DependencyGrid::DependencyGrid(int n, int m, int l, int np, int nun)
	:
	grid_(n, m, l, np, nun, nun),
	n_(n),
	m_(m),
	l_(l),
	np_(np),
	nun_(nun)
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
// / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / //
//                                                                           //
//                                                                           //
// Atom implementation                                                       //
//                                                                           //
//                                                                           //
// / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / //
//=============================================================================
Atom::Atom(int n, int m, int l, int np)
	:
	atom_(n, m, l, np),
	n_(n),
	m_(m),
	l_(l),
	np_(np)
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
