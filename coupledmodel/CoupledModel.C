#include "CoupledModel.H"
#include "Ocean.H"
#include "Atmosphere.H"
#include "SuperVector.H"

#include <vector>
#include <memory>
#include <functional>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

CoupledModel::CoupledModel(Teuchos::RCP<Ocean> ocean,
						   std::shared_ptr<Atmosphere> atmos,
						   Teuchos::RCP<Teuchos::ParameterList> params)
	:
	ocean_(ocean),
	atmos_(atmos),
	useExistingState_ (params->get("Use existing state", false)),
	solvingScheme_    (params->get("Solving scheme", 'G')),
	useScaling_       (params->get("Use scaling", false)),
	iterGS_           (params->get("Max GS iterations", 10)),
	toleranceGS_      (params->get("GS tolerance", 1e-1)),
	useHash_          (params->get("Use hashing", true)),
	syncHash_(-1),
	rhsHash_(-1),
	jacHash_(-1),
	syncCtr_(0),
	buildPrecEvery_   (params->get("Rebuild preconditioner stride", 1)),
	idrSolver_(*this),
	gmresSolver_(*this),
	idrInitialized_(false),
	gmresInitialized_(false),
	idrSolveCtr_(0)
{
	stateView_ =
		std::make_shared<SuperVector>(ocean_->getState('V')->getOceanVector(),
									  atmos_->getState('V')->getAtmosVector() );
	
	solView_ =
 		std::make_shared<SuperVector>(ocean_->getSolution('V')->getOceanVector(),
									  atmos_->getSolution('V')->getAtmosVector() );
	
	rhsView_ =
 		std::make_shared<SuperVector>(ocean_->getRHS('V')->getOceanVector(),
									  atmos_->getRHS('V')->getAtmosVector() );

	// Get the contribution of the atmosphere to the ocean in the Jacobian
	B_     = ocean_->getAtmosBlock();
	rowsB_ = ocean_->getSurfaceTRows();
	
	// Get the contribution of the ocean to the atmosphere in the Jacobian
	C_     = atmos_->getOceanBlock();
	
	// Output parameters
	INFO(*params);
	
	// Synchronize landmask
	synchronizeLandmask();
}

//------------------------------------------------------------------
void CoupledModel::synchronize()
{	
	// A synchronization is only necessary when the states have actually
	// changed, so we compute, compare and store a hash
	if (useHash_)
	{
		std::size_t hash = stateView_->hash();
		if (hash == syncHash_)
			return;
		else
			syncHash_ = hash;
	}

	TIMER_START("CoupledModel: synchronize...");

	syncCtr_++; // Keep track of synchronizations
	
	// Copy the atmosphere from the current combined state
	std::vector<double> atmos(*(stateView_->getAtmosVector()));

	// Copy the surface temperature from the ocean model
	std::vector<double> sst(*(ocean_->getSurfaceT()));
	
	// Set the atmosphere in the ocean
	ocean_->setAtmosphere(atmos);

	// Set the SST in the atmosphere
	atmos_->setOceanTemperature(sst);

	// Now that there is a new state we can also recompute the preconditioners
	if (syncCtr_ % buildPrecEvery_ == 0)
		ocean_->recomputePreconditioner();
	
	TIMER_STOP("CoupledModel: synchronize...");
}

void CoupledModel::synchronizeLandmask()
{
	// Get the landmask used in the ocean
	std::shared_ptr<std::vector<int> > landm = ocean_->getLandMask();

	// Get rid of every layer in the landmask except for the surface 
	int n = ocean_->getNdim();
	int m = ocean_->getMdim();
	int l = ocean_->getLdim();
	landm->erase(landm->begin(), landm->begin() + l*(m+2)*(n+2));
	landm->erase(landm->begin() + (m+2)*(n+2), landm->end());
	
	// Put the surface landmask in the atmosphere
	atmos_->setLandMask(landm);
}

//------------------------------------------------------------------
void CoupledModel::computeJacobian()
{
	// Check whether the combined (state,par) vector has changed,
	// if so continue, if not return
	if (useHash_)
	{
		std::size_t hash = getHash();
		if (hash == jacHash_)
			return;
		else
			jacHash_ = hash;
	}
	
	// Synchronize the states
	if (solvingScheme_ != 'D') { synchronize(); }
	
	ocean_->computeJacobian();	// Ocean
	atmos_->computeJacobian();	// Atmosphere
}

//------------------------------------------------------------------
void CoupledModel::computeRHS()
{
	// Check whether the combined (state,par) vector has changed,
	// if so continue, if not return
	if (useHash_)
	{
		std::size_t hash = getHash();
		if (hash == rhsHash_)
			return;
		else
			rhsHash_ = hash;
	}
	
	// Synchronize the states
	if (solvingScheme_ != 'D') { synchronize(); }
	
	ocean_->computeRHS();	// Ocean
	atmos_->computeRHS(); 	// Atmosphere
}

//====================================================================
void CoupledModel::initializeIDR()
{
	Teuchos::RCP<Teuchos::ParameterList> solverParams_ =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams_.ptr());
	
	idrSolver_.setParameters(solverParams_);
	idrSolver_.setSolution(getSolution('V'));
	idrSolver_.setRHS(getRHS('V'));
	clearSPFreq_ = solverParams_->get("Clear search space frequency", 2);
	idrInitialized_ = true;
}

//====================================================================
void CoupledModel::initializeGMRES()
{
	Teuchos::RCP<Teuchos::ParameterList> solverParams_ =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams_.ptr());
	
	gmresSolver_.setParameters(solverParams_);
	gmresInitialized_ = true;
}

//------------------------------------------------------------------
void CoupledModel::solve(std::shared_ptr<SuperVector> rhs)
{
	// Start solve
	TIMER_START("CoupledModel: solve...");
	if (useScaling_)
		ocean_->scaleProblem(Teuchos::rcp(rhs.get(), false));
	
	if (solvingScheme_ == 'D') // fully decoupled solve
	{
		ocean_->solve(Teuchos::rcp(rhs.get(), false));
		atmos_->solve(rhs);
	}
	else if (solvingScheme_ == 'B') // backward block GS solve
		blockGSSolve(rhs);
	else if (solvingScheme_ == 'I') // IDR on complete matrix
		IDRSolve(rhs);
	else if (solvingScheme_ == 'G') // GMRES on complete matrix
		GMRESSolve(rhs);
	else
		WARNING("(CoupledModel::Solve()) Invalid mode!",
				__FILE__, __LINE__);

	if (useScaling_)
		ocean_->unscaleProblem(Teuchos::rcp(rhs.get(), false));
	
	// Update the profile after a solve
	printProfile(profile);
	TIMER_STOP("CoupledModel: solve...");
}

//------------------------------------------------------------------
void CoupledModel::IDRSolve(std::shared_ptr<SuperVector> rhs)
{
	if (!idrInitialized_)
		initializeIDR();
		

	// Clear the search space every X calls:
	if (idrSolveCtr_ % clearSPFreq_ == 0)
		idrSolver_.clearSearchSpace();
		
	idrSolveCtr_++;
 	solView_->zero();
	idrSolver_.setSolution(getSolution('V'));
	idrSolver_.setRHS(rhs);

	INFO("CoupledModel: IDR solve...");
	idrSolver_.solve();
	INFO("CoupledModel: IDR solve... done");
		
	int iters    = idrSolver_.getNumIters();
	double nrm   = idrSolver_.explicitResNorm();
	double res   = computeResidual(rhs);

	INFO("CoupledModel IDR, i = " << iters << " exp res norm: " << nrm);
	INFO("CoupledModel residual = " << res);		
	TRACK_ITERATIONS("CoupledModel IDR iterations...", iters);
}

//------------------------------------------------------------------
void CoupledModel::GMRESSolve(std::shared_ptr<SuperVector> rhs)
{
	if (!gmresInitialized_)
		initializeGMRES();
		
	solView_->zero();
	gmresSolver_.setSolution(getSolution('V'));
	gmresSolver_.setRHS(rhs);

	INFO("CoupledModel: GMRES solve...");
	gmresSolver_.solve();
	INFO("CoupledModel: GMRES solve... done");
	
	int iters    = gmresSolver_.getNumIters();
	double nrm   = gmresSolver_.residual();
	double res   = computeResidual(rhs);

	INFO("CoupledModel GMRES, i = " << iters << " exp res norm: " << nrm);
	INFO("CoupledModel residual = " << res);		
	TRACK_ITERATIONS("CoupledModel GMRES iterations...", iters);
	TRACK_RESIDUAL("CoupledModel GMRES residual...", res);
}

//------------------------------------------------------------------
void CoupledModel::blockGSSolve(std::shared_ptr<SuperVector> rhs)
{
	// ***************************************************************
 	// Notation: J = [A,B;C,D], x = [x1;x2], b = [b1;b2]
	//           M = [A, 0; 0, D], E = [0, 0; -C, 0], F = [0, -B; 0, 0]
	//
	// Symmetric block GS: (M-F)*x^{k+1/2} = E*x^{k} + b
	//                     (M-E)*x^{k+1)   = F*x^{k+1/2) + b
	//
	// This leads to iteratively solving   D*x2 = -C*x1 + b2
	//                                     A*x1 = -B*x2 + b1
	// 
	// After the iteration we do a final solve with  D*x2 = -C*x1 + b2
	//  (because it's cheap)
	// ***************************************************************

	double residual;
	double old_residual = computeResidual(rhs);
	
    // Initialize solution [x1;x2] = 0
 	std::shared_ptr<SuperVector> x = getSolution('C', 'C');
	x->zero();

	// Clear saved Krylov space in ocean solver
	ocean_->clearSearchSpace();
	
	// Start iteration
	int i;
	for (i = 1; i <= iterGS_; ++i)
	{
		// Create -C*x1 + b2
		x->linearTransformation(*C_, *rowsB_, 'O', 'A');
		x->update(1, *rhs, -1);

		// Solve D*x2 = -C*x1 + b2
		atmos_->solve(x);

		// Retrieve solution
		x = getSolution('C','C');

		// Create -B*x2 + b1
		x->linearTransformation(*B_, *rowsB_, 'A', 'O');
		x->update(1, *rhs, -1);

		// Solve A*x1 = -B*x2 + b1
		ocean_->solve(Teuchos::rcp(x.get(), false));

		// Retrieve solution
		x = getSolution('C','C');
		
		// Calculate residual
		residual = computeResidual(rhs);
		INFO("CoupledModel: blockGS, i = " << i
			 << ", ||b-Jx||/||b|| = " << residual << ", tol = " << toleranceGS_);

		if (residual > old_residual)
			WARNING("INCREASING RESIDUAL!", __FILE__, __LINE__);
		
		if (residual < toleranceGS_)
			break;
	}

	// Do a final solve with D
	// Create -C*x1 + b2 and solve D*x2 = -C*x1 + b2
	x->linearTransformation(*C_, *rowsB_, 'O', 'A');
	x->update(1, *rhs, -1);
	atmos_->solve(x);

	TRACK_ITERATIONS("CoupledModel: blockGS iterations...", i);
	
	if (i == iterGS_)
		WARNING("GS tolerance not reached...", __FILE__, __LINE__);

}

//------------------------------------------------------------------
std::shared_ptr<SuperVector>
CoupledModel::applyMatrix(SuperVector const &v)
{
	TIMER_START("CoupledModel: apply matrix...");
	std::shared_ptr<SuperVector> x = std::make_shared<SuperVector>(v);
	SuperVector y(v);
	SuperVector z(v);
	
	x->linearTransformation(ocean_->getJacobian());   // A*x1
	x->linearTransformation(atmos_->getJacobian());   // D*x2

	y.linearTransformation(*B_, *rowsB_, 'A', 'O');  // B*x2
	y.zeroAtmos();
	z.linearTransformation(*C_, *rowsB_, 'O', 'A');  // C*x1
	z.zeroOcean();

	x->update(1, y, 1);                              // A*x1 + B*x2
	x->update(1, z, 1);                              // D*x2 + C*x1
	TIMER_STOP("CoupledModel: apply matrix...");
	return x;
}

//------------------------------------------------------------------
void CoupledModel::applyMatrix(SuperVector const &v, SuperVector &out, char mode)
{
	TIMER_START("CoupledModel: apply matrix...");

	out.zero();	// Initialize output	
	// Fill the ocean and atmos part of output
	ocean_->applyMatrix(v, out);  // A*x1
	atmos_->applyMatrix(v, out);  // D*x2

	if (mode == 'C')
	{
		// Make temporary copies of v to store linear transformations
		SuperVector y(v);
		SuperVector z(v);

		// Perform mappings
		y.linearTransformation(*B_, *rowsB_, 'A', 'O');  // B*x2
		z.linearTransformation(*C_, *rowsB_, 'O', 'A');  // C*x1

		// Just to be sure
		y.zeroAtmos();        
		z.zeroOcean();        

		out.update(1,y,1);  // A*x1 + B*x2
		out.update(1,z,1);  // D*x2 + C*x1
	}	
	TIMER_STOP("CoupledModel: apply matrix...");
}

//------------------------------------------------------------------
std::shared_ptr<SuperVector>
CoupledModel::applyPrecon(SuperVector const &v)
{
	TIMER_START("CoupledModel: apply preconditioner...");
	SuperVector tmp1 = *(ocean_->applyPrecon(v));
	SuperVector tmp2 = *(atmos_->applyPrecon(v));
	TIMER_STOP("CoupledModel: apply preconditioner...");	
	return std::make_shared<SuperVector>(tmp1.getOceanVector(),
										 tmp2.getAtmosVector());
}

//------------------------------------------------------------------
void CoupledModel::applyPrecon(SuperVector const &v, SuperVector &out)
{
	TIMER_START("CoupledModel: apply preconditioner...");	
	out.zero();	// Initialize output
	ocean_->applyPrecon(v, out);
	atmos_->applyPrecon(v, out);

	TIMER_STOP("CoupledModel: apply preconditioner...");	
}

//------------------------------------------------------------------
double CoupledModel::computeResidual(std::shared_ptr<SuperVector> rhs)
{
	TIMER_START("CoupledModel: compute residual...");
	
	SuperVector r(*solView_);
	applyMatrix(*solView_, r);
	
	r.update(1, *rhs, -1);                           //  b-Jx	
	double relResidual = r.norm() / rhs->norm();     // ||b-Jx||/||b||
	
	TIMER_STOP("CoupledModel: compute residual...");
	return relResidual;
}

//----------------------------------------------------------------------------
double CoupledModel::computeResidual(SuperVector const &rhs,
									 SuperVector const &x, char mode)
{
	SuperVector r(x);
	applyMatrix(x, r, mode);
	r.update(1, rhs, -1);	
	return r.norm('V');     // return ||b-Jx||
}

//------------------------------------------------------------------
std::shared_ptr<SuperVector> CoupledModel::getSolution(char mode)
{
	// obtain solution based on mode
	if (mode == 'V') // View
		return solView_;
	else if (mode == 'C') // Copy
	{
		return std::make_shared<SuperVector>(
			ocean_->getSolution('C')->getOceanVector(),
			atmos_->getSolution('C')->getAtmosVector() );
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return std::shared_ptr<SuperVector>();
	}
}

//------------------------------------------------------------------
std::shared_ptr<SuperVector> CoupledModel::getSolution(char mode1,
													   char mode2)
{
	return std::make_shared<SuperVector>(
		ocean_->getSolution(mode1)->getOceanVector(),
		atmos_->getSolution(mode2)->getAtmosVector() );
}

//------------------------------------------------------------------
std::shared_ptr<SuperVector> CoupledModel::getState(char mode)
{
	if (mode == 'V') // View
		return stateView_;
	else if (mode == 'C') // Copy
	{
		return std::make_shared<SuperVector>(
			ocean_->getState('C')->getOceanVector(),
			atmos_->getState('C')->getAtmosVector() );
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return std::shared_ptr<SuperVector>();
	}
}

//------------------------------------------------------------------
std::shared_ptr<SuperVector> CoupledModel::getRHS(char mode)
{
	if (mode == 'V') // View
		return rhsView_;
	else if (mode == 'C') // Copy
	{
		return std::make_shared<SuperVector>(
			ocean_->getRHS('C')->getOceanVector(),
			atmos_->getRHS('C')->getAtmosVector() );
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return std::shared_ptr<SuperVector>();
	}
}

//------------------------------------------------------------------
void CoupledModel::setState(std::shared_ptr<SuperVector> state)
{
	ocean_->setState(Teuchos::rcp(state.get(), false));
    atmos_->setState(state);
}

//------------------------------------------------------------------
void CoupledModel::setRHS(std::shared_ptr<SuperVector> rhs)
{
	ocean_->setRHS(Teuchos::rcp(rhs.get(), false));
    atmos_->setRHS(rhs);
}

//------------------------------------------------------------------
double CoupledModel::getPar()
{
	// The parameters should remain equal among the models
	// Different continuation parameters for different models
	// is not defined (yet).
	double par_ocean = ocean_->getPar();
	double par_atmos = atmos_->getPar();
	if (std::abs(par_ocean - par_atmos) > 1e-8)
	{
		WARNING("par_ocean != par_atmos" << " ocean: " << par_ocean
				<< " atmos: " << par_atmos
				<< " putting maximum in both models ", __FILE__, __LINE__);
		double max = std::max(par_ocean, par_atmos);
		ocean_->setPar(max);
		atmos_->setPar(max);
		return max;
	}
	return par_ocean;
}

//------------------------------------------------------------------
void CoupledModel::setPar(double value)
{
	ocean_->setPar(value);
	atmos_->setPar(value);
}

//------------------------------------------------------------------
double CoupledModel::getParDestination()
{
	// The parameters should remain equal among the models
	// Different continuation parameters for different models
	// is not defined (yet).
	double parDest_ocean = ocean_->getParDestination();
	double parDest_atmos = atmos_->getParDestination();
	if (std::abs(parDest_ocean - parDest_atmos) > 1e-8)
	{
		WARNING("parDest_ocean != parDest_atmos !!",
				__FILE__, __LINE__);
 		INFO("ocean: " << parDest_ocean);
		INFO("atmos: " << parDest_atmos);
		return -1;
	}
	return parDest_ocean;
}

//------------------------------------------------------------------
void CoupledModel::preProcess()
{
	ocean_->preProcess();
	atmos_->preProcess();
}

//------------------------------------------------------------------
void CoupledModel::postProcess()
{
	// If the solver is completely decoupled, this is the right
	// moment to synchronize
	if (solvingScheme_ == 'D')
		synchronize();

	// Let the models do their own post-processing
	ocean_->postProcess();
	atmos_->postProcess();

	// Print the contents of the profile
	printProfile(profile);
}

//------------------------------------------------------------------
// This is a copy of the code in main.C --> should be factorized
void CoupledModel::printProfile(ProfileType profile)	
{
	if (timerStack.empty() == false)
		WARNING("Unequal amount of TIMER_START and TIMER_STOP uses",
				__FILE__, __LINE__);
	
	std::ostringstream profilefile("profile_output");   // setting up a filename
	std::ofstream file(profilefile.str().c_str());      // setup output file

	// Set format flags
	file << std::left;

	// Define line format
#ifndef LINE
# define LINE(s1, s2, s3, s4, s5, s6, s7, s8, s9)						\
	{																	\
		int sp = 3;  int it = 5;  int id = 5;							\
		int db = 12; int st = 45;										\
		file << std::setw(id) << s1	<< std::setw(sp) << s2				\
			 << std::setw(st) << s3 << std::setw(sp) << s4				\
			 << std::setw(db) << s5	<< std::setw(sp) << s6				\
			 << std::setw(it) << s7	<< std::setw(sp) << s8				\
			 << std::setw(db) << s9	<< std::endl;						\
	}
#endif

	// Header
	LINE("", "", "", "", "cumul.", "", "calls", "", "average");
	
	// Display timings of the separate models, summing
	int counter = 0;
	for (auto const &map : profile)
		if (map.first.compare(0,8,"_NOTIME_") != 0)
		{
			counter++;
			std::stringstream s;
			s << " (" << counter << ")";
			LINE(s.str(), "", map.first, ":", map.second[0], "",
				 map.second[1], "", map.second[2]);
		}

	// Newline
	file << std::endl;
	
	// Display iteration information
	for (auto const &map : profile)
		if (map.first.compare(0,8,"_NOTIME_") == 0 )
		{
			counter++;
			std::stringstream s;
			s << " (" << counter << ")";
			LINE(s.str(), "", map.first.substr(8), ":", map.second[0], "",
				 map.second[1], "", map.second[2]);
		}
}

//------------------------------------------------------------------
std::size_t CoupledModel::getHash()
{
	TIMER_START("CoupledModel: hash...");
	// create hash function
	std::hash<double> double_hash;
	
	// first we obtain the hash of stateView
	std::size_t seed = stateView_->hash();

	// then we extend and manipulate it with the current parameter
	seed ^= double_hash(getPar()) + (seed << 6);
	TIMER_STOP("CoupledModel: hash...");
	return seed;
}

//------------------------------------------------------------------
void CoupledModel::test()
{
	INFO("CoupledModel views...");
	INFO("state:  " << stateView_->norm());
	INFO("length: " << stateView_->length());
	INFO("sol:    " << solView_->norm());
	INFO("rhs:    " << rhsView_->norm());

	std::cout << "CoupledModel: stateView..." << std::endl;
	std::cout << " length: " << stateView_->length()  << std::endl;
	std::cout << " norm:   " << stateView_->norm()    << std::endl;
}

//-----------------------------------------------------------------------------
template<typename T>
void CoupledModel::write(std::vector<T> &vector, const std::string &filename)
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
		std::cout << "WARNING (Atmosphere::write): vector is empty"
				  << __FILE__ <<  __LINE__ << std::endl;
}
