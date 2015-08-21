#include "CoupledModel.H"
#include "Ocean.H"
#include "Atmosphere.H"
#include "SuperVector.H"

#include <vector>
#include <memory>
#include <functional>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

CoupledModel::CoupledModel(Teuchos::RCP<Ocean> ocean,
						   std::shared_ptr<Atmosphere> atmos,
						   Teuchos::RCP<Teuchos::ParameterList> params)
	:
	ocean_(ocean),
	atmos_(atmos),
	syncHash_(-1),
	rhsHash_(-1),
	jacHash_(-1)
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
	rowsB_ = ocean_->getSSTRows();

	// Do some outputting to clutter the working directory
	write(*B_, "coupled_B.txt");
	write(*rowsB_, "coupled_rowsB.txt");

	// Get the contribution of the ocean to the atmosphere in the Jacobian
	C_     = atmos_->getOceanBlock();
	write(*C_, "coupled_C.txt");

	// Get parameters and flags from file, see xml for documentation
	solvingScheme_ = params->get("Solving scheme", 'E');
	kNeumann_      = params->get("Order of Neumann approximation", 1);
	iterSOR_       = params->get("SOR iterations", 2);
	relaxSOR_      = params->get("SOR relaxation", 1.0);
	useHash_       = params->get("Use hashing", true);
	
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

	INFO("CoupledModel: start sync");
	TIMER_START("CoupledModel: synchronize...");
	
	// Copy the atmosphere from the current combined (SuperVector) state
	std::vector<double> atmos(*(stateView_->getAtmosVector()));

	// Copy the SST from the ocean model
	std::vector<double> sst(*(ocean_->getSST()));
	
	// Set the atmosphere in the ocean
	ocean_->setAtmosphere(atmos);

	// Set the SST in the atmosphere
	atmos_->setOceanTemperature(sst);

	TIMER_STOP("CoupledModel: synchronize...");
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
	if ( (solvingScheme_ == 'E') ||
		 (solvingScheme_ == 'S')   ) { synchronize(); }
	
	// Ocean
	ocean_->computeJacobian();

	// Atmosphere
	atmos_->computeJacobian();
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
	if ( (solvingScheme_ == 'E') ||
		 (solvingScheme_ == 'S')   ) { synchronize(); }
	
	// Ocean
	ocean_->computeRHS();

	// Atmosphere
	atmos_->computeRHS();
}

//------------------------------------------------------------------
void CoupledModel::solve(std::shared_ptr<SuperVector> rhs)
{
	// Start solve
	if (solvingScheme_ == 'D') // fully decoupled solve
	{
		ocean_->solve(Teuchos::rcp(rhs.get(), false));
		atmos_->solve(rhs);
	}
	else if (solvingScheme_ == 'S') // backward block SOR solve
	{
		blockSORSolve(rhs);
	}
	else if (solvingScheme_ == 'E') // elimination based solve
	{
		eliminationSolve(rhs);
	}
	else
		WARNING("(CoupledModel::Solve()) Invalid mode!",
				__FILE__, __LINE__);
}

//------------------------------------------------------------------
// void CoupledModel::symBlockGSSolve(std::shared_ptr<SuperVector> rhs)
// {
// 	// Notation: J = [A,B;C,D], x = [x1;x2], b = [b1;b2]
// 	// Initialize solution
// 	std::shared_ptr<SuperVector> x = getSolution('C', 'C');
	
// }

//------------------------------------------------------------------
void CoupledModel::blockSORSolve(std::shared_ptr<SuperVector> rhs)
{
	// Notation: J = [A,B;C,D], x = [x1;x2], b = [b1;b2]	
	// Set relaxation parameter
	double w = relaxSOR_; 
	
	// Initialize solution, get a copy of what's already in there
	std::shared_ptr<SuperVector> x = getSolution('C', 'C');
	std::shared_ptr<SuperVector> y = getSolution('C', 'C');
	// Set to zero
	x->zero();
	y->zero();
	
	for (int i = 0; i < iterSOR_; ++i)
	{
		// -map the ocean to the atmosphere using matrix block C
		// -apply relaxation and add rhs
		x->linearTransformation(*C_, *rowsB_, 'O', 'A');
		x->update(w, *rhs, -w);

		if (std::abs(1-w) > 1e-5)
		{
			// - map the atmosphere to itself with D
			// - apply relaxation and add to x
			y->linearTransformation(atmos_->getJacobian());
			x->update(1-w, *y, 1);
		}

		// Solve D*x2 = -w*C*x1 + (1-w)*D*x2  + w*b2
		atmos_->solve(x);

		// Retrieving the solution (twice)
		x = getSolution('C', 'C');
		y = getSolution('C', 'C');		

		// - map the atmosphere to the ocean using matrix block B
		// - apply the relaxation and add the rhs
		x->linearTransformation(*B_, *rowsB_, 'A', 'O');
		x->update(w, *rhs, -w);

		if (std::abs(1-w) > 1e-5)
		{			
			// - map the ocean to itself with A
			// - apply relaxation and add to x
			y->linearTransformation(ocean_->getJacobian());
			x->update(1-w, *y, 1);
		}
		
		// Solve A*x1 = -w*B*x2 + (1-w)*A*x1 + w*b1
		ocean_->solve(Teuchos::rcp(x.get(), false));

		// Retrieving the solution
		x = getSolution('C', 'C');
		y = getSolution('C', 'C');
	}

    // Map the ocean to the atmosphere using matrix block C
	x->linearTransformation(*C_, *rowsB_, 'O', 'A');
	x->update(1, *rhs, -1);

	// Solve D*x2 = -C*x1 + b2
	atmos_->solve(x);
}

//------------------------------------------------------------------
void CoupledModel::eliminationSolve(std::shared_ptr<SuperVector> rhs)
{
	//.......................................................
	// Elimination based solve:
	// Let J = [A,B;C,D], x = [x1;x2], b = [b1;b2]
	//     A: ocean,
	//     D: atmosphere,
	//     B: influence of atmosphere on ocean
	//     C: influence of ocean on atmosphere
	//
	// We solve this system in an elimination based fashion. The inverse of the
	//  Schur complement is approximated using a Neumann series expansion:
	//  inv(A - B*inv(D)*C) = inv(I-inv(A)*B*inv(D)*C)*inv(A)
	//                \approx (I+inv(A)*B*inv(D)*C)*inv(A)
	// This approach results in a sequence of solves.
	//.......................................................

	// Relaxation for the coupling in the rhs --> leave it?
	double lambda = 1;

	// D*w1 = b2 ............................................
	atmos_->solve(rhs);

	// We extract the solution w1 from the atmosphere. The solution in the ocean
	// is also obtained so that its properties are available when performing a
	// linear transformation.
	std::shared_ptr<SuperVector> w1 = getSolution('C', 'C');

	// btmp = b1 - lambda*B*w1 ..............................

	// Compute B*w1:
	//   linear transformation from atmosphere to ocean
	//   from STL vector to Epetra vector)
	w1->linearTransformation(*B_, *rowsB_, 'A', 'O');
	w1->update(1, *rhs, -lambda);

	// A*w2 = btmp
	std::shared_ptr<SuperVector> btmp(w1);
	ocean_->solve(Teuchos::rcp(btmp.get(), false));

	if (kNeumann_ > 0)
	{
		// ......................................................
		// Neumann k = 1
		// extract solution from ocean solve -> w2
		// We need 2 copies:
		//  -one to manipulate and get Cw2 and
		//  -one to store and use in the calculation of x1

		std::shared_ptr<SuperVector> w2 = getSolution('C', 'C');

		// call this one Cw2
		std::shared_ptr<SuperVector> Cw2 = getSolution('C', 'C');

		// Compute linear transformation from ocean to atmosphere: C*w2
		Cw2->linearTransformation(*C_, *rowsB_, 'O', 'A');

		// D*w3 = C*w2 ..........................................
		atmos_->solve(Cw2);

		// extract solution from atmosphere solve -> w3
		std::shared_ptr<SuperVector> w3 = getSolution('C', 'C');

		// linear transformation from atmosphere to ocean
		w3->linearTransformation(*B_, *rowsB_, 'A','O');
		std::shared_ptr<SuperVector> Bw3(w3); // call it Bw3

		//  A*w4 = B*w3
		ocean_->solve(Teuchos::rcp(Bw3.get(), false));

		if (kNeumann_ > 1)
		{
			// ......................................................
			// Neumann k = 2
			// We need 2 copies:
			//  -one to manipulate and get Cw4 and
			//  -one to store and use in the calculation of x1
			std::shared_ptr<SuperVector> w4  = getSolution('C','C');
			std::shared_ptr<SuperVector> Cw4 = getSolution('C','C');
			Cw4->linearTransformation(*C_, *rowsB_, 'O', 'A');
			atmos_->solve(Cw4);

			std::shared_ptr<SuperVector> w5 = getSolution('C','C');
			w5->linearTransformation(*B_, *rowsB_, 'A', 'O');
			std::shared_ptr<SuperVector> Bw5(w5);

			ocean_->solve(Teuchos::rcp(Bw5.get(),false));

			// x1 = w2 + w4 + w6
			std::shared_ptr<SuperVector> x1 = getSolution('V','C');
			x1->update(1, *w4, 1);
			x1->update(1, *w2, 1);
		}
		else
		{
			// x1 = w2 + w4
			std::shared_ptr<SuperVector> x1 = getSolution('V','C');
			x1->update(1, *w2, 1);
		}
	}
	else
	{
		// x1 = w2 -> already in ocean_, no need for extraction
	}
	// btmp = b2 - C*x1: linear mapping from atmosphere to ocean
	// obtain a copy of the solution x1 = w4 + w2 in the ocean:
	// call it btmp
	btmp = getSolution('C','C');

	// Get C*x1 in the stdVector of btmp:
	btmp->linearTransformation(*C_, *rowsB_, 'O','A');

	// Compute btmp = b2 - lambda*C*x1:
	btmp->update(1, *rhs, -lambda);

	// Solve D*x2 = btmp:
	atmos_->solve(btmp);

	// By now the ocean model will contain x1 and the atmosphere model
	// will have x2 = inv(D)*btmp.
	// The getSolution routine will be responsible for extracting these vectors.
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
		WARNING("par_ocean != par_atmos !!",
				__FILE__, __LINE__);
 		std::cout << "ocean: " << par_ocean << std::endl;
		std::cout << "atmos: " << par_atmos << std::endl;
		return -1;
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
void CoupledModel::postConvergence()
{
	// If the solver is completely decoupled, this is the right
	// moment to synchronize
	if (solvingScheme_ == 'D')
		synchronize();

	// Let the models do their post-convergence processing
	ocean_->postConvergence();
	atmos_->postConvergence();
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
