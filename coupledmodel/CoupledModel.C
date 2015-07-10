
#include "CoupledModel.H"
#include "Ocean.H"
#include "Atmosphere.H"
#include "SuperVector.H"

#include <vector>
#include <memory>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

typedef std::vector<std::vector<bool> > Graph;
//------------------------------------------------------------------
CoupledModel::CoupledModel(Graph &couplings,
						   Teuchos::RCP<Epetra_Comm> comm)
	:
	couplings_(couplings),
	comm_(comm)
{
	// Create ocean object using the parallel communicator
	ocean_ = Teuchos::rcp(new Ocean(comm_));

	// Create atmosphere object
	atmosphere_ = std::make_shared<Atmosphere>();

	stateView_ =
		std::make_shared<SuperVector>(ocean_->getState('V')->getEpetraVector(),
								 atmosphere_->getState('V')->getStdVector() );

	solView_ =
 		std::make_shared<SuperVector>(ocean_->getSolution('V')->getEpetraVector(),
								 atmosphere_->getSolution('V')->getStdVector() );

	rhsView_ =
 		std::make_shared<SuperVector>(ocean_->getRHS('V')->getEpetraVector(),
								 atmosphere_->getRHS('V')->getStdVector() );

	// Get the contribution of the atmosphere to the ocean in the Jacobian
	B_     = ocean_->getAtmosBlock();
	rowsB_ = ocean_->getSSTRows();

	// Get the contribution of the ocean to the atmosphere in the Jacobian
	C_     = atmosphere_->getOceanBlock();
}

//------------------------------------------------------------------
void CoupledModel::synchronize(double relaxation)
{
	INFO("CoupledModel: synchronize...");
	ocean_->setAtmosphere(*(stateView_->getStdVector()), relaxation);

	// returns a copy of the sst in the ocean
	atmosphere_->setOceanTemperature(*(ocean_->getSST()), relaxation);
	INFO("CoupledModel: synchronize... done");
}

//------------------------------------------------------------------
void CoupledModel::computeJacobian()
{
	// Ocean
	ocean_->computeJacobian();

	// Atmosphere
	atmosphere_->computeJacobian();
}

//------------------------------------------------------------------
void CoupledModel::computeRHS()
{
	// Ocean
	ocean_->computeRHS();

	// Atmosphere
	atmosphere_->computeRHS();
}

//------------------------------------------------------------------
void CoupledModel::solve(std::shared_ptr<SuperVector> rhs)
{
#if 0
	// Ocean
	ocean_->solve(Teuchos::rcp(rhs.get(), false));

	// Atmosphere
	atmosphere_->solve(rhs);
#endif
#if 1
	//--------------------------------------------------------------------------
	// Alternative solve:
	// Let J = [A,B;C,D], x = [x1;x2], b = [b1;b2]
	// We solve this system in an elimination based fashion. The inverse of the
	//  Schur complement is approximated using a Neumann series expansion:
	//  inv(A - B*inv(D)*C) = inv(I-inv(A)*B*inv(D)*C)*inv(A)
	//                \approx (I+inv(A)*B*inv(D)*C)*inv(A)
	// This approach results in a sequence of solves.

	// Relaxation for the coupling in the rhs
	double lambda = 1;
	
	// D*w1 = b2
	atmosphere_->solve(rhs);

	// We extract the solution w1 from the atmosphere. The solution in the ocean
	// is also obtained so that its properties are available when performing a
	// linear transformation.
	std::shared_ptr<SuperVector> w1 = getSolution('C', 'C');
	
    // btmp = b1 - lambda*B*w1 linear transformation from atmosphere to ocean
	//                             (from STL vector to Epetra vector)
	w1->linearTransformation(*B_, *rowsB_, 'A', 'O');
	w1->update(1, *rhs, -lambda);

	// A*w2 = btmp 
	std::shared_ptr<SuperVector> btmp(w1);
	ocean_->solve(Teuchos::rcp(btmp.get(), false));
	
	// extract solution from ocean solve -> w2
	// call it Cw2
	std::shared_ptr<SuperVector> Cw2 = getSolution('C', 'C');
	
	// linear transformation from ocean to atmosphere: C*w2
	Cw2->linearTransformation(*C_, *rowsB_, 'O', 'A'); 

	// get another copy of w2 for later use
	std::shared_ptr<SuperVector> w2 = getSolution('C', 'C');

#if 1
	// D*w3 = C*w2
	atmosphere_->solve(Cw2);	
	
	// extract solution from atmosphere solve -> w3
	std::shared_ptr<SuperVector> w3 = getSolution('C', 'C');

	// linear transformation from atmosphere to ocean
	w3->linearTransformation(*B_, *rowsB_, 'A','O');
	std::shared_ptr<SuperVector> Bw3(w3); // call it Bw3
	
	//  A*w4 = B*w3
	ocean_->solve(Teuchos::rcp(Bw3.get(), false));

#endif	
#if 0

//	k = 2  
	
#endif	
	// extract solution VIEW from ocean -> w4
	// and a COPY from the atmosphere for later use
	// call it x1   
	std::shared_ptr<SuperVector> x1 = getSolution('V','C');
	
	// update the solution in the ocean: x1 = w2 + w4
	x1->update(1, *w2, 1);
	
	// btmp = b2 - C*x1: linear mapping from atmosphere to ocean
	// obtain a copy of the solution x1 = w4 + w2 in the ocean:
	// call it btmp
	btmp = getSolution('C','C');

	// Get C*x1 in the stdVector of btmp:
	btmp->linearTransformation(*C_, *rowsB_, 'O','A');

	// Compute btmp = b2 - lambda*C*x1:
	btmp->update(1, *rhs, -lambda);

    // Solve D*x2 = btmp:
	atmosphere_->solve(btmp);

	// By now the ocean model will contain x1 = w2 + w4 and the atmosphere model
	// will have x2 = inv(D)*btmp.
	// The getSolution routine will be responsible for extracting these vectors.	
	//--------------------------------------------------------------------------
#endif
}

//------------------------------------------------------------------
std::shared_ptr<SuperVector> CoupledModel::getSolution(char mode)
{
	if (mode == 'V') // View
		return solView_;
	else if (mode == 'C') // Copy
	{
		return std::make_shared<SuperVector>(
			ocean_->getSolution('C')->getEpetraVector(),
			atmosphere_->getSolution('C')->getStdVector() );
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return nullptr;
	}
}

//------------------------------------------------------------------
std::shared_ptr<SuperVector> CoupledModel::getSolution(char mode1,
													   char mode2)
{
	return std::make_shared<SuperVector>(
		ocean_->getSolution(mode1)->getEpetraVector(),
		atmosphere_->getSolution(mode2)->getStdVector() );
}

//------------------------------------------------------------------
std::shared_ptr<SuperVector> CoupledModel::getState(char mode)
{
	if (mode == 'V') // View
		return stateView_;
	else if (mode == 'C') // Copy
	{
		return std::make_shared<SuperVector>(
			ocean_->getState('C')->getEpetraVector(),
			atmosphere_->getState('C')->getStdVector() );
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return nullptr;
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
			ocean_->getRHS('C')->getEpetraVector(),
			atmosphere_->getRHS('C')->getStdVector() );
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return nullptr;
	}
}

//------------------------------------------------------------------
void CoupledModel::setState(std::shared_ptr<SuperVector> state)
{
	ocean_->setState(Teuchos::rcp(state.get(), false));
    atmosphere_->setState(state);
}

//------------------------------------------------------------------
void CoupledModel::setRHS(std::shared_ptr<SuperVector> rhs)
{
	ocean_->setRHS(Teuchos::rcp(rhs.get(), false));
    atmosphere_->setRHS(rhs);
}

//------------------------------------------------------------------
double CoupledModel::getPar()
{
	// The parameters should remain equal among the models
	// Different continuation parameters for different models
	// is not defined (yet).
	double par_ocean = ocean_->getPar();
	double par_atmos = atmosphere_->getPar();
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
	atmosphere_->setPar(value);
	synchronize(1.0);
}

//------------------------------------------------------------------
double CoupledModel::getParDestination()
{
	// The parameters should remain equal among the models
	// Different continuation parameters for different models
	// is not defined (yet).
	double parDest_ocean = ocean_->getParDestination();
	double parDest_atmos = atmosphere_->getParDestination();
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
void CoupledModel::dumpState()
{
	ocean_->dumpState();
	atmosphere_->dumpState();
}

//------------------------------------------------------------------
void CoupledModel::test()
{
	std::cout << "CoupledModel: stateView..." << std::endl;
	std::cout << " length: " << stateView_->length()  << std::endl;
	std::cout << " norm:   " << stateView_->norm()    << std::endl;
}
