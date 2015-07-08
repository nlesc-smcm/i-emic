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

	test();
}

//------------------------------------------------------------------
void CoupledModel::synchronize()
{
	INFO("CoupledModel: synchronize...");
	ocean_->setAtmosphere(*(stateView_->getStdVector()));

	// returns a copy of the sst in the ocean
	std::vector<double> sst = ocean_->getSST();
	atmosphere_->setOceanTemperature(sst);
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
	// Ocean
	ocean_->solve(Teuchos::rcp(rhs.get(), false));

	// Atmosphere
	atmosphere_->solve(rhs);

	// // Alternative solve:
	// // Let J = [A,B;C,D], x = [x1;x2], b = [b1;b2]
	// // D*w1 = b2
	// atmosphere_->solve(rhs);
	// std::shared_ptr<SuperVector> w1 =
	// 	std::make_shared<SuperVector>(atmosphere_->getSolution('C')->getStdVector());

    // // bt = b1 - B*w1 linear mapping from atmosphere to ocean
	// //               (from STL vector to Epetra vector)
	// w1->linearTransformation(B_,'A','O');
	// w1->update(1, *rhs, -1);
	
	// // A*w2 = bt 
	// std::shared_ptr<SuperVector> bt(w1);
	// ocean_->solve(bt);

	// // D*w3 = C*w2
	// std::shared_ptr<SuperVector> w3 =
	// 	std::make_shared<SuperVector>(ocean_->getSolution('C')->getEpetraVector());
	
	
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
	synchronize();
}

//------------------------------------------------------------------
void CoupledModel::test()
{
	std::cout << "CoupledModel: stateView..." << std::endl;
	std::cout << " length: " << stateView_->length()  << std::endl;
	std::cout << " norm:   " << stateView_->norm()    << std::endl;
}
