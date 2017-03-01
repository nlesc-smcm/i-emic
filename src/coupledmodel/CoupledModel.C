#include "CoupledModel.H"
#include "Ocean.H"
#include "AtmospherePar.H"
#include "Combined_MultiVec.H"

#include <vector>
#include <memory>
#include <functional>

#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

//==================================================================
// constructor
CoupledModel::CoupledModel(std::shared_ptr<Ocean> ocean,
						   std::shared_ptr<AtmospherePar> atmos,
						   Teuchos::RCP<Teuchos::ParameterList> params)
	:
	ocean_(ocean),
	atmos_(atmos),
	
	stateView_(std::make_shared<Combined_MultiVec>
			   (ocean->getState('V'), atmos->getState('V'))),	
	solView_(std::make_shared<Combined_MultiVec>
			 (ocean->getSolution('V'), atmos->getSolution('V'))),
	rhsView_(std::make_shared<Combined_MultiVec>
			 (ocean->getRHS('V'), atmos->getRHS('V'))),
	
	parName_          (params->get("Continuation parameter",
								   "Combined Forcing")),
	solvingScheme_    (params->get("Solving scheme", 'G')),

	iterGS_           (params->get("Max GS iterations", 10)),
	toleranceGS_      (params->get("GS tolerance", 1e-1)),
	syncCtr_          (0)
{
	// Let the sub-models know our continuation parameter
	ocean_->setParName(parName_);
	atmos_->setParName(parName_);
	
	// Communicate surface landmask
	LandMask mask = ocean_->getLandMask();
	atmos_->setLandMask(mask.local); // only RCP<Epetra_IntVector> 

	// // Setup coupling blocks
	// std::vector<double> C12values;
	// std::vector<int>    C12rows;
	// ocean_->getAtmosBlock(C12values, C12rows); 
	
	// std::vector<double> C21values;
	// std::vector<int>    C21rows;
	// atmos_->getOceanBlock(C21values, C21rows);

	// // A->O and O->A coupling blocks
	// C12_ = CouplingBlock("AO", C12values, C12rows, C21rows);
	// C21_ = CouplingBlock("OA", C21values, C21rows, C12rows);
	
	// C12_.info();
	// C21_.info();

	// // Get the contribution of the atmosphere to the ocean in the Jacobian	
	// B_     = std::make_shared<std::vector<double> >(C12values);
	// rowsB_ = std::make_shared<std::vector<int> >(C12rows);
	// // Get the contribution of the ocean to the atmosphere in the Jacobian
	// C_     = std::make_shared<std::vector<double> >(C21values);
	
	// Output parameters
	INFO(*params);

	// Synchronize state
	synchronize();
	
}

//------------------------------------------------------------------
// Here we throw Epetra_Vectors around. Inside the models we should check
// that the maps are correct, which in fact checks that the domain decompositions
// are compatible.
void CoupledModel::synchronize()
{

	TIMER_START("CoupledModel: synchronize...");

	syncCtr_++; // Keep track of synchronizations
	
	// Get the atmosphere temperature
	Teuchos::RCP<Epetra_Vector> atmos = atmos_->getT();

	// Get sst restricted state vector from ocean model
	Teuchos::RCP<Epetra_Vector> sst = ocean_->getSurfaceT();
	
	// Set the atmosphere in the ocean
	ocean_->setAtmosphere(atmos);
	
	// Set the SST in the atmosphere
	atmos_->setOceanTemperature(sst);
	
	TIMER_STOP("CoupledModel: synchronize...");
}

//------------------------------------------------------------------
void CoupledModel::computeJacobian()
{
	TIMER_START("CoupledModel: compute Jacobian");
	
	// Synchronize the states
	if (solvingScheme_ != 'D') { synchronize(); }
	
	ocean_->computeJacobian();	// Ocean
	atmos_->computeJacobian();	// Atmosphere
	
	TIMER_STOP("CoupledModel: compute Jacobian");
}

//------------------------------------------------------------------
void CoupledModel::computeRHS()
{
	TIMER_START("CoupledModel compute RHS");
		
	// Synchronize the states in the fully coupled case
	if (solvingScheme_ != 'D') { synchronize(); }
	
	ocean_->computeRHS();	// Ocean
	atmos_->computeRHS(); 	// Atmosphere
	
	TIMER_STOP("CoupledModel compute RHS");
}

//====================================================================
void CoupledModel::initializeFGMRES()
{
	INFO("CoupledModel: initialize FGMRES...");

	Teuchos::RCP<Teuchos::ParameterList> solverParams_ =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams_.ptr());	
	
	INFO("CoupledModel: initialize FGMRES done");
}

//------------------------------------------------------------------
void CoupledModel::solve(std::shared_ptr<Combined_MultiVec> rhs)
{
	// Start solve
	TIMER_START("CoupledModel: solve...");
	
	if (solvingScheme_ == 'D') // fully decoupled solve
	{
		ocean_->solve(rhs->First());
		atmos_->solve(rhs->Second());
	}
	else if (solvingScheme_ == 'B') // backward block GS solve
		blockGSSolve(rhs);
	else if (solvingScheme_ == 'F') // FGMRES (Belos) on complete matrix
		FGMRESSolve(rhs);
	else
 		WARNING("(CoupledModel::Solve()) Invalid mode!",
				__FILE__, __LINE__);

	INFO("CoupledModel residual = " << computeResidual(rhs));
	
	// Update the profile after a solve
	// printProfile(profile);
	TIMER_STOP("CoupledModel: solve...");
}

//------------------------------------------------------------------
void CoupledModel::FGMRESSolve(std::shared_ptr<Combined_MultiVec> rhs)
{}

//------------------------------------------------------------------
// not implemented
void CoupledModel::blockGSSolve(std::shared_ptr<Combined_MultiVec> rhs)
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

	// double residual;
	// double old_residual = computeResidual(rhs);
	
    // Initialize solution [x1;x2] = 0
 	// std::shared_ptr<Combined_MultiVec> x = getSolution('C');
	// x->PutScalar(0.0);

	// Start iteration
	// int i;
	// for (i = 1; i <= iterGS_; ++i)
	// {
	// 	// Create -C*x1 + b2
	// 	x->linearTransformation(*C_, *rowsB_, 'O', 'A');
	// 	x->update(1, *rhs, -1);

	// 	// Solve D*x2 = -C*x1 + b2
	// 	atmos_->solve(x);

	// 	// Retrieve solution
	// 	x = getSolution('C');

	// 	// Create -B*x2 + b1
	// 	x->linearTransformation(*B_, *rowsB_, 'A', 'O');
	// 	x->update(1, *rhs, -1);

	// 	// Solve A*x1 = -B*x2 + b1
	// 	ocean_->solve(Teuchos::rcp(x.get(), false));

	// 	// Retrieve solution
	// 	x = getSolution('C');
		
	// 	// Calculate residual
	// 	residual = computeResidual(rhs);
	// 	INFO("CoupledModel: blockGS, i = " << i
	// 		 << ", ||b-Jx||/||b|| = " << residual << ", tol = " << toleranceGS_);

	// 	if (residual > old_residual)
	// 		WARNING("INCREASING RESIDUAL!", __FILE__, __LINE__);
		
	// 	if (residual < toleranceGS_)
	// 		break;
	// }

	// // Do a final solve with D
	// // Create -C*x1 + b2 and solve D*x2 = -C*x1 + b2
	// x->linearTransformation(*C_, *rowsB_, 'O', 'A');
	// x->update(1, *rhs, -1);
	// atmos_->solve(x);

	// TRACK_ITERATIONS("CoupledModel: blockGS iterations...", i);
	
	// if (i == iterGS_)
	// 	WARNING("GS tolerance not reached...", __FILE__, __LINE__);

}

//------------------------------------------------------------------
// 	out = [A B; C D] * [v1; v2]
void CoupledModel::applyMatrix(Combined_MultiVec const &v,
							   Combined_MultiVec &out, char mode)
{
	TIMER_START("CoupledModel: apply matrix...");

	// Initialize output	
	out.PutScalar(0.0);	  

	// Apply the diagonal blocks
	ocean_->applyMatrix(*v.First(),  *out.First());   // A*v1
	atmos_->applyMatrix(*v.Second(), *out.Second());  // D*v2

	if (mode == 'C') 
	{
		// Obtain temporary vector
		Combined_MultiVec z(v);
		z.PutScalar(0.0);

		// Apply coupling blocks --> Parallelize
		// C12_.applyMatrix(v.Second(), z.First());
		// C21_.applyMatrix(v.First(), z.Second());

		INFO("CoupledModel::applyMatrix not fully implemented yet!!");

		out.Update(1.0, z, 1);  
	}	
	TIMER_STOP("CoupledModel: apply matrix...");
}

//------------------------------------------------------------------
void CoupledModel::applyPrecon(Combined_MultiVec const &v,
							   Combined_MultiVec &out, char mode)
{
	TIMER_START("CoupledModel: apply preconditioner2...");	

	out.PutScalar(0.0);	// Initialize output

	if ((mode == 'C') && (iterGS_ != 0))
	{
		INFO("CoupledModel::applyPrecon coupled preconditioning not implemented yet!!");
	}
	else
	{
		ocean_->applyPrecon(*v.First(),  *out.First() );
		atmos_->applyPrecon(*v.Second(), *out.Second());
	}

	TIMER_STOP("CoupledModel: apply preconditioner2...");	
}

//------------------------------------------------------------------
double CoupledModel::computeResidual(std::shared_ptr<Combined_MultiVec> rhs)
{
	std::shared_ptr<Combined_MultiVec> r =
		std::make_shared<Combined_MultiVec>(*solView_);
	
	applyMatrix(*solView_, *r);

	double rhsNorm = Utils::norm(rhs);
	
	r->Update(1, *rhs, -1); //  b-Jx
	r->Scale(1.0 / rhsNorm);
	double relResidual = Utils::norm(r); // ||b-Jx||/||b||
	
	return relResidual;
}

//------------------------------------------------------------------
std::shared_ptr<Combined_MultiVec> CoupledModel::getSolution(char mode)
{
	// obtain solution based on mode
	if (mode == 'V') // View
		return solView_;
	else if (mode == 'C') // Copy
	{
		return std::make_shared<Combined_MultiVec>(
			ocean_->getSolution('C'),
			atmos_->getSolution('C'));
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return std::shared_ptr<Combined_MultiVec>();
	}
}

//------------------------------------------------------------------
std::shared_ptr<Combined_MultiVec> CoupledModel::getState(char mode)
{
	if (mode == 'V') // View
		return stateView_;
	else if (mode == 'C') // Copy
	{
		return std::make_shared<Combined_MultiVec>(
			ocean_->getState('C'),
			atmos_->getState('C'));
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return std::shared_ptr<Combined_MultiVec>();
	}
}

//------------------------------------------------------------------
std::shared_ptr<Combined_MultiVec> CoupledModel::getRHS(char mode)
{
	if (mode == 'V') // View
		return rhsView_;
	else if (mode == 'C') // Copy
	{
		return std::make_shared<Combined_MultiVec>(
			ocean_->getRHS('C'),
			atmos_->getRHS('C'));
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return std::shared_ptr<Combined_MultiVec>();
	}
}

//------------------------------------------------------------------
double CoupledModel::getPar()
{
	double par_ocean = ocean_->getPar(parName_);
	double par_atmos = atmos_->getPar(parName_);

	// In the case that the internal parameters are not the same, 
	// we return the maximum. This happens when we perform continuations
	// in parameters that do not exist in all models. 
	return std::max(par_ocean, par_atmos);
}

//------------------------------------------------------------------
void CoupledModel::setPar(double value)
{
	ocean_->setPar(parName_, value);
	atmos_->setPar(parName_, value);
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
}
