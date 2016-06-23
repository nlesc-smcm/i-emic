//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef TOPO_H
#define TOPO_H

#include "TopoDecl.H"
#include "GlobalDefinitions.H"

#include <vector>
#include <sstream>
#include <math.h>

//==================================================================
// constructor
template<typename Model, typename ParameterList>
Topo<Model, ParameterList>::
Topo(Model model, ParameterList pars)
	:
	model_  (model),
	pars_   (pars),
	nMasks_ (pars->get("Number of mask files", 0)),
	delta_  (pars->get("Delta", 0.0))	
{
	INFO("Topo constructor...");
	loadMaskFileNames();
	loadMasks();

	model_->setLandMask(landMasks_[2]);

	// We initialize our vectors with the ones in the model.
	stateView_ = model_->getState('V');
	solView_   = model_->getSolution('V');
	rhs_       = model_->getRHS('C');

	// Setup indices k,a and b
	k_ = static_cast<int>(std::floor(delta_));
	a_ = std::vector<int>(nMasks_);
	b_ = std::vector<int>(nMasks_);		
	a_[0] = 0;
	b_[0] = 1;
	for (int k = 0; k != nMasks_-1; ++k)
	{
		a_[k+1] = a_[k] + 2 * ((k + 1) % 2);
		b_[k+1] = b_[k] + 2 * ( k % 2 );
	}
	
	INFO("Topo constructor... done");
}

//==================================================================
// constructor
template<typename Model, typename ParameterList>
Topo<Model, ParameterList>::
~Topo()
{
	INFO("Topo destructor");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::loadMaskFileNames()
{
	INFO(" Topo: loading mask files...");
	
	landMaskFileNames_ = std::vector<std::string>(nMasks_);
	std::stringstream fileID;
	for (int i = 0; i != nMasks_; ++i)
	{
		fileID << "Mask file " << i;		
		landMaskFileNames_[i] = pars_->get(fileID.str().c_str(), "");
		assert(landMaskFileNames_[i] != "");
		INFO("  " << fileID.str() << ": " << landMaskFileNames_[i]);
		fileID.str("");
		fileID.clear();
	}
	
	INFO(" Topo: loading mask files... done");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::loadMasks()
{
	INFO(" Topo: loading mask arrays...");

	landMasks_ = std::vector<LandMask>(nMasks_);
	
	for (int i = 0; i != nMasks_; ++i)
	{
		INFO("   loading " << landMaskFileNames_[i]);
		landMasks_[i] = model_->getLandMask(landMaskFileNames_[i]);
	}
	
	INFO(" Topo: loading mask arrays... done");
}

//==================================================================
template<typename Model, typename ParameterList>
double Topo<Model, ParameterList>::getPar()
{
	return delta_;
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::setPar(double par)
{
	delta_ = par;
}

//==================================================================
template<typename Model, typename ParameterList>
typename Topo<Model, ParameterList>::VectorPtr
Topo<Model, ParameterList>::getState(char mode)
{
	if (mode == 'V')
	{
		return stateView_;
	}
	else if (mode == 'C')
	{
		VectorPtr stateCopy = model_->getState('C');
		// copy construction
		*stateCopy = *stateView_;
		return stateCopy;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}
}

//==================================================================
template<typename Model, typename ParameterList>
typename Topo<Model, ParameterList>::VectorPtr
Topo<Model, ParameterList>::getRHS(char mode)
{
	if (mode == 'V')
	{
		return rhs_;
	}
	else if (mode == 'C')
	{
		VectorPtr rhsCopy = model_->getRHS('C');
		// copy construction
		*rhsCopy = *rhs_;
		return rhsCopy;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//==================================================================
template<typename Model, typename ParameterList>
typename Topo<Model, ParameterList>::VectorPtr
Topo<Model, ParameterList>::getSolution(char mode)
{
	if (mode == 'V')
	{
		return solView_;
	}
	else if (mode == 'C')
	{
		VectorPtr solCopy = model_->getSolution('C');
		// copy construction
		*solCopy = *solView_;
		return solCopy;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::computeRHS()
{
	INFO("Topo: computeRHS...");
	// Get our current landmask index
	k_ = static_cast<int>(std::floor(delta_));

	// Calculate homotopy factors
	double facA = pow(cos(M_PI * delta_ / 2), 2);
	double facB = pow(sin(M_PI * delta_ / 2), 2);		

	// Set landmask (A) in the model
	model_->setLandMask(landMasks_[a_[k_]]);

	// Evaluate and compute RHS
	model_->computeRHS();

	// Let our rhs_ point to a copy of the RHS in the model
	*rhs_ = *model_->getRHS('C');
	INFO("   norm RHS (A) = " << rhs_->norm() << " facA = " << facA);

	// Multiply RHS (A) with its homotopy factor
	rhs_->scale(facA);

	// Set landmask (B) in the model
	model_->setLandMask(landMasks_[b_[k_]]);

	// Evaluate and compute RHS
	model_->computeRHS();
	
	// Obtain copy of the RHS in the model
	VectorPtr rhstmp = model_->getRHS('C');
	INFO("   norm RHS (B) = " << rhstmp->norm() << " facB = " << facB);

	// Multiply RHS (B) with its homotopy factor
	rhstmp->scale(facB);

	// Add RHS(B) to RHS(A)
	rhs_->update(1.0, *rhstmp, 1.0);

	INFO("   norm facB*B+facA*A = " << rhs_->norm());

	// To be safe we restore the model to the k-th landmask
	// and evaluate it.
	model_->setLandMask(landMasks_[k_]);
	model_->computeRHS();
	INFO("Topo: computeRHS... done");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::computeJacobian()
{
	// Get our current landmask index
	k_ = static_cast<int>(std::floor(delta_));

	// For now we only compute Jacobians based on the k-th
	// landmask.

	// Set k-th landmask in the model
	model_->setLandMask(landMasks_[k_]);
	model_->computeJacobian();	
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::solve(VectorPtr b)
{
	// Let the model do its thing
	model_->setLandMask(landMasks_[k_]);
	model_->solve(b);
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::preProcess()
{
	model_->preProcess();
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::postProcess()
{
	model_->postProcess();
}

#endif
