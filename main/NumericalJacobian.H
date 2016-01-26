#ifndef NUMERICALJACOBIAN_H
#define NUMERICALJACOBIAN_H

#include "GlobalDefinitions.H"
#include <vector>

template<typename Model, typename VectorPtr>
class NumericalJacobian
{

	Model     model_;
    VectorPtr state_;
	std::vector<VectorPtr> tangents_;
	double h_; // finite difference increment

public:
	NumericalJacobian(Model &model, VectorPtr state)
		:
		model_    (model),
		state_    (state),
		tangents_ (std::vector<VectorPtr>()), 
		h_        (1e-8)
		{}

	void compute()
		{
			INFO("Computing numerical Jacobian");
			INFO("  this may take a while...  ");
			model_.computeRHS();
		    VectorPtr Fcopy = model_.getRHS('C');
			VectorPtr Fview = model_.getRHS('V');
			for (int j = 0; j != state_->length(); ++j)
			{
				state_->updateElement(j, h_, 1.0);
				model_.computeRHS();
				Fview->update(-1.0/h_, *Fcopy, 1.0/h_);
				Fview->print("Fcol" + std::to_string(j));
				state_->updateElement(j, -h_, 1.0);
				if ( j % 100  == 0)
					INFO("      " << j << " of " << state_->length());
			}
			INFO("Computing numerical Jacobian done");
		}

	
};


#endif
