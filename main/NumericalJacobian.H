#ifndef NUMERICALJACOBIAN_H
#define NUMERICALJACOBIAN_H

#include "GlobalDefinitions.H"
#include <vector>

template<typename Model, typename VectorPtr>
class NumericalJacobian
{

	Model     model_;
    VectorPtr state_;

	// Compressed column storage for the matrix
	std::vector<double> values_;
	std::vector<int>    rowIdx_;
	std::vector<int>    colPtr_;
	double h_; // finite difference increment

public:
	NumericalJacobian(Model &model, VectorPtr state)
		:
		model_    (model),
		state_    (state),
		values_   (std::vector<double>()),
		rowIdx_   (std::vector<int>()),
		colPtr_   (std::vector<int>()),
		h_        (1e-8)
		{}

	void compute()
		{
			INFO("Computing numerical Jacobian");
			INFO(" ONLY DO THIS FOR SMALL PROBLEMS AND ON A SINGLE CORE!!");

			values_.clear();
			rowIdx_.clear();
			colPtr_.clear();
			
			model_.computeRHS();
			
		    VectorPtr Fcopy = model_.getRHS('C');
			VectorPtr Fview = model_.getRHS('V');
			
			double nrm;
			for (int j = 0; j != state_->length(); ++j)
			{
				nrm = state_->norm();
				state_->updateElement(j, h_, 1.0);
				model_.computeRHS();
				Fview->update(-1.0/h_, *Fcopy, 1.0/h_);

				// values_.push_append(Fview->nonzeros())
				
				state_->updateElement(j, -h_, 1.0);
				
				if (std::abs(nrm - state_->norm()))
					ERROR("Something went wrong in updateElement",
						  __FILE__, __LINE__);
				
				if ( j % 100  == 0)
					INFO("      " << j << " of " << state_->length());
			}
			INFO("Computing numerical Jacobian done");
		}
};


#endif
