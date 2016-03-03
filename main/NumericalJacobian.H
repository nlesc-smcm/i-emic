#ifndef NUMERICALJACOBIAN_H
#define NUMERICALJACOBIAN_H

#include "GlobalDefinitions.H"
#include <vector>

template<typename Model, typename Vector>
class NumericalJacobian
{
	
	Model  &model_;
    Vector &state_;

	// Compressed column storage for the matrix
	std::vector<double> values_;
	std::vector<int>    rowIdx_;
	std::vector<int>    colPtr_;
	double h_; // finite difference increment

public:
	NumericalJacobian(Model &model, Vector &state)
		:
		model_    (model),
		state_    (state),
		values_   (std::vector<double>()),
		rowIdx_   (std::vector<int>()),
		colPtr_   (std::vector<int>()),
		h_        (1e-4)
		{}

	void compute()
		{
			INFO("Computing numerical Jacobian");
			INFO(" ONLY DO THIS FOR SMALL PROBLEMS AND ON A SINGLE CORE!!");
			
			values_.clear();
			rowIdx_.clear();
			colPtr_.clear();
			
			model_.computeRHS();
			
		    Vector Fcopy = *model_.getRHS('C');
			Vector Fview = *model_.getRHS('V');
			
			double nrm = state_.norm();
			int elctr = 0; // element counter
			for (int j = 0; j != state_.length(); ++j)
			{
				state_.updateElement(j, h_, 1.0);
				model_.computeRHS();
				Fview.update(-1.0, Fcopy, 1.0);
				Fview.scale(1.0/h_);

				colPtr_.push_back(elctr);
				for (int i = 0; i != Fview.length(); ++i)
					if (std::abs(Fview[i]) > 1e-5)
					{
						values_.push_back(Fview[i]);
						rowIdx_.push_back(i);
						elctr++;
					}
				
				state_.updateElement(j, -h_, 1.0);
				
				if (std::abs(nrm - state_.norm()))
					ERROR("Something went wrong in updateElement",
						  __FILE__, __LINE__);
				
				if ( j % 100  == 0)
					INFO("      " << j << " of " << state_.length());
			}
			colPtr_.push_back(elctr);
			INFO("Computing numerical Jacobian done");
		}
	
	void print(std::string const &basename)
		{
			std::ofstream co_str, ico_str, beg_str;
			co_str.open(basename + "_co");
			for (auto &i: values_)
				co_str << i << '\n';			
			co_str.close();
			ico_str.open(basename + "_ico");
			for (auto &i: rowIdx_)
				ico_str << i << '\n';			
			ico_str.close();
			beg_str.open(basename + "_beg");
			for (auto &i: colPtr_)
				beg_str << i << '\n';			
			beg_str.close();
		}
};


#endif
