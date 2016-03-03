#ifndef NUMERICALJACOBIAN_H
#define NUMERICALJACOBIAN_H

#include "GlobalDefinitions.H"
#include <vector>

template<typename Model, typename Vector>
class NumericalJacobian
{
	// Compressed column storage for the matrix
	std::vector<double> values_;
	std::vector<int>    rowIdx_;
	std::vector<int>    colPtr_;
	double h_;   // finite difference increment

public:
	NumericalJacobian()
		:
		values_   (std::vector<double>()),
		rowIdx_   (std::vector<int>()),
		colPtr_   (std::vector<int>()),
		h_        (1)
		{
			INFO("NumericalJacobian: constructor");
		}

	void compute(Model &model, Vector &state)
		{
			INFO("Computing numerical Jacobian\n");
			INFO(" ONLY DO THIS FOR SMALL PROBLEMS AND ON A SINGLE CORE!!");
			INFO("  OTHERWISE I DON'T KNOW WHAT WILL HAPPEN!!\n");
			
			values_.clear();
			rowIdx_.clear();
			colPtr_.clear();
			
			model.computeRHS();
			
		    Vector Fcopy = *model.getRHS('C');
			Vector Fview = *model.getRHS('V');
			
			double nrm = state.norm();
			INFO("    state norm: " << nrm);
			INFO("      rhs norm: " << Fview.norm());
			int elctr = 0; // element counter
			for (int j = 0; j != state.length(); ++j)
			{
				state.updateElement(j, h_, 1.0);
				model.computeRHS();
				state.updateElement(j, -h_, 1.0);
				
				Fview.update(-1.0/h_, Fcopy, 1.0/h_);

				//Fview.print("F" + std::to_string(j));
				
				colPtr_.push_back(elctr);
				for (int i = 0; i != Fview.length(); ++i)
					if (std::abs(Fview[i]) > 1e-5)
					{
						values_.push_back(Fview[i]);
						rowIdx_.push_back(i);
						elctr++;
					}
				
				
				if (std::abs(nrm - state.norm()) > 1e-8)
					ERROR("Something went wrong in updateElement",
						  __FILE__, __LINE__);
				
				if ( j % 100  == 0)
					INFO("      " << j << " of " << state.length());
			}
			
			colPtr_.push_back(elctr);
			INFO("   " << std::abs(nrm - state.norm()));
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
