#ifndef NUMERICALJACOBIAN_H
#define NUMERICALJACOBIAN_H

#include "GlobalDefinitions.H"
#include <vector>

template<typename ModelPtr, typename VectorPtr>
class NumericalJacobian
{
	// Compressed column storage for the matrix
	std::vector<double> values_;
	std::vector<int>    rowIdx_;
	std::vector<int>    colPtr_;
	double h_;   // finite difference increment
    double tol_; // tolerance

public:
	NumericalJacobian()
		:
		values_   (std::vector<double>()),
		rowIdx_   (std::vector<int>()),
		colPtr_   (std::vector<int>()),
		h_        (1e2),
        tol_      (1e-5)
		{
			INFO("NumericalJacobian: constructor");
		}

	void compute(ModelPtr model, VectorPtr state)
		{
			INFO("Computing numerical Jacobian\n");
			INFO(" ONLY DO THIS FOR SMALL PROBLEMS AND ON A SINGLE CORE!!");
			INFO("  OTHERWISE I DON'T KNOW WHAT WILL HAPPEN!!\n");

			values_.clear();
			rowIdx_.clear();
			colPtr_.clear();
			
			model->computeRHS();
			
		    VectorPtr Fcopy = model->getRHS('C');
			VectorPtr Fview = model->getRHS('V');

			double nrm = Utils::norm(state);
			int elctr  = 0; // element counter

			for (int j = 0; j != state->GlobalLength(); ++j)
			{
                if ( j % 6  == 0)
                {
                    std::cout << std::endl << "*-*-*-*- Num.Jac. column " << j
                              << " of " << state->GlobalLength()
                              << "-*-*-*-*" << std::endl;
					INFO(std::endl << "*-*-*-*- Num.Jac. column " << j
                         << " of " << state->GlobalLength()
                         << "-*-*-*-*" << std::endl);
                }

                (*state)[j] += h_;

				Fview->PutScalar(0.0);
				model->computeRHS();

                (*state)[j] -= h_;
                
				assert(std::abs(nrm - Utils::norm(state)) < 1e-12);
				
				Fview->Update(-1.0/h_, *Fcopy, 1.0/h_);
				
				colPtr_.push_back(elctr);

				for (int i = 0; i != Fview->GlobalLength(); ++i)
					if (std::abs((*Fview)[i]) > tol_)
					{
						values_.push_back((*Fview)[i]);
						rowIdx_.push_back(i);
						elctr++;
					}								
				
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

    void setTolerance(double tol) { tol_ = tol; }

    void seth(double h) { h_ = h; }
};


#endif
