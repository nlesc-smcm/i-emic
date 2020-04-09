#ifndef NUMERICALJACOBIAN_H
#define NUMERICALJACOBIAN_H

#include "GlobalDefinitions.H"
#include "Utils.H"

#include <vector>

template<typename ModelPtr, typename VectorPtr>
class NumericalJacobian
{
public:
    struct CCS // compressed column storage
    {
        std::vector<double> co;
        std::vector<int>   ico;
        std::vector<int>   beg;
    };
    
private:
	// Compressed column storage for the matrix. Here, this is more
	// convenient than CRS.
	std::vector<double>     co_; // co
	std::vector<int>        ico_; // ico
	std::vector<int>        beg_; // beg

	double h_;   // finite difference increment
    double tol_; // tolerance


public:
	NumericalJacobian()
		:
		co_   (std::vector<double>()),
		ico_   (std::vector<int>()),
		beg_   (std::vector<int>()),
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

			co_.clear();
			ico_.clear();
			beg_.clear();

			model->computeRHS();

		    VectorPtr Fcopy = model->getRHS('C');
			VectorPtr Fview = model->getRHS('V');

			double nrm = Utils::norm(state);
			int elctr  = 0; // element counter

			for (int j = 0; j != state->MyLength(); ++j)
			{
                if ( j % 30  == 0)
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

				beg_.push_back(elctr);

				for (int i = 0; i != Fview->MyLength(); ++i)
					if (std::abs((*Fview)[i]) > tol_)
					{
						co_.push_back((*Fview)[i]);
						ico_.push_back(i);
						elctr++;
					}
			}

			beg_.push_back(elctr);
			INFO("Computing numerical Jacobian done");
		}

	void print(std::string const &basename)
		{
			std::ofstream co_str, ico_str, beg_str;
			co_str.open(basename + "_co");
			for (auto &i: co_)
				co_str << std::scientific << std::setprecision(12) << i << '\n';
			co_str.close();
			ico_str.open(basename + "_ico");
			for (auto &i: ico_)
				ico_str << std::scientific << std::setprecision(12) << i << '\n';
			ico_str.close();
			beg_str.open(basename + "_beg");
			for (auto &i: beg_)
				beg_str << std::scientific << std::setprecision(12) << i << '\n';
			beg_str.close();
		}

    void setTolerance(double tol) { tol_ = tol; }

    void seth(double h) { h_ = h; }

    double sumValues()
        {
            double sum = 0;
            for (auto &el: co_)
                sum += el;
            return sum;
        }
    
    void fillCCS(CCS &ccs)
        {
            ccs.co  = co_;
            ccs.ico = ico_;
            ccs.beg = beg_;
        }
};


#endif
