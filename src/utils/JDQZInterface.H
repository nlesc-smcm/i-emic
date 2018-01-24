#ifndef JDQZINTERFACE_H
#define JDQZINTERFACE_H

#include "GlobalDefinitions.H"

//! Class to interface one of our models to the JDQZ++ eigenvalue solver.

template<typename Model, typename VectorType>
class JDQZInterface
{
public:
	using Vector = VectorType;

private:
	//! Model (pointer)
	Model model_;

	//! Problem size
	size_t n_;

    //! Temporary vector
    VectorType tmp_;
	

public:
    //! constructor 
	JDQZInterface(Model model, VectorType v) :
		model_(model), n_(v.length()), tmp_(v) {}

    //! destructor
    ~JDQZInterface()
        {
            INFO("JDQZInterface destructor called...");
        }
	
 	//! Subroutine to compute r = Aq
	void AMUL(VectorType const &q, VectorType &r)
		{
			model_->applyMatrix(q.real, r.real);
			model_->applyMatrix(q.imag, r.imag);
		}

	//! Subroutine to compute r = Bq
	void BMUL(VectorType const &q, VectorType &r)
		{
            model_->applyMassMat(q.real, r.real);
            model_->applyMassMat(q.imag, r.imag);
		}

	//! Subroutine to compute q = K^-1 q
	void PRECON(VectorType &q)
		{
            tmp_.zero();
			model_->applyPrecon(q.real, tmp_.real);
			model_->applyPrecon(q.imag, tmp_.imag);
            q = tmp_;
		}
	
	size_t size() { return n_; }
};

#endif
