#ifndef JDQZINTERFACE_H
#define JDQZINTERFACE_H

#include "ComplexVector.H"

// Class to interface one of our models to the JDQZ++ eigenvalue solver.

template<typename Model, typename VectorType>
class JDQZInterface
{
public:
	using Vector = VectorType;

private:
	// Model (pointer)
	Model model_;

	// Problem size
	size_t n_;
	

public:
	
	JDQZInterface(Model model, VectorType v) :
		model_(model), n_(v.length()) {}
	
 	// Subroutine to compute r = Aq
	void AMUL(VectorType const &q, VectorType &r)
		{
			model_->applyMatrix(q.real, r.real);
			model_->applyMatrix(q.imag, r.imag);
		}

	// Subroutine to compute r = Bq
	void BMUL(VectorType const &q, VectorType &r)
		{
			r.real = q.real;
			r.imag = q.imag;
		}

	// Subroutine to compute q = K^-1 q
	void PRECON(VectorType &q)
		{
			// model_->applyPrecon(q.real(), q.real());
			// model_->applyPrecon(q.imag(), q.imag());
		}
	
	size_t size() { return n_; }	
};

#endif
