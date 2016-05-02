#ifndef JDQZINTERFACE_H
#define JDQZINTERFACE_H

#include "SuperVector.H"

template<typename Model>
class JDQZInterface
{
public:
	// Define vector type
	using Vector = ComplexSuperVector;
	
private:

	// Model (pointer)
	Model model_;

	// Problem size
	size_t n_;
	

public:
	
	JDQZInterface(Model model, Vector v) :
		model_(model), n_(v.length()) {}
	
 	// Subroutine to compute r = Aq
	void AMUL(Vector const &q, Vector &r)
		{
			model_->applyMatrix(q.real, r.real);
			model_->applyMatrix(q.imag, r.imag);
		}

	// Subroutine to compute r = Bq
	void BMUL(Vector const &q, Vector &r)
		{
			r.real = q.real;
			r.imag = q.imag;
		}

	// Subroutine to compute q = K^-1 q
	void PRECON(Vector &q)
		{
			// model_->applyPrecon(q.real(), q.real());
			// model_->applyPrecon(q.imag(), q.imag());
		}
	
	size_t size() { return n_; }	
};

#endif
