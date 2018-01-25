#ifndef COMPLEX_VEC_H
#define COMPLEX_VEC_H

#include "Utils.H"
#include "GlobalDefinitions.H"

template<typename Vector>
class ComplexVector
{
	Vector tmp_;

public:
	Vector real;
	Vector imag;

	//! constructor accepting a single Vector, which will form the
	//! real part
	ComplexVector(Vector const &re)
		:
		tmp_(re),
		real(re),
		imag(re)
		{
			imag.PutScalar(0.0);
		}

	//! constructor accepting a real and an imaginary part
	ComplexVector(Vector const &re, Vector const &im)
		:
		tmp_(re),
		real(re),
		imag(im)
		{
			assert(real.GlobalLength() == imag.GlobalLength());
		}

	//! constructor accepting a real and an imaginary part
	ComplexVector(ComplexVector const &other)
		:
		tmp_(other.real),
		real(other.real),
		imag(other.imag)
		{
			assert(real.GlobalLength() == imag.GlobalLength());
		}

    //! destructor
    ~ComplexVector()
        {}
    
		
	//! assignment operator
	void operator=(ComplexVector const &other)
		{
			real = other.real;
			imag = other.imag;
		}

	//! obtain global length
	int length() const {return real.GlobalLength();}

	//! calculate 2-norm
	double norm() const
		{
			return sqrt(dot(*this).real());			
		}

	//! calculate dot product
	std::complex<double> dot(ComplexVector const &other) const
		{
			double rePart = Utils::dot(&real, &other.real) +
				Utils::dot(&imag, &other.imag);
			double imPart = Utils::dot(&real, &other.imag) -
				Utils::dot(&imag, &other.real);
			std::complex<double> result(rePart, imPart);
			return result;			
		}

	//! this = a*x + this
	void axpy(std::complex<double> a, ComplexVector const &x)
		{
			assert(this->length() == x.length());

			real.Update( a.real(), x.real, 1.0);
			real.Update(-a.imag(), x.imag, 1.0);

			imag.Update( a.imag(), x.real, 1.0);
			imag.Update( a.real(), x.imag, 1.0);

		}

	//! this = a*x + b*this
	void axpby(std::complex<double> a, ComplexVector const &x,
			   std::complex<double> b)
		{
			assert(this->length() == x.length());
			scale(b);
			axpy(a, x);			
		}

	//! scale vector with complex number
	void scale(std::complex<double> a)
		{
			tmp_ = real; // copy
			real.Update(-a.imag(), imag,  a.real());
			imag.Update( a.imag(), tmp_,  a.real());			
		}

	//! zero-out vector
	void zero()
		{
			real.PutScalar(0.0);
			imag.PutScalar(0.0);
		}

	//! randomize real part
	//! zero imaginary part
	void random()
		{
			real.Random();
		    imag.PutScalar(0.0);
		}
	
	void print(){ Utils::print(&real, "real"); Utils::print(&imag, "imag"); }
	
};

#endif
