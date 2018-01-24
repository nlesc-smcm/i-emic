#ifndef JDQZOPERATORINTERFACE_H
#define JDQZOPERATORINTERFACE_H

// Class to interface the three key (Epetra) operators to JDQZ++.

template<typename Matrix, typename MassMatrix,
         typename Precon, typename VectorType>
class JDQZOperatorInterface
{
    //! Matrix
    Matrix     A_;

    //! Mass matrix
    MassMatrix B_;

    //! Preconditioner
    Precon     M_;

    //! Problem size
    size_t n_;

    //! Temporary vector
    VectorType tmp_;

public:
    //! constructor
    JDQZOperatorInterface(Matrix A, MassMatrix B,
                          Precon M, VectorType v)
        A_(A), B_(B), M_(M), tmp_(v), n_(v.length()) {}

    // Subroutine to compute r = Aq
	void AMUL(VectorType const &q, VectorType &r)
		{
			A_->Apply(q.real, r.real);
			A_->Apply(q.imag, r.imag);
		}

	// Subroutine to compute r = Bq
	void BMUL(VectorType const &q, VectorType &r)
		{
            r.real.Multiply(1.0, *B_, q.real, 0.0);
            r.imag.Multiply(1.0, *B_, q.imag, 0.0);
		}

	// Subroutine to compute q = K^-1 q
	void PRECON(VectorType &q)
		{
            tmp_.zero();
            M_->ApplyInverse(q.real, tmp_.real);
            M_->ApplyInverse(q.imag, tmp_.imag);
            q = tmp_;
		}
	
	size_t size() { return n_; }

};


#endif
