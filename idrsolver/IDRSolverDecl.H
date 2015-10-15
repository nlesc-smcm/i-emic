#ifndef IDRSolverDecl_H
#define IDRSolverDecl_H

// Templated types are assumed to be shared_pointers: we use -> in calls to
// their members.

// Model should be a class with members:
//    -applyMatrix(Vector v), performing matrix vector product t=Av
//    -applyPrecon(Vector v), applying the operation v = P^{-1} x
// Note that Model should be compatible with Vector

// Vector should be a class with members:
//    -update(double scalarA, Vector A, double scalarThis, performing
//      this = scalarA * A + scalarThis * this
//    -norm()
//    -copy construction


template<typename ModelPointer, typename VectorPointer>
class IDRSolver
{
	// We require the pointers to be of shared_ptr/RCP type
	using  Model  = typename ModelPointer::element_type;	
	using  Vector = typename VectorPointer::element_type;
	
	ModelPointer  model_;
	VectorPointer x_;
	VectorPointer b_;

	// number close to machine precision
	double mp_;
	double tol_;
	double angle_;
	bool   smoothing_;
	int    s_;
	int    maxit_;

	// smoothing vectors
	Vector xs_;
	Vector rs_;

	// shadow space?
	std::vector<Vector> P_;

public:
	// constructor
	IDRSolver(ModelPointer model, VectorPointer x0);
	
	// destructor
	~IDRSolver();
	
	VectorPointer solve(VectorPointer b);

	double calc_omega(Vector const &t, Vector const &s, double angle)
	
	void test();

private:

	void createP();
	
};


#endif
