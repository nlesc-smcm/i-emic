#ifndef GMRESSOLVERDECL_H
#define GMRESSOLVERDECL_H

// Some templated types are assumed to be shared_pointers: we use -> in calls to
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

template<typename Model, typename VectorPointer>
class GMRESSolver
{
	using Vector    = typename VectorPointer::element_type;
	using STLVector = typename std::vector<double>;
	using Matrix    = typename std::vector<STLVector>;

	Model &model_;  // We hold a reference to the model

	// Since we move around these vectors we use pointers
	VectorPointer x_; // state
	VectorPointer b_; // rhs

	bool haveInitSol_;
	bool haveRHS_;

	double tol_;     // tolerance
	double resid_;   // scaled residual norm

	int maxit_;      // max # iterations
	int m_;          // # iterations before restart
	int iter_;       // iteration counter

	bool leftPrec_;  // use left preconditioning
	
	int verbosity_;
	
public:
	// constructor
	GMRESSolver(Model &model);
	GMRESSolver(Model &model, VectorPointer x, VectorPointer b);
	// destructor
	~GMRESSolver();

	int solve();

	// get members
	VectorPointer getSolution() { return x_; }
	VectorPointer getRHS()      { return b_; }

	// set members
	void setSolution(VectorPointer x) { x_ = x; haveInitSol_ = true;}
	void setRHS(VectorPointer b) { b_ = b; haveRHS_ = true;}
	
	template<typename ParListPtr>
	void setParameters(ParListPtr pars);

	double residual();
	int getNumIters() { return iter_; }

private:
	void GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn);
	void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn);
	void Update(int k, Matrix &H, STLVector &s, std::vector<Vector> &V);
	void printIterStatus();
};

#endif
