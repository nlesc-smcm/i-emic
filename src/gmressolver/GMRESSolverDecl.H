#ifndef GMRESSOLVERDECL_H
#define GMRESSOLVERDECL_H

// Some templated types are assumed to be shared_pointers: we use -> in calls to
// their members.
//
// Model should be a class with members:
//    -applyMatrix(Vector v, Vector x), performing matrix vector product x=A*v
//    -applyPrecon(Vector v, Vector x), applying the operation x=inv(P)*v
// Model should be compatible with Vectors
//
// Vector should be a class with:
//    -update(double scalarA, Vector A, double scalarThis), performing
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

	// Since we move these vectors around we use (shared) pointers
	VectorPointer x_; // state
	VectorPointer b_; // rhs

	Vector xCopy_;    // local copy of the state

	STLVector y_;
	
	bool haveInitSol_;
	bool haveRHS_;

	char minimizeScheme_;   // 'B' backward solve after Givens rotations
	                        // 'Q' QR solve using lapack routine

	bool flexible_;         // When using iterative solvers in the (right) preconditioner
	                        // the user should use FlexibleGMRES.

	bool computeExplResid_; // Choose to compute explicit residual (can be expensive)
	
	double tol_;            // tolerance
	double resid_;          // scaled residual norm
	double explResid_;      // explicit residual norm

	int maxit_;             // max # iterations
	int m_;                 // # iterations before restart
	int iter_;              // iteration counter
	
	bool prec_;             // use preconditioning
	bool leftPrec_;         // use left preconditioning
	                        // (default is right preconditioning)
	
	int verbosity_;         // 0 is quiet
	
public:
	// constructors
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
	void Update(int last, Matrix &H, STLVector &s, std::vector<Vector> &V);
	void UpdateCopy(int last, Matrix &H, STLVector &s, std::vector<Vector> &V);

	void backSolve(int m, Matrix &H, STLVector &s);
	void LLSSolve(int m, Matrix &H, STLVector &s); // uses lapack
	
	void printIterStatus();
	
	void   compute_y(int last, Matrix &H, STLVector &s);
	double compute_r(int last, Matrix &H, STLVector &s);
	double compute_explicit_residual(int last, Matrix &H, STLVector &s, std::vector<Vector> &V);
};

#endif
