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


template<typename Model, typename VectorPointer>
class IDRSolver
{
	// We require the pointers to vectors to be of shared_ptr/RCP type
	using  Vector = typename VectorPointer::element_type;
	
	Model &model_;
	
	VectorPointer x_;
	VectorPointer b_;

	bool haveInitSol_;
	bool haveRHS_;
	bool smoothing_;
	bool replacement_;
	bool trueres_;
	int  verbosity_;

	// Specify whether or not to save the search space
	bool inispace_;

	// number close to machine precision
	double mp_;
	double tol_;
	double tolb_;
	double normr_;
	double angle_;
	
	// 
	int    s_;	
	int    maxit_;
	int    replacements_;
	int    iter_;

	// smoothing vectors
	Vector xs_;
	Vector rs_;

	// Shadow space
	std::vector<Vector> P_;

	// Initial search space
	std::vector<Vector> U_init_;

	// vector with residual norms
	std::vector<double> resvec_;

public:
	// constructor
	IDRSolver(Model &model);
	
	// constructor
	IDRSolver(Model &model,
			  VectorPointer x0,
			  VectorPointer b);
	
	// destructor
	~IDRSolver();

	template<typename ParListPtr>
	void setParameters(ParListPtr pars);
	
	int solve();

	double calc_omega(Vector const &t, Vector const &s, double angle);

	// get members
	VectorPointer getSolution() { return x_; }
	VectorPointer getRHS()      { return b_; }

	// set members
	void setSolution(VectorPointer x) { x_ = x; haveInitSol_ = true;}
	void setRHS(VectorPointer b) { b_ = b; haveRHS_ = true;}

	void clearSearchSpace()
		{
			if (verbosity_ > 0)
				std::cout << "IDR: Clearing search space!!" << std::endl;
			U_init_.clear();
		}
	void printResVec();
	void printIterStatus();
	
	void test();

	double explicitResNorm();
	double implicitResNorm() { return normr_; }

	int getNumIters() { return iter_; }

private:

	void createP();
	void writeVector(std::vector<double> &vector,
					 const std::string &filename);
	
};


#endif
