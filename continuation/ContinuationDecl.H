#ifndef CONTINUATIONDECL_H
#define CONTINUATIONDECL_H

//! Pseudo-arclength continuation class using an
//! Euler predictor and a Newton corrector.

template<typename ModelPtr, typename VectorPtr>
class Continuation
{
	ModelPtr  model_;
	
	VectorPtr state_;
	VectorPtr stateDot_;
	VectorPtr rhs_;
	VectorPtr dFdPar_;
	
	double par_;
	double parDot_;
	double ds_;
	double epsilon_;

	bool initialTangent_;

	int newtonIterations_;
public:
	Continuation(ModelPtr model);	
	void Step();
private:
	void CreateInitialTangent();
	void CreateTangent();
	void ComputeDFDPar();
	void EulerPredictor();
	void NewtonCorrector();
};
#endif
