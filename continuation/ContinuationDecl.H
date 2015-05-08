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
	VectorPtr stateDir_;
	VectorPtr rhs_;
	VectorPtr dFdPar_;
	
	double par_;
	double parDot_;
	double parDir_;
	double ds_;
	double epsilon_;

	bool destinationReached_;
	bool initialTangent_;
	bool backTracking_;

	int step_;
	int maxSteps_;
	int iter_;
	int newtonIterations_;
	int backTrack_;
	int numBackTrackingSteps_;
	
	double newtonTolerance_;
	double normRHS_;
	double normRHStest_;
	
public:
	Continuation(ModelPtr model);
	void Run();
private:
	void Step();
	void CreateInitialTangent();
	void CreateTangent();

	// Compute the derivative of F with respect to
	// the continuation parameter. Mode governs the
	// computation of the RHS in the model.
	// Modes: 'F' --> force computing RHS
	//        'N' --> let Newton iterations decide
	void ComputeDFDPar(char mode = 'N');
	void EulerPredictor();
	void NewtonCorrector();
	void RunBackTracking();

	// Detect special points.
	// Modes: 'D' --> converge on destination
	//        'P' --> converge on first sign switch in parDot
	void Detect(char mode = 'D');
};
#endif
