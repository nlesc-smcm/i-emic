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
	double dsMin_;
	double dsMax_;
	double dsScale1_;
	double dsScale2_;
	double dsOld_;
	double epsilon_;

	bool destinationReached_;
	bool initialTangent_;
	bool backTracking_;

	int step_;
	int maxSteps_;
	int iter_;
	int newtonIterations_;
	int minAllowedIterations_;
	int maxAllowedIterations_;
	int backTrack_;
	int numBackTrackingSteps_;
	int parDotSign_;
	
	double newtonTolerance_;
	double normRHS_;
	double normRHStest_;
	
public:
	Continuation(ModelPtr model);
	void Run();
private:
	void Step();
	void Store();
	void Restore();
	void CreateInitialTangent();
	void CreateTangent();

	// Compute the derivative of F with respect to
	// the continuation parameter. Mode governs the
	// computation of the RHS in the model.
	// Modes: 'F' --> force computing RHS
	//        'A' --> do not force computing RHS
	void ComputeDFDPar(char mode = 'A');
	void EulerPredictor();
	int  NewtonCorrector();
	int  RunBackTracking();

	// Detect special points.
	// Modes: 'D' --> converge on destination
	//        'P' --> converge on first sign switch in parDot
	void Detect(char mode = 'D');
	void AdjustStep();
	void Reset();
	
};
#endif
