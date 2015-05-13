#ifndef CONTINUATIONDECL_H
#define CONTINUATIONDECL_H

//! Pseudo-arclength continuation class using an
//! Euler predictor and a Newton corrector.

template<typename ModelPtr, typename VectorPtr>
class Continuation
{
	ModelPtr  model_;
	
	// A copy of the state in the model
	VectorPtr stateCopy_;
	VectorPtr stateView_;
	VectorPtr stateDot_;
	VectorPtr rhsCopy_;
	VectorPtr rhsView_;
	VectorPtr dFdPar_;
	
	double par_;
	double parDot_;
	double ds_;
	double dsMin_;
	double dsMax_;
	double zeta_;
	double scale1_;
	double scale2_;
	double epsilon_;

	bool destinationReached_;
	bool backTracking_;

	int step_;
	int maxSteps_;
	int iter_;
	int maxNewtonIterations_;
	int minAllowedIterations_;
	int maxAllowedIterations_;
	int backTrack_;
	int numBackTrackingSteps_;
	int parDotSign_;
	
	double newtonTolerance_;
	double normRHS_;
	double normRHStest_;

	// See Store() and Restore() for its use
	struct Storage
	{
		VectorPtr state0;  // prev state
		VectorPtr state00; // prev prev state
		double par0;       // prev par
		double par00;      // prev prev par
		double ds0;        // prev ds
		double ds00;       // prev prev ds
	};

	Storage storage_;
	
public:
	Continuation(ModelPtr model);
	void Run();
private:
	void Step();
	void Store();
	void Restore();

	// Create initial Euler tangent
	void CreateInitialTangent();

	// Modes:  'S' Secant
	//         'E' Euler
	void CreateTangent(char mode);

	// Normalize the tangent with an additional scaling
	// for the state component.
	void Normalize(VectorPtr stateDot, double scale,
				   double parDot);

	// Compute the derivative of F with respect to
	// the continuation parameter. Mode governs the
	// computation of the RHS in the model.
	// Modes: 'F' --> force compute RHS
	//        'A' --> do not force compute RHS
	void ComputeDFDPar(char mode = 'A');
	void EulerPredictor();
	int  NewtonCorrector();
	int  RunBackTracking(VectorPtr stateDir, double parDir);

	// Detect special points.
	// Modes: 'D' --> converge on destination
	//        'P' --> converge on first sign switch in parDot
	void Detect(char mode = 'D');
	void AdjustStep();
	void Reset();
	void Info();
	
};
#endif
