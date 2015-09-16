//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef CONTINUATIONDECL_H
#define CONTINUATIONDECL_H
#include <vector>

//! Pseudo-arclength continuation class using an
//! Euler predictor and a Newton corrector.
//! 
//! The templated type ParameterList should be a pointer 
//! to a parameter storing object having a .get() method.
//! ...
template<typename Model, typename ParameterList>
class Continuation
{
	using Vector = typename Model::element_type::VectorPtr;
	Model  model_;
	
	Vector stateCopy_; // copy of the state in the model
	Vector stateView_; // view of the state
	Vector rhsCopy_;   // copy of the rhs in the model
	Vector rhsView_;   // view of the rhs in the model
	Vector solView_;   // view of the solution in the model

	Vector stateDot_;  // derivative of the state w.r.t. arclength par.
	
	Vector dFdPar_;    // derivative of the rhs w.r.t. continuation par.
	
	double par_;       // continuation parameter
	double parDot_;    // derivative of cont. parameter w.r.t. arclength par.

	std::vector<double> parHist_;        // parameter history
	std::vector<double> stateNormHist_;  // state history
	
	double ds_;    // continuation step-size
	double dsMin_; // minimum continuation step-size
	double dsMax_; // maximum continuation step-size

	double zeta_;  // additional scaling for stateDot

	double scale1_; // conservative scaling for adaptive step
	double scale2_; // drastic scaling for  resetting adaptive step

	double epsilon_; // variation used for finite difference

	bool destinationReached_;   // status flag when destination reached

	bool backTracking_;         // perform backtracking when Newton fails
	int  backTrack_;            // backtracking counter
	int  numBackTrackingSteps_; // number of backtracking steps

	int step_;                 // continuation step counter
	int maxSteps_;             // maximum number of continuation steps:
	
   	int newtonIter_;              // newton iteration counter
	int maxNewtonIterations_;     // set the max for the Newton iterations
	int minAllowedIterations_;    // if less then we increase the step size
	int maxAllowedIterations_;    // if more then we decrease the step size
	double newtonTolerance_;      // set the tolerance in the Newton corrector
  bool rejectFailedNewton_;

	int parDotSign_;           // used for detecting sign switch
	
	double normRHS_;
	double normRHStest_;

	// See Store() and Restore() for its use
	struct Storage
	{
		Vector state0;  // prev state
		Vector state00; // prev prev state
		double par0;       // prev par
		double par00;      // prev prev par
		double ds0;        // prev ds
		double ds00;       // prev prev ds
	};

	Storage storage_;
	
public:
	// constructor
	Continuation(Model model, ParameterList pars);

	// destructor
	~Continuation();

	// run continuation
	void run();

	// test
	void test();
private:
	void step();
	void store();
	void restore();
	void analyzeHist();
	void testCopyView();

	// Create initial Euler tangent
	void createInitialTangent();

	// Modes:  'S' Secant
	//         'E' Euler
	void createTangent(char mode);

	// Normalize the tangent with an additional scaling
	// for the state component.
	void normalize(Vector stateDot, double scale,
				   double parDot);

	// Compute the derivative of F with respect to
	// the continuation parameter. Mode governs the
	// computation of the RHS in the model.
	// Modes: 'F' --> force compute RHS
	//        'A' --> do not force compute RHS
	void computeDFDPar(char mode = 'A');
	void eulerPredictor();
	int  newtonCorrector();
	int  runBackTracking(Vector stateDir, double parDir);

	// Detect special points.
	// Modes: 'D' --> converge on destination
	//        'P' --> converge on first sign switch in parDot
	void detect(char mode = 'D');
	void adjustStep();
	void reset();
	void info();
	void modelInfo();
	
};
#endif
