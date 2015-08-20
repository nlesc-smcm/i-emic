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
	
	// copy of the state in the model
	Vector stateCopy_;
    // view of the state
	Vector stateView_;
	Vector stateDot_;
	Vector rhsCopy_;
    // view of the rhs in the model
	Vector rhsView_;
	// view of the sol in the model
	Vector solView_;
	
	Vector dFdPar_;
	
	double par_;
	double parDot_;

	// continuation history
	std::vector<double> parHist_;
	std::vector<double> stateNormHist_;
	
    // initial continuation step-size
	double ds_;
    // minimum continuation step-size
	double dsMin_;
    // maximum continuation step-size
	double dsMax_;
    // initialize additional scaling
	double zeta_;
    // conservative scaling for adaptive step
	double scale1_;
    // drastic scaling for  resetting adaptive step
	double scale2_;
    // variation used for finite difference
	double epsilon_;

	bool destinationReached_;
    // perform backtracking when Newton fails
	bool backTracking_;

	int step_;
    // maximum number of continuation steps:
	int maxSteps_;
    // newton iteration counter
	int newtonIter_;
    // set the max for the Newton iterations
	int maxNewtonIterations_;
    // if less then we increase the step size
	int minAllowedIterations_;
    // if more then we decrease the step size
	int maxAllowedIterations_;
	// backtracking counter
	int backTrack_;
    // number of backtracking steps
	int numBackTrackingSteps_;
	// used for detecting sign switch
	int parDotSign_;
	
    // set the tolerance in the Newton corrector
	double newtonTolerance_;
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
