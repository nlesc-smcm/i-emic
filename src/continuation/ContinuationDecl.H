//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef CONTINUATIONDECL_H
#define CONTINUATIONDECL_H
#include <vector>

//! Pseudo-arclength continuation class using an
//! Euler predictor and a Newton corrector.
//! 
//! The templated Model type should be a pointer to a model
//! with a specific set of member functions:
//! ...
//!
//! The templated type ParameterList should be a pointer 
//! to a parameter storing object having a .get() method.
//!
//! --> At this point the differences between the responsibilities of
//! routines run() and step() are vague, we should improve this.

template<typename Model, typename ParameterList>
class Continuation
{
	using Vector = typename Model::element_type::VectorPtr;
	Model model_;
	
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

	double parDest_; // parameter destination

	double zeta_;  // additional scaling for stateDot

	double scale1_; // conservative scaling for adaptive step
	double scale2_; // drastic scaling for  resetting adaptive step

	double epsilon_; // variation used for finite difference

	bool destinationReached_;   // status flag when destination reached
	bool secant_;               // status flag giving detect() power over ds

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
	double destinationTolerance_; // set the tolerance for the final convergence
	double tanScaling_;           // scaling of the state tangent
	bool rejectFailedNewton_;     // set to false if you feel lucky
	double practicalTol_;         // actual tolerance used for the problem
	bool cycleTolerance_;         // if true we cycle between multiples of tol
	bool usePracticalTol_;        // if true we scale the tolerance
	bool scaleRHS_;               // experimental rhs scaling

	// If the parameter dependence in the rhs is constant we can reuse statedot in the Newton corrector,
	// effectively disregarding updates in the Jacobian. This means we do a partial Newton-chord iteration.
	bool newtChordHybr_;
	
	int parDotSign_;              // used for detecting sign switch
	int signMonitor_;             // track the sign of monitor in detect()
	
	double normRHS_;              // ||F(x_old)||
	double normRHStest_;          // ||F(x_new)||

	// See Store() and Restore() for its use
	struct Storage
	{
		Vector state0;     // prev state
		Vector state00;    // prev prev state
		double par0;       // prev par
		double par00;      // prev prev par
		double ds0;        // prev ds
		double ds00;       // prev prev ds
	};

	Storage storage_;
	
public:

	// default constructor
	Continuation(){};
	
	// constructor
	Continuation(Model model, ParameterList pars);

	// destructor
	~Continuation();

	// run continuation
	void run();

	// test
	void test();
private:
        int  step();
	void store();
	void restore();
	void analyzeHist();
	void testCopyView();

	// Create initial Euler tangent
	void createInitialTangent();

	// Modes:  'S' Secant
	//         'E' Euler
	void createTangent(char mode);
	void computeTolerance();

	// Normalize the tangent with an additional scaling
	// for the state component.
	void normalize();

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
