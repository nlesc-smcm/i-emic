/*======================================================================
   Author T.E.Mulder 2015/16/17/18 (t.e.mulder@uu.nl)
   ======================================================================*/
#ifndef CONTINUATIONDECL_H
#define CONTINUATIONDECL_H

#include <vector>

#include "ComplexVector.H"
#include "JDQZInterface.H"

#ifdef HAVE_JDQZPP
#include "jdqz.hpp"
#else
template<typename Matrix>
class JDQZ;
#endif

//! Pseudo-arclength continuation class using an
//! Euler predictor and a Newton corrector.
//!
//! The templated Model type should be a pointer to a model
//! with a specific set of member functions:
//!
//!  VectorPtr getState()
//!  VectorPtr getRHS()
//!  VectorPtr getSolution()
//!  ...
//!
//!  double getPar()
//!  void   setPar()
//!  ...
//!
//!  void preProcess()
//!  void postProcess()
//!  ...
//!
//!  void computeRHS()
//!  void computeJacobian()
//!  void solve()
//!  ...
//!
//! A Model should maintain its own Vector, which we expect
//! to be of pointer type as well and named VectorPtr.
//! A Vector is expected to have an Epetra_MultiVector style interface:
//!
//!  void    Update()         update
//!  void    Print()          print to ostream
//!  int     GlobalLength()   length of vector
//!  double  Norm()           2-norm
//!  double  Dot()            dot
//!  ...
//!

template<typename Model>
class Continuation
{
    using VectorPtr  = typename Model::element_type::VectorPtr;
    using Vector     = typename VectorPtr::element_type;
    using JDQZsolver = JDQZ<JDQZInterface<Model,
                                          ComplexVector<Vector> > >;
    Model model_;
    Teuchos::ParameterList paramList_;

    //! copy of the state in the model
    VectorPtr stateCopy_;
    //! view of the state
    VectorPtr stateView_;
    //! copy of the rhs in the model
    VectorPtr rhsCopy_;
    //! view of the rhs in the model
    VectorPtr rhsView_;
    //! view of the solution in the model
    VectorPtr solView_;

    //! derivative of the state w.r.t. arclength par.
    VectorPtr stateDot_;

    //! derivative of the rhs w.r.t. continuation par.
    VectorPtr dFdPar_;

    //! continuation parameter name
    std::string parName_;
    //! continuation parameter value
    double par_;
    //! starting value continuation parameter
    double startingPar_;
    //! derivative of cont. parameter w.r.t. arclength par.
    double parDot_;

    //! parameter history
    std::vector<double> parHist_;
    //! state history
    std::vector<double> stateNormHist_;

    //! continuation step-size
    double ds_;
    //! initial continuation step-size
    double dsInit_;
    //! minimum continuation step-size
    double dsMin_;
    //! maximum continuation step-size (absolute)
    double dsMax_;
    //! step size at beginning of secant process
    double dsStart_;

    //! parameter destinations
    std::vector<double> destinations_;       // list we adjust
    std::vector<double> destinationsBackup_; // list to restore to

    //! maximum allowed number of destinations
    static const int maxNumDest_ = 10;

    //! zeta is an additional scaling for the state tangent
    double zeta_;

    //! conservative scaling for adaptive step
    double scale1_;
    //! drastic scaling for adaptive step reset
    double scale2_;

    //! variation used for finite difference
    double epsilon_;

    //! status flag when aborting
    bool abortFlag_;
    //! disable adjustStep()
    bool fixStepSize_;
    //! status flag when final destination reached
    bool reachedLastDest_;
    //! status flag giving detect() power over ds
    bool secant_;

    //! perform backtracking when Newton fails
    bool backTracking_;
    //! backtracking counter
    int  backTrack_;
    //! number of backtracking steps
    int  numBackTrackingSteps_;

    //! allowed increase in backtracking
    double backTrackIncrease_;

    //! continuation step counter
    int step_;
    //! maximum number of continuation steps:
    int maxSteps_;

    //! reset counter
    int resetCounter_;

    //! newton iteration counter
    int newtonIter_;
    //! newton iteration counter
    int sumNewtonIter_;
    //! set the max for the Newton iterations
    int maxNewtonIterations_;
    //! set the min for the Newton iterations
    int minNewtonIterations_;
    //! if less/more then we increase/decrease the step-size
    double optNewtonIterations_;
    //! set the tolerance in the Newton corrector
    double newtonTolerance_;
    //! set the tolerance for the final convergence
    double destinationTolerance_;
    //! let model class determine a simple monitor
    bool userDetect_;
    //! specify type of special point detection:
    //!      'D' converge at destination 0,1,... (default)
    //!      'P' converge at turning point
    char detectMode_;
    //! scaling of the state tangent
    double tanScaling_;
    //! normalization strategy 'O' old or 'N' new
    char normalizeStrategy_;
    //! eigenvalue analysis: 'N' never,
    //!                      'E' at the end of a run,
    //!                      'P' at every converged point (during post-processing).
    char eigenvalueAnalysis_;
    //! set to false if you feel lucky
    bool rejectFailedNewton_;
    //! give up at minimum step size
    bool giveUpAtdsMin_;

    //! If the parameter dependence in the rhs is constant we can reuse statedot
    //! in the Newton corrector, partially disregarding updates in the Jacobian.
    //! This means we do a partial Newton-chord iteration.
    bool newtChordHybr_;

    //! Specify the tangent type in the body of the continuation
    //! E: Euler
    //! S: Secant
    char tangentType_;

    //! Specify the residual test for convergence in the Newton process
    //! R: ||R||_2,            2-norm of R (augmented F)
    //! D: ||(dx,dpar)||_inf,  inf-norm of the update (dx,dpar)
    char residualTest_;

    //! E: Euler, solve J dxdpar = -dF
    //! A: Assign  dx = -DF (when J approx I)
    char initialTangent_;

    //! print important vectors for manual inspection
    bool printImportantVectors_;

    //! post processing:
    //!  "at every point"
    //!  "at final point"
    std::string postProcess_;

    //! Specify a bound on the RHS computed from the predicted state
    //! If it exceeds the bound we choose a smaller step size-ds
    double predictorBound_;

    //! used for detecting sign switch
    int parDotSign_;

    //! track the sign of monitor in detect()
    std::vector<int> signMonitor_;
    //! ||F(x_old)||
    double normRHS_;
    //! ||F(x_new)||
    double normRHStest_;

    //! See Store() and Restore() for its use
    struct Storage
    {
        //! prev state
        VectorPtr state0;
        //! prev prev state
        VectorPtr state00;
        //! prev d/ds state
        VectorPtr stateDot0;
        //! prev par
        double par0;
        //! prev prev par
        double par00;
        //! prev ds
        double ds0;
        //! prev prev ds
        double ds00;
        //! prev pardot
        double parDot0;
    };

    Storage storage_;

    std::shared_ptr<JDQZsolver> jdqz_;

public:

    //! default constructor
    Continuation(){};

    //! constructor
    Continuation(Model model);
    Continuation(Model model, Teuchos::RCP<Teuchos::ParameterList> pars);
    Continuation(Model model, Teuchos::ParameterList& pars);

    //! destructor
    ~Continuation();

    //! run continuation
    int run();

    //! test
    void test();

    const Teuchos::ParameterList& getParameters();
    void setParameters(Teuchos::ParameterList&);

    static Teuchos::ParameterList getDefaultInitParameters();
    static Teuchos::ParameterList getDefaultParameters();

private:

    int  step();
    void store();
    void restore();
    void analyzeHist();
    void testCopyView();
    void initialize();

    //! Create initial Euler tangent
    void createInitialTangent();

    //! Modes:  'S' Secant
    //!         'E' Euler
    void createTangent(char mode);

    //! Normalize the tangent with an additional scaling
    //! for the state component.
    void normalize();

    //! Compute the derivative of F with respect to
    //! the continuation parameter. Mode governs the
    //! computation of the RHS in the model.
    //! Modes: 'F' : force compute RHS
    //!        'A' : do not force compute RHS
    void computeDFDPar(char mode = 'A');

    int  eulerPredictor();
    int  newtonCorrector();

    int  runBackTracking(VectorPtr stateDir, double parDir);

    //! Detect special points.
    void detect();

    //! let model::monitor define stopping criterion
    void userDetect();

    void adjustStep();
    void reset();
    void info();
    void finalize();
    void modelInfo();

    //! solve generalized eigenvalue problem
    void eigenSolver();

    //! write essential continuation data to datafile
    void writeData(bool describe = false);
};

#endif
