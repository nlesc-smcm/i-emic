/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef IMPLICIT_TIME_STEPPER_H
#define IMPLICIT_TIME_STEPPER_H

#include "AbstractTimeStepper.H"

class OceanModel;
class ThetaStepperEvaluator;
namespace NOX {
	namespace Epetra {
		class Group;
		class LinearSystem;
	}
}

//! class implementing Theta-method time stepping (Crank-Nicholson etc.)

/*! TODO: finish documentation
  This class is currently not very general and does not conform to the 
  steppers from the Rythmos package.
*/
class ImplicitTimeStepper : public AbstractTimeStepper
{
  
public: 
    
    //! Construct from OceanModel and NOX data structures
    ImplicitTimeStepper(Teuchos::RCP<OceanModel> model,
                        Teuchos::ParameterList& stepperParams,
                        Teuchos::ParameterList& nlParams);

    //! Destructor
    virtual ~ImplicitTimeStepper();
    
    //! reset (call this when a time step failed)
    void reset();
    
    //! advance one step
    bool Step(const Epetra_Vector &x_old, double t_old, 
			        Epetra_Vector &x_new, double dt);

	//! returns true if an estimate of the local truncation error is available
   
	/*! the estimate is computed only after the second step after a reset,
	  so before adapting the timestep etc. one should check if the estimate
	  is actually available.
	*/
	bool hasErrorEstimate(void) const;
   
	//! return an estimate of the 1-norm of the LTE or -1 if not available
	double getErrorEstimate(void) const {return d_lteEst;}

	//! return an estimate of difficulty for the last step performed
	//! (0: easy, 1.0: at the limits, > 1.0: failed)
	double getEffortEstimate(void) const {return d_effortEst;}
   
	//! stability condition: none needed
	double getStabCond(void) const {return -1.0;}

	//! returns the F(u) where u is the most recent approximation to the solution
   
	/*! when the last time step u_{n+1} failed to converge, this is F(u_n) instead.
	 */
	const Epetra_Vector& getF(void) {return *d_f_n;}
    
protected:
   
	//! The ocean model. 
	Teuchos::RCP<OceanModel> d_model;  

	//! Input parameters.
	Teuchos::ParameterList& d_paramList;

	//! Nonlinear solver parameters.
	Teuchos::ParameterList& d_nlParams;
        
    //! The Theta Model (computes F and J for Newton)
	//! is a ThetaStepperEvaluator object.
    Teuchos::RCP<ThetaStepperEvaluator> d_thetaModel;
    
    //! Linear system.
    Teuchos::RCP<NOX::Epetra::LinearSystem> d_linSys;
    
    //! the NOX group (to set starting guess etc.))
    Teuchos::RCP<NOX::Epetra::Group> d_curGroup;

    //! convergence tests for NOX solver
    Teuchos::RCP<NOX::StatusTest::Generic> d_convTests;
    
    //! the nonlinear solver
    Teuchos::RCP<NOX::Solver::Generic> d_nlSolver;

    //! std::string identifying the scheme to be used
    
    /*! supported values:   
	  - "Theta"     
	  - "Adaptive Theta"  
    */
    std::string d_scheme;

    //! output filestream
    Teuchos::RCP<std::ostream> d_out;    

    //! \name parameters for the theta method
    //@{ 
    
    //! weighting factor between forward and backward Euler
    double d_theta;
    
    //! predictor scheme ("Constant" or "Secant")
    std::string d_predictor;
    
    //! estimate of the local truncation error (1-norm)
    double d_lteEst;

    //! estimate of the effort required for last step performed
    //! (0: easy, 1: hard)
    double d_effortEst;
    
    //! do we compute LTE estimates?
    bool d_haveNormEst;

    //! preconditioner reuse: reuse for how many _timesteps_?
    int d_maxPrecAge, d_curPrecAge;
    
    //! is true if the last step failed. Then we need to compute f(u_n) at the 
    //! start of the next step and will not have an error estimate for u_{n+1}.
    bool d_restarted;
    
    //@}
    
    
private: 

    //! after computing u^{n+1} from u^n, we approximate the time-derivative
    //! by a secant, u' = (u^{n+1}-u^n)/dt_n. This will then be used in the 
    //! next time-step as a predictor before the Newton solve. We use the   
    //! secant predictor ustart = u_n + dt * uprime
    Teuchos::RCP<Epetra_Vector> d_uprime, d_ustart;
    
    //! t-derivative at previous time step, approximated   
    //! by F(u_{n-1}), and at the current/next time level. 
    //! Note that these are computed beforehand, i.e.      
    //! at the end of time-step n-1 (at t=t_n) f_n is      
    //! computed and can be used for the error estimate    
    //! in step n (at t_n+1)
    Teuchos::RCP<Epetra_Vector> d_f_nm1, d_f_n, d_f_np1;
    
    //! previous time step-size \delta_t_{n-1} = t_n-t_n-1,
    //! and current time step-size \delta_t_n=t_n+1 - t_n
    double d_dt_nm1, d_dt_n;
    
    //! estimate the local truncation error
    
    /*! we estimate the LTE using finite differences to approximate
	  the second and third derivative of u_n
    
	  \TODO: finish documentation
    */

	//! See legacy code for old implementation without adaptive stepsizes.
	//! Current version uses proper finite differences for variable dt.
    double estimateError(double theta,
						 const Epetra_Vector& rh0,
						 const Epetra_Vector& rh1,
						 const Epetra_Vector& rh2,
						 double dt_old, double dt) const;
	
};

#endif
