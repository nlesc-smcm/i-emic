/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef ABSTRACT_TIME_STEPPER_H
#define ABSTRACT_TIME_STEPPER_H

#include "LOCA_Epetra_Interface_TimeDependent.H"

class Epetra_Vector;

//! Abstract base class for time integration methods.

/*! TODO: finish documentation
*/
class AbstractTimeStepper 
  {
    
  public: 
    
    //! Destructor
    virtual ~AbstractTimeStepper(){;}
    
    //! reset (call this when a time step failed)
    virtual void reset() = 0;
    
    //! advance one step
    virtual bool Step(const Epetra_Vector& x_old, double t_old,
                            Epetra_Vector& x_new, double dt) = 0;

   //! returns true if an estimate of the local truncation error is available
   
   /*! the estimate is computed only after the second step after a reset,
       so before adapting the timestep etc. one should check if the estimate
       is actually available.
   */
   virtual bool hasErrorEstimate(void) const = 0;
   
   //! return an estimate of the 1-norm of the LTE or -1 if not available
   virtual double getErrorEstimate(void) const = 0;

   //! return an estimate of the effort made in last step performed.
   //! 0: really easy, 1: really hard, >1 step failed.
   //! the default implementation always returns 0.5.
   virtual double getEffortEstimate(void) const {return 0.5;};
	  
   //! get a bound on the time step. If such a bound is not available, -1 is returned
   virtual double getStabCond(void) const = 0;

   //! returns F(u) where u is the most recent approximation to the solution
   
   /*! when the last time step u_{n+1} failed to converge, this is F(u_n) instead.
   */
   virtual const Epetra_Vector& getF(void) = 0;
	  
  };

#endif
