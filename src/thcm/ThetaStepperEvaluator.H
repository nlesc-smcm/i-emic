/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef THETA_STEPPER_EVALUATOR_H
#define THETA_STEPPER_EVALUATOR_H

#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "LOCA_Epetra_Interface_TimeDependent.H"
#include "Epetra_Vector.h"

//! this class defines the nonlinear system to be solved by 
//! Newton's method if you want to do 'Theta' type time-    
//! integration, i.e. Theta=1:   backward Euler               
//!                   Theta=0:   forward Euler                
//!                   Theta=0.5: Crank-Nicholson            
//!                   etc.                                  
class ThetaStepperEvaluator : public NOX::Epetra::Interface::Required,
							  public NOX::Epetra::Interface::Jacobian
  {

public:
  
  //! construct with given model evaluator
  ThetaStepperEvaluator(Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> model, 
         Teuchos::RCP<Epetra_CrsMatrix> jacPtr, double t_=0.0, double theta=0.5, double dt=1.0);
  
  //! destructor
  virtual ~ThetaStepperEvaluator();
  
  //! \name get/set functions
  
  //!@{
  
  //! theta will be set for matrix rows where B!=0.
  //! if B==0, the local theta will be 1.
  void set_theta(double theta);

  //! we allow using a different theta for every variable
  void set_theta(const Epetra_Vector& theta){*theta_=theta;}

  //! //==??
  void set_dt(double dt){dt_=dt;}
  
  //! get f(x_n,t_n) after setting x_n by caling reset()
  const Epetra_Vector& get_f(){return *f_old_;}
  
  //!@}
  
  //! reset stepper: this should be called before the Newton
  //! solve to make sure the correct x_old is used.         
  bool reset(double t_n, const Epetra_Vector& x_n);
            
  //! \name NOX Required and Jacobian interface
  
  //!@{

  //! compute F(u) = (1/dt)B (u-u_n) + Theta*f(u) + (1-Theta) f(u_n)
  //! TODO: might it be better to use B udot + f(1/2(u+u_n)) ?
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F,
                            const FillType fillFlag);
  //!
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //!@}

protected:

  //! pointer to model
  Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> model_;
  
  //! jacobian operator
  Teuchos::RCP<Epetra_CrsMatrix> jacPtr_;

  //! pointer to mass matrix
  Teuchos::RCP<Epetra_CrsMatrix> B_;
  
  //! diagonal of B in a vector;
  Teuchos::RCP<Epetra_Vector> diagB_;

  //! theta parameter
  Teuchos::RCP<Epetra_Vector> theta_;
  
  //! time level t
  double t_;
  
  //! step size for next step
  double dt_;
  
  //! old solution vector x_n
  Teuchos::RCP<Epetra_Vector> x_old_;
  
  //! approximation of derivative, x_new - x_old
  Teuchos::RCP<Epetra_Vector> x_dot_;
  
  //! f(u_n), last evaluation of f
  Teuchos::RCP<Epetra_Vector> f_old_, f_new_;

  };

#endif
