/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

#include "iostream"
#include "Epetra_Comm.h"
#include "EpetraExt_MultiComm.h"
#include "globdefs.H"
#include "ThetaStepperEvaluator.H"

#include "EpetraExt_MatrixMatrix.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

//! construct with given model evaluator
ThetaStepperEvaluator::ThetaStepperEvaluator(
	Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> model, 
	Teuchos::RCP<Epetra_CrsMatrix> jacPtr, 
	double t, double theta, double dt)
    :
	model_(model), jacPtr_(jacPtr), dt_(dt), t_(t)
{
    DEBUG("ThetaStepperEvaluator: constructor");
    const Epetra_Map& map = jacPtr->RowMap();
    theta_ = Teuchos::rcp(new Epetra_Vector(map));
    x_dot_ = Teuchos::rcp(new Epetra_Vector(map));
    x_old_ = Teuchos::rcp(new Epetra_Vector(map));
    f_old_ = Teuchos::rcp(new Epetra_Vector(map));
    f_new_ = Teuchos::rcp(new Epetra_Vector(map));        
    reset(t_,*x_old_);
    // evaluate model to get B matrix
    DEBUG("ThetaStepperEvaluator - compute mass matrix...");
    model_->computeShiftedMatrix(0.0,1.0, *x_old_, *jacPtr);
    B_ = Teuchos::rcp(new Epetra_CrsMatrix(*jacPtr));
    diagB_ = Teuchos::rcp(new Epetra_Vector(map));
    CHECK_ZERO(B_->ExtractDiagonalCopy(*diagB_));
    // construct default theta's: use given scalar theta where
    // B is nonzero, and 1 where B is zero.
    this->set_theta(theta);
    DEBUG("ThetaStepperEvaluator: constructor done");
}
     
     
//! destructor
ThetaStepperEvaluator::~ThetaStepperEvaluator()
{
}

void ThetaStepperEvaluator::set_theta(double theta)
{
    for (int i=0;i<theta_->MyLength();i++)
	{
		if ( (*diagB_)[i] == 0.0) (*theta_)[i] = 1.0;
		else                 (*theta_)[i] = theta;
	}
}


// reset stepper: f has to be re-evaluated, x_old and t_ are set
bool ThetaStepperEvaluator::reset(double t_n, const Epetra_Vector &x_n)
{
    t_      = t_n;
    *x_old_ = x_n;

    x_dot_->PutScalar(0.0);

    model_->setXdot(*x_dot_, t_);
	
    bool success = model_->computeF(*x_old_, *f_old_, NOX::Epetra::Interface::Required::Residual);

    return success;
}   


bool ThetaStepperEvaluator::computeF(const Epetra_Vector& x, Epetra_Vector& F,
									 const FillType fillFlag)
{
	bool stat = true;
    
    x_dot_->PutScalar(0.0);
    model_->setXdot(*x_dot_, t_+dt_);
	
	stat = model_->computeF(x, *f_new_, fillFlag);
	
    // x_dot = (x - x_old)/dt
    CHECK_ZERO(x_dot_->Update(1.0/dt_, x, -1.0/dt_, *x_old_, 0.0));
    CHECK_ZERO(B_->Multiply(false, *x_dot_, F));
    
    // set F = 1/dtB(u-u_old) + (1-theta)f(u_old) + theta f(u)
    for (int i = 0; i < F.MyLength(); i++)
	{
		F[i] += (1 - (*theta_)[i]) * (*f_old_)[i] 
			    + (*theta_)[i]*(*f_new_)[i];
	}

    return stat;                 
}

bool ThetaStepperEvaluator::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
    x_dot_->PutScalar(0.0);
    // set time (x_dot doesn't matter for matrix)
    model_->setXdot(*x_dot_, t_+dt_);
    
    // instead we compute A, 'left scale' it by theta_, add 1/dt*B.
    
    // we put the matrix into model->getJacobian, that should be the same
    // as Jac, in general.
    Teuchos::RCP<Epetra_CrsMatrix> A = jacPtr_;
    if (&Jac!=A.get()) Error("expecting shared pointer",__FILE__,__LINE__);
    bool success = model_->computeShiftedMatrix(1.0,0.0,x,Jac);
    const Epetra_Vector& alpha = *theta_;
    double beta = 1/dt_;
    CHECK_ZERO(A->LeftScale(alpha));
    CHECK_ZERO(EpetraExt::MatrixMatrix::Add(*B_,false,beta,*A,1.0));
    return success;
}

