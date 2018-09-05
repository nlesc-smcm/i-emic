#ifndef THETA_H
#define THETA_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "GlobalDefinitions.H"

//! Here we inherit a templated model and adjust the rhs and jac
//! computation to create a time (theta) stepping problem.

template<typename Model>
class Theta : public Model
{
public:
    using VectorPtr     = typename Model::VectorPtr;
    using MatrixPtr     = typename Model::MatrixPtr;
    using ParameterList = typename Model::ParameterList;

private: 
    //! theta method parameter
    double theta_;

    //! time step size
    double timestep_;

    VectorPtr oldState_;
	VectorPtr xDot_;
	VectorPtr Bxdot_;
	VectorPtr oldRhs_;

public:
    //-------------------------------------------------------
    //! constructor
    Theta(Teuchos::RCP<Epetra_Comm> comm, ParameterList params)
        :
        Model(comm, params),
        theta_     (1.0),
        timestep_  (1.0e-3)
        {
            // Initialize a few datamembers
            oldState_ = Model::getState('C');
            xDot_     = Model::getState('C');
            Bxdot_    = Model::getState('C');
            oldRhs_   = Model::getRHS('C');
        }

    //!-------------------------------------------------------
    void store()
        {
            if (oldState_ == Teuchos::null or
                !oldState_->Map().SameAs(Model::state_->Map()))
            {
                oldState_ =
                    Teuchos::rcp(new Epetra_Vector(Model::state_->Map()));
            }
	
            *oldState_ = *Model::state_;

            if (oldRhs_ == Teuchos::null or
                !oldRhs_->Map().SameAs(Model::rhs_->Map()))
            {
                oldRhs_ =
                    rcp(new Epetra_Vector(Model::rhs_->Map()));
            }
	
            *oldRhs_   = *Model::rhs_;
        }

    //!-------------------------------------------------------
	void restore()
        {
            if (Model::state_ == Teuchos::null or
                !Model::state_->Map().SameAs(oldState_->Map()))
            {
                Model::state_ =
                    Teuchos::rcp(new Epetra_Vector(oldState_->Map()));
            }
	
            *Model::state_ = *oldState_;
	
            if (Model::rhs_ == Teuchos::null or
                !Model::rhs_->Map().SameAs(oldRhs_->Map()))
            {
                Model::rhs_ =
                    Teuchos::rcp(new Epetra_Vector(oldRhs_->Map()));
            }
	
            *Model::rhs_ = *oldRhs_;

        }

    //!-------------------------------------------------------
	void setTimestep(double timestep){ timestep_ = timestep; }

    //!-------------------------------------------------------
	void setTheta(double theta)
        {
            // Check theta
            if (theta <= 0 || theta > 1)
            {
                ERROR("Incorrect theta:" << theta,
                      __FILE__, __LINE__);
            }
            theta_ = theta;        
        }

    //!-------------------------------------------------------
    //! Compute theta method RHS
    //!
    //! The final theta timestepping RHS is given by
	//! (B dx/dt)/theta + F(x) + (theta-1)/theta * F(x_old)
    
	void computeRHS()
        {
            // Check theta
            if (theta_ <= 0 || theta_ > 1)
            {
                WARNING("Theta: Incorrect theta: " << theta_,
                        __FILE__, __LINE__);
            }	

            // compute ordinary discretization
            Model::computeRHS();

            // compute mass matrix
            Model::computeMassMat();

            // calculate xdot = (xnew - xold) / (dt*theta)
            double s = 1.0 / timestep_ / theta_;
            xDot_->Update(s , *Model::state_, -s ,
                              *oldState_, 0.0);

            // Scale xdot with the values in the mass matrix
            Model::applyMassMat(*xDot_, *Bxdot_);

            // Combine terms
            double f = (theta_ - 1)/theta_;
            Model::rhs_->Update(f,    *oldRhs_, 1.0);
            Model::rhs_->Update(1.0,  *Bxdot_, 1.0);
            // Model::rhs_->Scale(-1.0); // why            
        }

    //!-------------------------------------------------------
    //! compute derivative of theta method rhs:
    //!  B / (theta * dt) + J.
    
	void computeJacobian()
        {

            Model::computeMassMat();
            VectorPtr diagB = Model::getMassMat();
            
            // Compute the ordinary Jacobian using the current state
            Model::computeJacobian();
            
            // Get the number of local elements in the matrix
            int numMyElements =	Model::jac_->Map().NumMyElements();
            
            // Get a list of the global element IDs owned by the calling proc
            int *myGlobalElements = Model::jac_->Map().MyGlobalElements();

            // Add to the Jacobian the values B[i] on the diagonal.
            // The result is B / (theta * dt) + J.
            double value;
            for (int i = 0; i != numMyElements; ++i)
            {
                value = (*diagB)[i] / timestep_ / theta_;
                Model::jac_->SumIntoGlobalValues(myGlobalElements[i], 1,
                                                 &value, myGlobalElements + i);
            }
            Model::jac_->FillComplete();
        }
};

#endif
