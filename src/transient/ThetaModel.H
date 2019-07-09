#ifndef THETAMODEL_H
#define THETAMODEL_H

#include "GlobalDefinitions.H"

//! Here we inherit a templated model and adjust the rhs and jac
//! computation to create a time (theta) stepping problem.

template<typename Model>
class ThetaModel : public Model
{
public:
    using VectorPtr      = typename Model::VectorPtr;
    using ConstVectorPtr = typename Model::ConstVectorPtr;

protected:
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
    template<typename ParameterList>
    ThetaModel(Teuchos::RCP<Epetra_Comm> comm,
               ParameterList model_params,
               ParameterList params)
        :
        Model(comm, model_params),
        theta_(params->get("theta", 1.0)),
        timestep_(1.0e-3)
        {
            // Initialize a few datamembers
            oldState_ = Model::getState('C');
            xDot_     = Model::getState('C');
            Bxdot_    = Model::getState('C');
            oldRhs_   = Model::getRHS('C');
        }

    template<typename ParameterList>
    ThetaModel(Model const &model, ParameterList params)
        :
        Model(model),
        theta_(params->get("theta", 1.0)),
        timestep_(1.0e-3)
        {
            // Initialize a few datamembers
            oldState_ = Model::getState('C');
            xDot_     = Model::getState('C');
            Bxdot_    = Model::getState('C');
            oldRhs_   = Model::getRHS('C');
        }

    virtual ~ThetaModel() {}

    //!-------------------------------------------------------
    virtual void initStep(double timestep)
        {
            timestep_ = timestep;

            *oldState_ = *Model::getState('V');

            Model::preProcess();

            Model::computeRHS();
            *oldRhs_ = *Model::getRHS('V');
        }

    virtual void setState(ConstVectorPtr state)
        {
            if (Model::getState('V') != state)
                *Model::getState('V') = *state;
        }

    //!-------------------------------------------------------
    //! Compute theta method RHS
    //!
    //! The final theta timestepping RHS is given by
    //! M * u_n + dt * theta * F(u_(n+1)) + dt * (1-theta) * F(u_n) - M * u_(n+1) = 0
    virtual void computeRHS()
        {
            // Check theta
            if (theta_ < 0 || theta_ > 1)
            {
                WARNING("ThetaModel: Incorrect theta: " << theta_,
                        __FILE__, __LINE__);
            }

            // Compute ordinary discretization
            Model::computeRHS();

            // Compute mass matrix
            Model::computeMassMat();

            // Compute M * u_n - M * u_(n+1)
            CHECK_ZERO(xDot_->Update(-1.0, *Model::getState('V'), 1.0,
                                     *oldState_, 0.0));
            Model::applyMassMat(*xDot_, *Bxdot_);

            // Compute dt * theta * F(u_(n+1)) + dt * (1-theta) * F(u_n)
            CHECK_ZERO(Model::getRHS('V')->Update(
                           timestep_ * (1-theta_), *oldRhs_, timestep_ * theta_));

            // Combine the terms
            CHECK_ZERO(Model::getRHS('V')->Update(1.0,  *Bxdot_, 1.0));
        }

    //!-------------------------------------------------------
    //! compute derivative of theta method rhs:
    //! J2 = J - 1/(theta*dt) * M
    virtual void computeJacobian()
        {
            // Compute the ordinary Jacobian using the current state
            Model::computeJacobian();

            if (theta_ == 0)
                return;

            Model::computeMassMat();
            ConstVectorPtr diagB = Model::getMassMat('V');

            // Get the number of local elements in the matrix
            int numMyElements = Model::getJacobian()->Map().NumMyElements();

            // Get a list of the global element IDs owned by the calling proc
            int *myGlobalElements = Model::getJacobian()->Map().MyGlobalElements();

            // Add to the Jacobian the values B[i] on the diagonal.
            // The result is B / (theta * dt) + J.
            double value;
            for (int i = 0; i != numMyElements; ++i)
            {
                value = -(*diagB)[i] / timestep_ / theta_;
                CHECK_ZERO(Model::getJacobian()->SumIntoGlobalValues(
                               myGlobalElements[i], 1,
                               &value, myGlobalElements + i));
            }
            CHECK_ZERO(Model::getJacobian()->FillComplete());
        }

    //!-------------------------------------------------------
    //! J2 * x = 1/(theta*dt) * b
    virtual void solve(ConstVectorPtr rhs)
        {
            VectorPtr b = Model::getRHS('C');

            if (theta_ == 0.0)
            {
                ConstVectorPtr M = Model::getMassMat('V');
                for (int i = 0; i < b->MyLength(); i++)
                    (*Model::getSolution('V'))[i] = -(*b)[i] / (*M)[i];
                return;
            }

            CHECK_ZERO(b->Scale(1.0 / timestep_ / theta_));
            Model::solve(b);
        }
};

#endif
