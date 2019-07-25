#ifndef PROJECTEDTHETAMODEL_H
#define PROJECTEDTHETAMODEL_H

#include "ThetaModel.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"

#include "GlobalDefinitions.H"

Teuchos::RCP<Epetra_MultiVector> dot(
    Epetra_MultiVector const &x, Epetra_MultiVector const &y);

//! Here we inherit a templated model and adjust the rhs and jac
//! computation to create a time (theta) stepping problem.

template<typename Model>
class ProjectedThetaModel : public ThetaModel<Model>
{
protected:
    Teuchos::RCP<Epetra_Vector> smallState_;
    Teuchos::RCP<Epetra_Vector> largeRhs_;

    Teuchos::RCP<const Epetra_MultiVector> V_;
    Teuchos::RCP<Epetra_MultiVector> VMV_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> VAVmat_;
    Teuchos::RCP<Epetra_SerialDenseSolver> VAVsolver_;

public:
    //-------------------------------------------------------
    //! constructor
    template<typename ParameterList>
    ProjectedThetaModel(Teuchos::RCP<Epetra_Comm> comm,
                        ParameterList model_params,
                        ParameterList params,
                        Teuchos::RCP<const Epetra_MultiVector> const &V)
        :
        ThetaModel<Model>(comm, model_params, params),
        V_(V)
        {
            // Initialize a few datamembers
            smallState_ = restrict(*Model::getState('V'));
            largeRhs_ = Model::getRHS('C');
            ThetaModel<Model>::xDot_ = restrict(*Model::getState('V'));
            ThetaModel<Model>::Bxdot_ = restrict(*Model::getState('V'));
        }

    template<typename ParameterList>
    ProjectedThetaModel(Model const &model, ParameterList params,
                        Teuchos::RCP<const Epetra_MultiVector> const &V)
        :
        ThetaModel<Model>(model, params),
        V_(V)
        {
            // Initialize a few datamembers
            smallState_ = restrict(*Model::getState('V'));
            largeRhs_   = Model::getRHS('C');
            ThetaModel<Model>::xDot_ = restrict(*Model::getState('V'));
            ThetaModel<Model>::Bxdot_ = restrict(*Model::getState('V'));
        }

    virtual ~ProjectedThetaModel() {}

    //!-------------------------------------------------------
    virtual void initStep(double timestep)
        {
            ThetaModel<Model>::timestep_ = timestep;

            ThetaModel<Model>::oldState_ = restrict(*Model::state_);

            Model::preProcess();

            Model::rhs_ = largeRhs_;
            Model::computeRHS();
            ThetaModel<Model>::oldRhs_ = restrict(*Model::rhs_);

            TIMER_START("ProjectedThetaModel: Compute VMV");
            // Compute mass matrix
            Model::computeMassMat();

            auto tmp = Teuchos::rcp(new Epetra_MultiVector(*V_));
            Model::applyMassMat(*V_, *tmp);
            VMV_ = dot(*V_, *tmp);
            TIMER_STOP("ProjectedThetaModel: Compute VMV");

            TIMER_START("ProjectedThetaModel: Compute VAV");
            if (ThetaModel<Model>::theta_ == 0)
            {
                VAVmat_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                           Copy, VMV_->Values(), VMV_->Stride(),
                                           VMV_->MyLength(), VMV_->NumVectors()));
                // std::cout << *VAVmat_;
            }
            else
            {
                ThetaModel<Model>::computeJacobian();
                Model::applyMatrix(*V_, *tmp);
                auto VAV = dot(*V_, *tmp);
                VAVmat_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                           Copy, VAV->Values(), VAV->Stride(),
                                           VAV->MyLength(), VAV->NumVectors()));
            }
            VAVsolver_ = Teuchos::rcp(new Epetra_SerialDenseSolver());
            CHECK_ZERO(VAVsolver_->SetMatrix(*VAVmat_));
            CHECK_ZERO(VAVsolver_->Factor());
            TIMER_STOP("ProjectedThetaModel: Compute VAV");
        }

    virtual void setState(Teuchos::RCP<const Epetra_Vector> state)
        {
            Model::state_ = prolongate(*state);
            *smallState_ = *state;
        }

    //!-------------------------------------------------------
    //! Compute theta method RHS
    //!
    //! The final theta timestepping RHS is given by
    //! M * u_n + dt * theta * F(u_(n+1)) + dt * (1-theta) * F(u_n) - M * u_(n+1) = 0
    virtual void computeRHS()
        {
            // Check theta
            if (ThetaModel<Model>::theta_ < 0 || ThetaModel<Model>::theta_ > 1)
            {
                WARNING("ProjectedThetaModel: Incorrect theta: " << ThetaModel<Model>::theta_,
                        __FILE__, __LINE__);
            }

            Model::rhs_ = largeRhs_;

            // Compute ordinary discretization
            Model::computeRHS();

            Model::rhs_ = restrict(*Model::rhs_);

            // Compute M * u_n - M * u_(n+1)
            CHECK_ZERO(ThetaModel<Model>::xDot_->Update(
                           -1.0, *smallState_, 1.0, *ThetaModel<Model>::oldState_, 0.0));
            CHECK_ZERO(ThetaModel<Model>::Bxdot_->Multiply(
                           'N', 'N', 1.0, *VMV_, *ThetaModel<Model>::xDot_, 0.0));

            // Compute dt * theta * F(u_(n+1)) + dt * (1-theta) * F(u_n)
            CHECK_ZERO(Model::rhs_->Update(
                           ThetaModel<Model>::timestep_ * (1-ThetaModel<Model>::theta_),
                           *ThetaModel<Model>::oldRhs_,
                           ThetaModel<Model>::timestep_ * ThetaModel<Model>::theta_));

            // Combine the terms
            CHECK_ZERO(Model::rhs_->Update(1.0, *ThetaModel<Model>::Bxdot_, 1.0));
        }

    //!-------------------------------------------------------
    virtual void computeJacobian()
        {}

    //!-------------------------------------------------------
    //! J2 * x = 1/(theta*dt) * b
    virtual void solve(Teuchos::RCP<const Epetra_Vector> rhs)
        {
            TIMER_SCOPE("ProjectedThetaModel: Projected jacobian solve");
            auto y = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                      Copy, rhs->Values(), rhs->Stride(),
                                      rhs->MyLength(), rhs->NumVectors()));
            auto z = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                      Copy, rhs->Values(), rhs->Stride(),
                                      rhs->MyLength(), rhs->NumVectors()));

            if (ThetaModel<Model>::theta_ != 0.0)
            {
                CHECK_ZERO(z->Scale(1.0 / ThetaModel<Model>::timestep_ /
                                    ThetaModel<Model>::theta_));
            }
            else
            {
                CHECK_ZERO(z->Scale(-1.0));
            }
            // std::cout << *y << *z;

            CHECK_ZERO(VAVsolver_->SetVectors(*y, *z));
            CHECK_NONNEG(VAVsolver_->Solve());

            // std::cout << *y;

            Model::sol_ = Teuchos::rcp(
                new Epetra_Vector(
                    Copy, Epetra_LocalMap(y->M(), 0, *Model::Comm()),
                    y->A()));
        }

    //!-------------------------------------------------------
    Teuchos::RCP<Epetra_Vector> restrict(Epetra_MultiVector const &x) const
        {
            TIMER_SCOPE("ProjectedThetaModel: restrict");
            int m = V_->NumVectors();

            Epetra_LocalMap map(m, 0, x.Comm());
            Teuchos::RCP<Epetra_Vector> out = Teuchos::rcp(new Epetra_Vector(map));

            CHECK_ZERO(out->Multiply('T', 'N', 1.0, *V_, x, 0.0));
            return out;
        }

    //!-------------------------------------------------------
    Teuchos::RCP<Epetra_Vector> prolongate(Epetra_MultiVector const &x) const
        {
            TIMER_SCOPE("ProjectedThetaModel: prolongate");
            Teuchos::RCP<Epetra_Vector> out = Teuchos::rcp(new Epetra_Vector(V_->Map()));
            CHECK_ZERO(out->Multiply('N', 'N', 1.0, *V_, x, 0.0));
            return out;
        }

};
#endif
