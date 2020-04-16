#ifndef STOCHASTICTHETAMODEL_H
#define STOCHASTICTHETAMODEL_H

#include "ThetaModel.H"
#include "StochasticBase.H"

#include <Teuchos_RCP.hpp>

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

template<typename Model>
class StochasticThetaModel : public ThetaModel<Model>,
                             public StochasticBase
{
protected:
    //! Noise parameter
    double sigma_;

    //! Stochastic forcing
    Teuchos::RCP<Epetra_CrsMatrix> B_;

    //! Noise at the current timestep
    Teuchos::RCP<Epetra_Vector> G_;

public:
    //-------------------------------------------------------
    //! constructor
    template<typename ParameterList>
    StochasticThetaModel(Teuchos::RCP<Epetra_Comm> comm,
                         ParameterList model_params,
                         ParameterList params)
        :
        ThetaModel<Model>(comm, model_params, params),
        StochasticBase(*comm, params),
        sigma_(params->get("sigma", 1.0))
        {
            Model::computeForcing();
            B_ = Model::getForcing();
        }

    template<typename ParameterList>
    StochasticThetaModel(Model const &model, ParameterList params)
        :
        ThetaModel<Model>(model, params),
        StochasticBase(*model.Comm(), params),
        sigma_(params->get("sigma", 1.0))
        {
            Model::computeForcing();
            B_ = Model::getForcing();
        }

    void initStep(double timestep)
        {
            ThetaModel<Model>::initStep(timestep);

            // Compute noise for forcing which is independent per processor
            std::normal_distribution<double> distribution(0.0, 1.0);
            auto generator = std::bind(distribution, std::ref(*engine_));

            int m = B_->NumMyCols();
            if (!B_->ColMap().UniqueGIDs())
            {
                ERROR("The values of B are distributed", __FILE__, __LINE__);
            }

            Epetra_Vector pert(B_->ColMap());
            std::generate(pert.Values(), pert.Values() + m, generator);

            G_ = Model::getState('C');
            CHECK_ZERO(B_->Apply(pert, *G_));
            CHECK_ZERO(G_->Scale(sqrt(ThetaModel<Model>::timestep_) * sigma_));
        }

    //!-------------------------------------------------------
    //! Compute theta method RHS
    //!
    //! The final theta timestepping RHS is given by
    //! M * u_n + dt * theta * F(u_(n+1)) + dt * (1-theta) * F(u_n) - M * u_(n+1) + G * dW = 0
    void computeRHS()
        {
            ThetaModel<Model>::computeRHS();
            CHECK_ZERO(Model::rhs_->Update(1.0, *G_, 1.0));
        }
};

#endif
