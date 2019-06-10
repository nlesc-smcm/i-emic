#ifndef STOCHASTICPROJECTEDTHETAMODEL_H
#define STOCHASTICPROJECTEDTHETAMODEL_H

#include "ProjectedThetaModel.H"

#include <Teuchos_RCP.hpp>

#include "Epetra_MultiVector.h"

template<typename Model>
class StochasticProjectedThetaModel : public ProjectedThetaModel<Model>
{
public:
    using VectorPtr = typename Model::VectorPtr;

protected:
    //! Noise parameter
    double sigma_;

    //! Stochastic forcing
    Teuchos::RCP<Epetra_CrsMatrix> B_;

    //! Projected stochastic forcing
    Teuchos::RCP<Epetra_MultiVector> BV_;

    //! Noise at the current timestep
    Teuchos::RCP<Epetra_MultiVector> G_;

    //! RNG
    unsigned int noise_seed_;
    Teuchos::RCP<std::mt19937_64> engine_;

public:
    //-------------------------------------------------------
    //! constructor
    template<typename ParameterList>
    StochasticProjectedThetaModel(Teuchos::RCP<Epetra_Comm> comm, ParameterList params,
                                  Teuchos::RCP<Epetra_MultiVector> const &V)
        :
        ProjectedThetaModel<Model>(comm, params, V),
        sigma_(params->get("sigma", 1.0)),
        noise_seed_(params->get("noise seed", 0))
        {
            std::random_device::result_type seed = noise_seed_;
            if (seed == 0)
            {
                std::random_device rd;
                seed = rd();
            }
            std::seed_seq seeder{seed};
            engine_ = Teuchos::rcp(new std::mt19937_64(seeder));

            Model::computeForcing();
            B_ = Model::getForcing();

            BV_ = Teuchos::rcp(
                new Epetra_MultiVector(
                    B_->DomainMap(),
                    ProjectedThetaModel<Model>::V_->NumVectors()));
            CHECK_ZERO(B_->Multiply(
                           'T', *ProjectedThetaModel<Model>::V_, *BV_));
        }

    template<typename ParameterList>
    StochasticProjectedThetaModel(Model const &model, ParameterList params,
                                  Teuchos::RCP<Epetra_MultiVector> const &V)
        :
        ProjectedThetaModel<Model>(model, params, V),
        sigma_(params->get("sigma", 1.0)),
        noise_seed_(params->get("noise seed", 0))
        {
            std::random_device::result_type seed = noise_seed_;
            if (seed == 0)
            {
                std::random_device rd;
                seed = rd();
            }
            std::seed_seq seeder{seed};
            engine_ = Teuchos::rcp(new std::mt19937_64(seeder));

            Model::computeForcing();
            B_ = Model::getForcing();

            BV_ = Teuchos::rcp(
                new Epetra_MultiVector(
                    B_->DomainMap(),
                    ProjectedThetaModel<Model>::V_->NumVectors()));
            CHECK_ZERO(B_->Multiply(
                           'T', *ProjectedThetaModel<Model>::V_, *BV_));
        }

    void initStep(double timestep)
        {
            ProjectedThetaModel<Model>::initStep(timestep);

            // Compute noise for forcing which is independent per processor
            std::normal_distribution<double> distribution(0.0, 1.0);
            auto generator = std::bind(distribution, std::ref(*engine_));

            int m = BV_->MyLength();
            if (!BV_->Map().UniqueGIDs())
            {
                ERROR("The values of B are distributed", __FILE__, __LINE__);
            }

            Epetra_Vector pert(BV_->Map());
            std::generate(pert.Values(), pert.Values() + m, generator);

            Epetra_LocalMap Gmap(BV_->NumVectors(), 0, BV_->Map().Comm());
            G_ = Teuchos::rcp(new Epetra_MultiVector(Gmap, 1));
            CHECK_ZERO(G_->Multiply('T', 'N', 1.0, *BV_, pert, 0.0));
            CHECK_ZERO(G_->Scale(sqrt(ThetaModel<Model>::timestep_) * sigma_));
        }

    void setState(VectorPtr state)
        {
            ProjectedThetaModel<Model>::setState(state);
        }

    //!-------------------------------------------------------
    //! Compute theta method RHS
    //!
    //! The final theta timestepping RHS is given by
    //! M * u_n + dt * theta * F(u_(n+1)) + dt * (1-theta) * F(u_n) - M * u_(n+1) + G * dW = 0
    void computeRHS()
        {
            ProjectedThetaModel<Model>::computeRHS();
            CHECK_ZERO(Model::rhs_->Update(1.0, *G_, 1.0));
        }

    void computeJacobian()
        {}

    //!-------------------------------------------------------
    //! J2 * x = 1/(theta*dt) * b
    void solve(Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null)
        {
            ProjectedThetaModel<Model>::solve(rhs);
        }
};

#endif
