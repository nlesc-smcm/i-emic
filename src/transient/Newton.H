#ifndef NEWTON_H
#define NEWTON_H

#include <functional>

template<typename Model>
class Newton
{
public:
    using VectorPtr = typename Model::element_type::VectorPtr;
    using ConstVectorPtr = typename Model::element_type::ConstVectorPtr;

private:
    Model model_;

    std::function<ConstVectorPtr(
        ConstVectorPtr const &)> F_;
    std::function<ConstVectorPtr(
        ConstVectorPtr const &,
        ConstVectorPtr const &)> Jsol_;

    double tol_;
    int max_newton_steps_;

    bool converged_;
    int newton_steps_;
    double normdx_;
    double normF_;

    ConstVectorPtr Fx_;
public:
    template<typename ParameterList>
    Newton(Model model, ParameterList params);

    ConstVectorPtr run(
        ConstVectorPtr const &x0);

    bool converged() const;
    double normF() const;
    double normdx() const;
    int steps() const;
    ConstVectorPtr Fx() const;
};

template<typename Model>
template<typename ParameterList>
Newton<Model>::Newton(Model model, ParameterList params)
    :
    model_(model),
    tol_(params->get("Newton tolerance", 1e-8)),
    max_newton_steps_(params->get("maximum Newton iterations", 20))
{
    // Deterministic theta stepper:
    // M * u_n + dt * theta * F(u_(n+1)) + dt * (1-theta) * F(u_n) - M * u_(n+1) = 0
    // Noise is added in an explicit manner
    F_ = [this](ConstVectorPtr const &xnew) {
        TIMER_SCOPE("Newton: F");
        model_->setState(xnew);
        model_->computeRHS();
        return model_->getRHS('V');
    };

    // J2 = theta * dt * J - M, so J2*x = b
    // We write this as
    // J2 = J - 1/(theta*dt) * M, J2 * x = 1/(theta*dt) * b
    Jsol_ = [this](ConstVectorPtr const &xnew,
                    ConstVectorPtr const &b) {
        TIMER_SCOPE("Newton: Jacobian solve");
        model_->setState(xnew);
        model_->computeJacobian();
        model_->solve(b);
        return model_->getSolution('V');
    };
}

template<typename Model>
typename Newton<Model>::ConstVectorPtr Newton<Model>::run(
   ConstVectorPtr const &x0)
{
    TIMER_SCOPE("Newton: Newton");

    // Clone of the initial vector, which we only use here. This is
    // so the vector we return can not be changed when updating
    // the model. This is required in (T)AMS, because we store
    // every state.
    VectorPtr x = Utils::clone(x0);

    Fx_ = F_(x);
    normF_ = -1;
    converged_ = false;

    for (newton_steps_ = 0; newton_steps_ < max_newton_steps_; newton_steps_++)
    {
        ConstVectorPtr dx = Jsol_(x, Fx_);
        normdx_ = Utils::normInf(dx);

        CHECK_ZERO(x->Update(-1.0, *dx, 1.0));
        Fx_ = F_(x);
        normF_ = Utils::norm(Fx_);

        INFO("  Newton solver ------------------------------------");
        INFO("                            iter     = " << newton_steps_);
        INFO("                           ||F||2    = " << normF_);
        INFO("                           ||dx||inf = " << normdx_);
        INFO("\n");

        if (normdx_ < tol_ && normF_ < tol_)
        {
            converged_ = true;
            return x;
        }
        else if (normdx_ > 1e2)
        {
            WARNING("Norm exploding! ||dx||inf = "
                    << normdx_ << "\n", __FILE__, __LINE__);
            break;
        }
    }
    WARNING("Newton did not converge in " << newton_steps_
            << "steps with ||F|| = "
            << normF_ << "\n", __FILE__, __LINE__);
    return x;
}

template<typename Model>
bool Newton<Model>::converged() const
{
    return converged_;
}

template<typename Model>
double Newton<Model>::normF() const
{
    return normF_;
}

template<typename Model>
double Newton<Model>::normdx() const
{
    return normdx_;
}

template<typename Model>
int Newton<Model>::steps() const
{
    return newton_steps_;
}

template<typename Model>
typename Newton<Model>::ConstVectorPtr Newton<Model>::Fx() const
{
    return Fx_;
}

#endif
