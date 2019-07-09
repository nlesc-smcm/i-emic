#ifndef ADAPTIVETRANSIENT_H
#define ADAPTIVETRANSIENT_H

#include "Teuchos_RCP.hpp"

#include "Newton.H"
#include "Transient.hpp"

template<typename Model>
class AdaptiveTransient: public Transient<typename Model::element_type::ConstVectorPtr>
{
    using ConstVectorPtr = typename Model::element_type::ConstVectorPtr;

    Model model_;

    Teuchos::RCP<Newton<Model> > newton_;

    double tol_;
    int max_newton_steps_;
    bool adaptive_;

    int min_wanted_newton_steps_;
    int max_wanted_newton_steps_;

    double min_time_step_;
    double max_time_step_;

    double time_step_increase_;
    double time_step_decrease_;

    int nsteps_;
    int output_;

    double time_;
    int total_newton_steps_;
    bool init_wd_;

public:
    template<typename ParameterList>
    AdaptiveTransient(Model model, ParameterList params);
    template<typename ParameterList>
    AdaptiveTransient(Model model, ParameterList params,
                      ConstVectorPtr const &x0);

    int run();
    int total_newton_steps() { return total_newton_steps_; }

private:
    void writeData();
};

template<typename Model>
template<typename ParameterList>
AdaptiveTransient<Model>::
AdaptiveTransient(Model model, ParameterList params)
    :
    Transient<ConstVectorPtr>(),
    model_(model),
    newton_(Teuchos::rcp(new Newton<Model>(model, params))),
    adaptive_(params->get("adaptive time steps", false)),
    min_wanted_newton_steps_(params->get("minimum desired Newton iterations", 3)), 
    max_wanted_newton_steps_(params->get("maximum desired Newton iterations", 3)),
    min_time_step_(params->get("minimum time step", 1.0e-8)),
    max_time_step_(params->get("maximum time step", 1.0)),
    time_step_increase_(params->get("time step increase", 2.0)),
    time_step_decrease_(params->get("time step decrease", 2.0)),
    nsteps_(params->get("number of time steps", 10)),
    output_(params->get("HDF5 output frequency", 1)),
    total_newton_steps_(0),
    init_wd_(true)
{
    Transient<ConstVectorPtr>::set_parameters(*params);
}

template<typename Model>
template<typename ParameterList>
AdaptiveTransient<Model>::
AdaptiveTransient(Model model, ParameterList params,
                  ConstVectorPtr const &x0)
    :
    AdaptiveTransient(model, params)
{
    Transient<ConstVectorPtr>::x0_ = new ConstVectorPtr(x0);
}

//==================================================================
template<typename Model>
int AdaptiveTransient<Model>::run()
{
    ConstVectorPtr x = model_->getState('C');
    INFO("||x|| = " << Utils::norm(x));
    if (Transient<ConstVectorPtr>::x0_)
        x = *Transient<ConstVectorPtr>::x0_;

    Transient<ConstVectorPtr>::time_steps_ = 0;
    time_ = 0;

    bool test_step = ( nsteps_ < 0 ) ? true : Transient<ConstVectorPtr>::time_steps_ < nsteps_;
    while ( ( time_ < Transient<ConstVectorPtr>::tmax_ )
            && ( test_step ) )
    {
        INFO("----------------------------------------------------------");
        INFO("Timestepping:    t = "
             << time_ * Transient<ConstVectorPtr>::in_years_
             << " y");

        model_->setState(x);
        model_->initStep(Transient<ConstVectorPtr>::dt_);

        ConstVectorPtr y = newton_->run(x);

        if (!newton_->converged())
        {
            WARNING("Newton did not converge! ||F|| = "
                    << newton_->normF() << "\nRestoring model",
                    __FILE__, __LINE__);

            Utils::save(newton_->Fx(), "failed_rhs"); // Print failed residual

            if (Transient<ConstVectorPtr>::dt_ == min_time_step_ || !adaptive_)
            {
                ERROR("\n ================================ \n  " <<
                      " minimum timestep reached, exiting... " <<
                      "\n ================================ \n",
                      __FILE__, __LINE__);
                return 1;
            }

            INFO("    adjusting time step.. old dt = " << Transient<ConstVectorPtr>::dt_);
            Transient<ConstVectorPtr>::dt_ = std::max(Transient<ConstVectorPtr>::dt_ / time_step_decrease_, min_time_step_);
            INFO("    adjusting time step.. new dt = " << Transient<ConstVectorPtr>::dt_);
            continue;
        }

        Transient<ConstVectorPtr>::time_steps_++;
        time_ += Transient<ConstVectorPtr>::dt_;
        x = y;

        INFO("  Newton converged, next time step -----------------");
        INFO("                step = " << Transient<ConstVectorPtr>::time_steps_);
        INFO("                time = " << time_);
        INFO("           ||F||2    = " << newton_->normF());
        INFO("           ||dx||inf = " << newton_->normdx());
        INFO("\n");

        model_->postProcess();

        // additional save
        if ((output_ > 0) && ((Transient<ConstVectorPtr>::time_steps_ % output_) == 0))
        {
            std::stringstream outFile;
            outFile << "transient_" << std::setprecision(8)
                    << time_ << ".h5";

            model_->saveStateToFile(outFile.str());
        }

        writeData();

        // Timestep adjustments
        if (adaptive_ && newton_->steps() < min_wanted_newton_steps_)
            Transient<ConstVectorPtr>::dt_ = std::min(Transient<ConstVectorPtr>::dt_ * time_step_increase_, max_time_step_);
        else if (adaptive_ && newton_->steps() > max_wanted_newton_steps_)
            Transient<ConstVectorPtr>::dt_ = std::max(Transient<ConstVectorPtr>::dt_ / time_step_decrease_, min_time_step_);

        test_step = ( nsteps_ < 0 ) ? true : Transient<ConstVectorPtr>::time_steps_ < nsteps_;

        total_newton_steps_ += newton_->steps();
    }
    return 0;
}

//==================================================================
template<typename Model>
void AdaptiveTransient<Model>::writeData()
{
    // Data string
    std::ostringstream tdatastring;

    if (init_wd_) // write description of entries
    {
        tdatastring << std::setw(_FIELDWIDTH_)
                    << "# time_(y)"
                    << std::setw(_FIELDWIDTH_/2.)
                    << "step"
                    << std::setw(_FIELDWIDTH_)
                    << "Transient<ConstVectorPtr>::dt_(y)"
                    << std::setw(_FIELDWIDTH_)
                    << "|x|"
                    << std::setw(_FIELDWIDTH_/2.)
                    << "NR"
                    << model_->writeData(true);

        WRITETDATA(tdatastring.str());
        init_wd_ = false;
    }
    tdatastring.str("");
    tdatastring.clear();
    tdatastring.precision(_PRECISION_);

    tdatastring << std::scientific
                << std::setw(_FIELDWIDTH_)
                << time_ * Transient<ConstVectorPtr>::in_years_
                << std::setw(_FIELDWIDTH_/2)
                << Transient<ConstVectorPtr>::time_steps_
                << std::setw(_FIELDWIDTH_)
                << Transient<ConstVectorPtr>::dt_ * Transient<ConstVectorPtr>::in_years_
                << std::setw(_FIELDWIDTH_)
                << Utils::norm(model_->getState('V'))
                << std::setw(_FIELDWIDTH_/2) <<  newton_->steps()
                << model_->writeData();

    WRITETDATA(tdatastring.str());
}

#endif
