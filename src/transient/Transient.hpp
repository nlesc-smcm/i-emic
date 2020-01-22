#ifndef TRANSIENT_HPP
#define TRANSIENT_HPP

#include "TransientDecl.hpp"

#include "GlobalDefinitions.H"

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>

template<class T>
struct AMSExperiment {
    T x0;

    std::vector<T> xlist;
    std::vector<double> dlist;
    std::vector<double> tlist;

    double max_distance;
    double time;
    double initial_time;
    double return_time;

    bool initialized;
    bool converged;

    AMSExperiment()
        :
        x0(),
        xlist(),
        dlist(),
        tlist(),
        max_distance(0.0),
        time(0.0),
        initial_time(0.0),
        return_time(0.0),
        initialized(false),
        converged(false)
        {}

    static bool sort(AMSExperiment<T> const *e1, AMSExperiment<T> const *e2)
        {
            return e1->max_distance > e2->max_distance;
        };
};

template<class T>
struct GPAExperiment {
    T x;
    double weight;
    double probability;
    double distance;
    bool converged;
};

std::string mem2string(long long mem);

template<class T>
Transient<T>::Transient()
    :
    method_("Transient"),
    x0_(nullptr),
    mfpt_(-1),
    probability_(-1),
    engine_initialized_(false)
{}

template<class T>
Transient<T>::Transient(std::function<T(T const &, double)> time_step)
    :
    time_step_(time_step),
    method_("Transient"),
    x0_(nullptr),
    mfpt_(-1),
    probability_(-1),
    engine_initialized_(false)
{}

template<class T>
Transient<T>::Transient(std::function<T(T const &, double)> time_step,
                        T const &x0)
    :
    time_step_(time_step),
    method_("Transient"),
    x0_(new T(x0)),
    mfpt_(-1),
    probability_(-1),
    engine_initialized_(false)
{}

template<class T>
Transient<T>::Transient(std::function<T(T const &, double)> time_step,
                        std::function<double(T const &)> dist_fun,
                        int vector_length)
    :
    time_step_(time_step),
    dist_fun_(dist_fun),
    method_("TAMS"),
    x0_(nullptr),
    vector_length_(vector_length),
    mfpt_(-1),
    probability_(-1),
    engine_initialized_(false)
{}

template<class T>
Transient<T>::Transient(std::function<T(T const &, double)> time_step,
                        std::function<double(T const &)> dist_fun,
                        T const &x0,
                        int vector_length)
    :
    time_step_(time_step),
    dist_fun_(dist_fun),
    method_("TAMS"),
    x0_(new T(x0)),
    vector_length_(vector_length),
    mfpt_(-1),
    probability_(-1),
    engine_initialized_(false)
{}

template<class T>
Transient<T>::~Transient()
{
    if (engine_initialized_)
        delete engine_;

    if (x0_)
        delete x0_;
}

template<class T>
template<class ParameterList>
void Transient<T>::set_parameters(ParameterList &params)
{
    method_ = params.get("method", method_);
    dt_ = params.get("time step", 0.01);
    tmax_ = params.get("maximum time", 1000.0);

    in_days_ = params.get("timescale in days", 737.2685);
    in_years_ = params.get("timescale in years", in_days_ / 365.);
    dt_ = params.get("time step (in y)", dt_ * in_years_) / in_years_;
    tmax_ = params.get("maximum time (in y)", tmax_ * in_years_) / in_years_;

    // GPA parameters
    tstep_ = params.get("GPA time step", 1.0);
    beta_ = params.get("beta", 1.0);

    // (T)AMS + GPA parameters
    bdist_ = params.get("B distance", 0.05);
    dist_tol_ = params.get("distance tolerance", 0.0005);
    num_exp_ = params.get("number of experiments", 1000);

    // AMS parameters
    adist_ = params.get("A distance", 0.05);
    cdist_ = params.get("C distance", 2 * adist_);
    num_init_exp_ = params.get("number of initial experiments", num_exp_);
    if (num_init_exp_ < num_exp_)
        num_init_exp_ = num_exp_;

    // (T)AMS parameters
    maxit_ = params.get("maximum iterations", num_exp_ * 10);

    // Writing parameters
    read_ = params.get("read file", "");
    write_ = params.get("write file", "");
    write_final_ = params.get("write final state", true);
    write_steps_ = params.get("write steps", -1);
    write_time_steps_ = params.get("write time steps", -1);
}

template<class T>
T Transient<T>::transient(T x, double dt, double tmax) const
{
    for (double t = dt; t <= tmax; t += dt)
        x = std::move(time_step_helper(x, dt));
    return x;
}

template<class T>
double Transient<T>::transient_max_distance(
    T x, double dt, double tmax, double max_distance) const
{
    const double lim = max_distance - bdist_;
    for (double t = dt; t <= tmax; t += dt)
    {
        x = std::move(time_step_helper(x, dt));
        if (dist_fun_(x) > lim)
            return t;
    }
    return -1;
}

template<class T>
void Transient<T>::transient_start(
    T const &x0, double dt, double tmax,
    AMSExperiment<T> &experiment) const
{
    T x(x0);

    experiment.initial_time = 0.0;

    for (double t = dt; t <= tmax; t += dt)
    {
        x = std::move(time_step_helper(x, dt));

        double dist = dist_fun_(x);
        if (dist > cdist_)
        {
            experiment.xlist.push_back(x);
            experiment.dlist.push_back(dist);
            experiment.tlist.push_back(0);
            experiment.max_distance = dist;
            experiment.initial_time = t;
            experiment.initialized = true;
            break;
        }
    }
}

template<class T>
void Transient<T>::transient_ams(
    double dt, double tmax,
    AMSExperiment<T> &experiment) const
{
    T x(experiment.xlist.back());
    double t = experiment.tlist.back() + dt;
    double tend = t + tmax;

    double max_distance = experiment.max_distance;

    for (; t <= tend; t += dt)
    {
        x = std::move(time_step_helper(x, dt));

        double dist = dist_fun_(x);

        if (dist < adist_)
        {
            if (experiment.return_time < dt / 2.0)
                experiment.return_time = t;
            break;
        }
        else if (dist > 1 - bdist_)
        {
            experiment.converged = true;
            experiment.xlist.push_back(x);
            experiment.tlist.push_back(t);
            experiment.dlist.push_back(1.0);
            max_distance = 1.0;
            break;
        }
        if (dist > max_distance + dist_tol_)
        {
            experiment.xlist.push_back(x);
            experiment.tlist.push_back(t);
            experiment.dlist.push_back(dist);
            max_distance = dist;
        }
    }

    experiment.max_distance = max_distance;
    experiment.time = t;
}

template<class T>
void Transient<T>::transient_tams(
    double dt, double tmax,
    AMSExperiment<T> &experiment) const
{
    T x(experiment.xlist.back());
    double t = experiment.tlist.back() + dt;

    double max_distance = experiment.max_distance;

    for (; t <= tmax; t += dt)
    {
        x = std::move(time_step_helper(x, dt));

        double dist = dist_fun_(x);
        if (dist > 1 - bdist_)
        {
            experiment.converged = true;
            experiment.xlist.push_back(x);
            experiment.tlist.push_back(t);
            experiment.dlist.push_back(1.0);
            max_distance = 1.0;
            break;
        }
        if (dist > max_distance + dist_tol_)
        {
            experiment.xlist.push_back(x);
            experiment.tlist.push_back(t);
            experiment.dlist.push_back(dist);
            max_distance = dist;
        }
    }

    experiment.time = experiment.tlist.back();
    experiment.max_distance = max_distance;
}

template<class T>
void Transient<T>::transient_gpa(
    double dt, double tmax,
    GPAExperiment<T> &experiment) const
{
    T x(experiment.x);
    double dist = -1;

    for (double t = dt; t <= tmax; t += dt)
    {
        x = std::move(time_step_helper(x, dt));

        dist = dist_fun_(x);
        if (dist > 1 - bdist_)
            experiment.converged = true;
    }

    experiment.distance = dist;
    experiment.x = x;
}

template<class T>
void Transient<T>::naive(T const &x0) const
{
    std::vector<GPAExperiment<T>> experiments(num_exp_);

    int converged = 0;
    for (int i = 0; i < num_exp_; i++)
    {
        experiments[i].x = x0;
        experiments[i].converged = false;
        transient_gpa(dt_, tmax_, experiments[i]);

        if (experiments[i].converged)
            converged++;
    }

    probability_ = (double)converged / (double)num_exp_;

    INFO("Transition probability T=" << tmax_ << ": " << probability_);
}


template<class T>
double Transient<T>::ams_elimination(
    std::string const &method,
    std::vector<AMSExperiment<T>> &experiments,
    double dt, double tmax) const
{
    int converged = 0;

    std::vector<AMSExperiment<T> *> reactive_experiments;
    std::vector<AMSExperiment<T> *> unconverged_experiments;
    std::vector<AMSExperiment<T> *> unused_experiments;
    std::vector<AMSExperiment<T> *> minimal_experiments;

    for (int i = 0; i < num_exp_; i++)
        reactive_experiments.push_back(&experiments[i]);

    for (auto exp: reactive_experiments)
    {
        if (!exp->converged)
            unconverged_experiments.push_back(exp);
        else
            converged++;
        unused_experiments.push_back(exp);
    }

    std::sort(unconverged_experiments.begin(),
              unconverged_experiments.end(), AMSExperiment<T>::sort);

    for (int i = its_; i < maxit_; i++)
    {
        minimal_experiments.clear();
        if (unconverged_experiments.size() > 0 && unused_experiments.size() > 0)
        {
            AMSExperiment<T> *minimal_exp = unconverged_experiments.back();
            AMSExperiment<T> *exp = unconverged_experiments.back();
            while (unconverged_experiments.size() > 0 &&
                   exp->max_distance == minimal_exp->max_distance)
            {
                minimal_experiments.push_back(exp);
                unconverged_experiments.pop_back();
                unused_experiments.erase(
                    std::find(unused_experiments.begin(),
                              unused_experiments.end(), exp));
                if (unconverged_experiments.size() > 0)
                    exp = unconverged_experiments.back();
            }
        }

        if (minimal_experiments.size() == 0 || unused_experiments.size() == 0)
            continue;

        ell_.push_back(minimal_experiments.size());
        if (ell_.back() == 1)
        {
            INFO("Eliminating 1 trajectory.");
        }
        else
        {
            INFO("Eliminating " << ell_.back() << " trajectories.");
        }

        its_++;

        for (auto &exp: minimal_experiments)
        {
            double max_distance = exp->max_distance;
            int rnd_idx = randint(0, unused_experiments.size()-1);
            while (unused_experiments[rnd_idx]->max_distance <= exp->max_distance)
                rnd_idx = randint(0, unused_experiments.size()-1);

            AMSExperiment<T> *rnd_exp = unused_experiments[rnd_idx];
            if (rnd_exp->dlist.size() == 0)
            {
                ERROR("Experiment " << rnd_idx << " has size 0.", __FILE__, __LINE__);
            }

            int same_distance_idx = -1;
            while (++same_distance_idx < (int)rnd_exp->dlist.size() &&
                   rnd_exp->dlist[same_distance_idx] < exp->max_distance);

            if (same_distance_idx == (int)rnd_exp->dlist.size())
            {
                ERROR("Distance larger than " << exp->max_distance
                      << " not found in experiment with max distance "
                      << rnd_exp->max_distance << ".", __FILE__, __LINE__);
            }

            exp->xlist = std::vector<T>(
                rnd_exp->xlist.begin(), rnd_exp->xlist.begin() + same_distance_idx + 1);
            exp->dlist = std::vector<double>(
                rnd_exp->dlist.begin(), rnd_exp->dlist.begin() + same_distance_idx + 1);
            exp->tlist = std::vector<double>(
                rnd_exp->tlist.begin(), rnd_exp->tlist.begin() + same_distance_idx + 1);

            if (method == "AMS")
                transient_ams(dt, tmax, *exp);
            else if (method == "TAMS")
                transient_tams(dt, tmax, *exp);
            else
            {
                ERROR("Method " << method << " does not exist.", __FILE__, __LINE__);
            }

            if (exp->converged)
                converged++;
            else
                unconverged_experiments.push_back(exp);

            INFO(method << ": " << its_ << " / " << maxit_ << ", "
                 << converged << " / " << num_exp_
                 << " converged with max distance "
                 << max_distance << " -> "
                 << exp->max_distance << " and t="
                 << exp->initial_time + exp->time
                 << " for experiment "
                 << std::find(reactive_experiments.begin(),
                              reactive_experiments.end(),
                              exp) - reactive_experiments.begin());
        }

        for (auto &exp: minimal_experiments)
            unused_experiments.push_back(exp);

        std::sort(unconverged_experiments.begin(),
                  unconverged_experiments.end(), AMSExperiment<T>::sort);

        // Cleanup
        double min_max_distance = 1.0;
        for (auto exp: reactive_experiments)
            min_max_distance = std::min(min_max_distance, exp->max_distance);

        if (its_ % 10 == 0)
        {
            INFO("Starting cleanup");
            for (auto exp: unused_experiments)
            {
                int min_max_idx = -1;
                while (++min_max_idx < (int)exp->dlist.size() &&
                       exp->dlist[min_max_idx] < min_max_distance);

                if (min_max_idx > 0)
                {
                    exp->xlist = std::vector<T>(
                        exp->xlist.begin() + min_max_idx,
                        exp->xlist.end());
                    exp->dlist = std::vector<double>(
                        exp->dlist.begin() + min_max_idx,
                        exp->dlist.end());
                    exp->tlist = std::vector<double>(
                        exp->tlist.begin() + min_max_idx,
                        exp->tlist.end());
                }
            }
            INFO("Finished cleanup");
        }

        write_helper(experiments, its_);
    }
    INFO("");

    if (write_final_)
        write(write_, experiments);

    double alpha = (double)converged / (double)num_exp_;
    for (int l: ell_)
        alpha *= 1.0 - (double)l / (double)num_exp_;

    return alpha;
}

template<class T>
void Transient<T>::ams(T const &x0) const
{
    INFO("Using AMS with an estimated memory usage of: "
         << mem2string((long long)(1.0 / dist_tol_ * num_exp_ *
                                   vector_length_ * sizeof(double))));

    std::vector<AMSExperiment<T>> experiments(num_init_exp_);
    for (int i = 0; i < num_init_exp_; i++)
        experiments[i].x0 = x0;

    its_ = 0;
    time_steps_ = 0;
    ell_.clear();

    if (read_ != "")
        read(read_, experiments);

    int converged = 0;
    double tmax = 100 * tmax_;
    time_steps_previous_write_ = 0;

    for (int i = 0; i < num_init_exp_; i++)
    {
        if (experiments[i].initialized)
            continue;

        transient_start(x0, dt_, tmax, experiments[i]);

        if (experiments[i].xlist.size() == 0)
        {
            ERROR("Initialization failed", __FILE__, __LINE__);
        }

        transient_ams(dt_, tmax, experiments[i]);

        // Erase data that we do not need for later experiments
        if (i >= num_exp_)
        {
            experiments[i].xlist = std::vector<T>();
            experiments[i].dlist = std::vector<double>();
            experiments[i].tlist = std::vector<double>();
        }

        if (experiments[i].converged)
            converged++;

        INFO("Initialization: " << i+1 << " / " << num_init_exp_ << ", "
             << converged << " / " << num_init_exp_
             << " converged with t="
             << experiments[i].initial_time + experiments[i].time);

        write_helper(experiments, i+1);
    }
    INFO("");

    double alpha = ams_elimination("AMS", experiments, dt_, tmax);

    double total_tr = 0;
    double total_t1 = 0;
    double total_t2 = 0;
    int num_t1 = num_init_exp_;
    int num_t2 = 0;
    converged = 0;
    for (int i = 0; i < num_exp_; i++)
    {
        auto &exp = experiments[i];
        total_tr += exp.time;
        converged += exp.converged;
    }

    for (auto &exp: experiments)
    {
        total_t1 += exp.initial_time;
        total_t2 += exp.return_time;
        if (exp.return_time > dt_ / 2.0)
            num_t2++;
    }

    double meann = 1.0 / alpha - 1.0;
    mfpt_ = meann * (total_t1 / (double)num_t1 + total_t2 / (double)num_t2) +
        total_t1 / (double)num_t1 + total_tr / (double)converged;
    INFO("Alpha: " << alpha);
    INFO("Mean first passage time: " << mfpt_);

    probability_ = 1.0 -  exp(-1.0 / mfpt_ * tmax_);
    INFO("Transition probability T=" << tmax_ << ": " << probability_);
}

template<class T>
void Transient<T>::tams(T const &x0) const
{
    INFO("Using TAMS with an estimated memory usage of: "
         << mem2string((long long)(tmax_ / dt_ * num_exp_ *
                                   vector_length_ * sizeof(double))));

    std::vector<AMSExperiment<T>> experiments(num_exp_);
    for (int i = 0; i < num_exp_; i++)
        experiments[i].x0 = x0;

    its_ = 0;
    time_steps_ = 0;
    ell_.clear();

    if (read_ != "")
        read(read_, experiments);

    int converged = 0;
    time_steps_previous_write_ = 0;

    for (int i = 0; i < num_exp_; i++)
    {
        if (experiments[i].initialized)
            continue;

        experiments[i].xlist.push_back(x0);
        experiments[i].dlist.push_back(0);
        experiments[i].tlist.push_back(0);

        transient_tams(dt_, tmax_, experiments[i]);

        experiments[i].initialized = true;

        if (experiments[i].converged)
            converged++;

        INFO("Initialization: " << i+1 << " / " << num_exp_ << ", "
             << converged << " / " << num_exp_
             << " converged with t="
             << experiments[i].time);

        write_helper(experiments, i+1);
    }
    INFO("");

    probability_ = ams_elimination("TAMS", experiments, dt_, tmax_);

    INFO("Transition probability T=" << tmax_ << ": " << probability_);
}

template<class T>
void Transient<T>::gpa(T const &x0) const
{
    INFO("Using GPA with an estimated memory usage of: "
         << mem2string((long long)(num_exp_ * vector_length_ * sizeof(double))));

    std::vector<GPAExperiment<T>> experiments(num_exp_);

    time_steps_ = 0;

    auto W = [this](double x){return exp(beta_ * x);};

    for (int i = 0; i < num_exp_; i++)
    {
        experiments[i].x = x0;
        experiments[i].weight = 1.0;
        experiments[i].probability = 1.0;
        experiments[i].distance = 0.0;
        experiments[i].converged = false;
    }

    for (double t = tstep_; t <= tmax_; t += tstep_)
    {
        // Compute the mean weight
        double sum = 0.0;
        for (int i = 0; i < num_exp_; i++)
            sum += experiments[i].weight;
        double eta = 1.0 / (double)num_exp_ * sum;

        std::vector<GPAExperiment<T>> old_experiments = experiments;

        // Sample based on the weights
        for (int i = 0; i < num_exp_; i++)
        {
            double val = randreal(0.0, sum);
            double  cumsum = 0.0;
            for (int j = 0; j < num_exp_; j++)
            {
                cumsum += old_experiments[j].weight;
                if (cumsum >= val)
                {
                    experiments[i] = old_experiments[j];
                    break;
                }
                if (j == num_exp_-1)
                {
                    ERROR("Particle not found with sum " << sum
                          << " and random value " << val, __FILE__, __LINE__);
                }
            }
        }

        // Step until the next tstep and recompute weights
        int converged = 0;
        for (int i = 0; i < num_exp_; i++)
        {
            transient_gpa(dt_, tstep_, experiments[i]);
            experiments[i].weight = W(experiments[i].distance);
            experiments[i].probability *= eta / experiments[i].weight;

            if (experiments[i].converged)
                converged++;
        }

        INFO("GPA: " << converged << " / " << num_exp_
             << " converged with t=" << t << " and eta=" << eta);
    }
    INFO("");

    probability_ = 0.0;
    for (int i = 0; i < num_exp_; i++)
        if (experiments[i].converged)
            probability_ += experiments[i].probability;
    probability_ /= (double)num_exp_;

    INFO("Transition probability T=" << tmax_ << ": " << probability_);
}

template<class T>
void Transient<T>::set_random_engine(unsigned int seed)
{
    if (engine_initialized_)
        delete engine_;

    std::seed_seq seeder{seed};
    engine_ = new std::mt19937_64(seeder);

    randint_ = [this](int a, int b) {
        std::uniform_int_distribution<int> int_distribution(a, b);
        int val = int_distribution(*engine_);
        return val;
    };
    randreal_ = [this](int a, int b) {
        std::uniform_real_distribution<double> real_distribution(a, b);
        double val = real_distribution(*engine_);
        return val;
    };
    engine_initialized_ = true;
}

template<class T>
int Transient<T>::randint(int a, int b) const
{
    if (engine_initialized_)
        return randint_(a, b);

    static bool first = true;
    if (first)
    {
        first = false;
        WARNING("Random engine not initialized.", __FILE__, __LINE__);
    }

    static thread_local std::random_device rd;
    static thread_local std::seed_seq seeder{rd()};
    static thread_local std::mt19937_64 engine(seeder);
    std::uniform_int_distribution<int> int_distribution(a, b);
    return int_distribution(engine);
}

template<class T>
int Transient<T>::randreal(double a, double b) const
{
    if (engine_initialized_)
        return randreal_(a, b);

    static bool first = true;
    if (first)
    {
        first = false;
        WARNING("Random engine not initialized.", __FILE__, __LINE__);
    }

    static thread_local std::random_device rd;
    static thread_local std::seed_seq seeder{rd()};
    static thread_local std::mt19937_64 engine(seeder);
    std::uniform_real_distribution<double> real_distribution(a, b);
    return real_distribution(engine);
}

template<class T>
T Transient<T>::time_step_helper(T const &x, double dt) const
{
    time_steps_++;
    return std::move(time_step_(x, dt));
}

template<class T>
void Transient<T>::write_helper(std::vector<AMSExperiment<T> > const &experiments,
                                int its) const
{
    if (write_ == "")
        return;

    if (write_steps_ > 0 && its % write_steps_ == 0)
    {
        time_steps_previous_write_ = time_steps_;
        write(write_, experiments);
        return;
    }

    if (write_time_steps_ > 0 && time_steps_ - time_steps_previous_write_ >= write_time_steps_)
    {
        time_steps_previous_write_ = time_steps_;
        write(write_, experiments);
        return;
    }
}

template<typename T>
int Transient<T>::run()
{
    return run(*x0_);
}

template<typename T>
int Transient<T>::run(T const &x0)
{
    if (method_ == "AMS")
        ams(x0);
    else if (method_ == "TAMS")
        tams(x0);
    else if (method_ == "GPA")
        gpa(x0);
    else if (method_ == "Naive")
        naive(x0);
    else if (method_ == "Transient")
        transient(x0, dt_, tmax_);
    else
    {
        ERROR("Method " << method_ << " does not exist.", __FILE__, __LINE__);
        return -1;
    }
    return 0;
}

template<typename T>
void Transient<T>::read(
    std::string const &name,
    std::vector<AMSExperiment<T> > &experiments) const
{
    WARNING("Reading transient data not implemented.", __FILE__, __LINE__);
}

template<typename T>
void Transient<T>::write(
    std::string const &name,
    std::vector<AMSExperiment<T> > const &experiments) const
{
    WARNING("Writing transient data not implemented.", __FILE__, __LINE__);
}

template<class T>
double Transient<T>::get_probability()
{
    return probability_;
}

template<class T>
double Transient<T>::get_mfpt()
{
    return mfpt_;
}

#include "Trilinos_version.h"

#if TRILINOS_MAJOR_MINOR_VERSION > 121300

namespace Teuchos { template<class T> class RCP; }
class Epetra_Vector;

template<>
void Transient<Teuchos::RCP<const Epetra_Vector> >::read(
    std::string const &name,
    std::vector<AMSExperiment<Teuchos::RCP<const Epetra_Vector> > > &experiments) const;

template<>
void Transient<Teuchos::RCP<const Epetra_Vector> >::write(
    std::string const &name,
    std::vector<AMSExperiment<Teuchos::RCP<const Epetra_Vector> > > const &experiments) const;

#endif //TRILINOS_MAJOR_MINOR_VERSION

#endif
