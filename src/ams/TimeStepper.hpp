#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include "TimeStepperDecl.hpp"

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

std::string mem2string(long long mem)
{
    double value = mem;
    std::string unit = "B";
    if (std::abs(value) > 1.0e3) {value *= 1.0e-3; unit = "kB";}
    if (std::abs(value) > 1.0e3) {value *= 1.0e-3; unit = "MB";}
    if (std::abs(value) > 1.0e3) {value *= 1.0e-3; unit = "GB";}
    if (std::abs(value) > 1.0e3) {value *= 1.0e-3; unit = "TB";}

    std::ostringstream ss;
    ss << std::fixed;
    ss.precision(2);
    ss << value << " " << unit;
    return ss.str();
}

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &, double)> time_step)
    :
    time_step_(time_step),
    engine_initialized_(false)
{}

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &, double)> time_step,
                            std::function<double(T const &)> dist_fun,
                            int vector_length)
    :
    time_step_(time_step),
    dist_fun_(dist_fun),
    vector_length_(vector_length),
    engine_initialized_(false)
{}

template<class T>
TimeStepper<T>::~TimeStepper()
{
    if (engine_initialized_)
        delete engine_;
}

template<class T>
template<class ParameterList>
void TimeStepper<T>::set_parameters(ParameterList &params)
{
    dt_ = params.get("time step", 0.01);
    tstep_ = params.get("GPA time step", 1.0);
    tmax_ = params.get("maximum time", 1000.0);
    beta_ = params.get("beta", 1.0);
    adist_ = params.get("A distance", 0.05);
    bdist_ = params.get("B distance", adist_);
    cdist_ = params.get("C distance", 2 * adist_);
    dist_tol_ = params.get("distance tolerance", 0.0005);
    num_exp_ = params.get("number of experiments", 1000);
    num_init_exp_ = params.get("number of initial experiments", num_exp_);
    maxit_ = params.get("maximum iterations", num_exp_ * 10);
    read_ = params.get("read file", "");
    write_ = params.get("write file", "");
    write_final_ = params.get("write final state", true);
    write_steps_ = params.get("write steps", -1);
    write_time_steps_ = params.get("write time steps", -1);

    if (num_init_exp_ < num_exp_)
        num_init_exp_ = num_exp_;
}

template<class T>
T TimeStepper<T>::transient(T x, double dt, double tmax) const
{
    for (double t = dt; t <= tmax; t += dt)
        x = std::move(time_step_helper(x, dt));
    return x;
}

template<class T>
double TimeStepper<T>::transient_max_distance(
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
void TimeStepper<T>::transient_start(
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
void TimeStepper<T>::transient_ams(
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
void TimeStepper<T>::transient_tams(
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
void TimeStepper<T>::transient_gpa(
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
void TimeStepper<T>::naive(T const &x0) const
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

    double tp = (double)converged / (double)num_exp_;

    std::cout << "Transition probability T=" << tp << ": " << tp << std::endl;
}

template<class T>
void TimeStepper<T>::ams(T const &x0) const
{
    std::cout << "Using AMS with an estimated memory usage of: "
              << mem2string((long long)(1.0 / dist_tol_ * num_exp_ *
                                        vector_length_ * sizeof(double)))
              << std::endl;

    std::vector<AMSExperiment<T>> experiments(num_init_exp_);
    for (int i = 0; i < num_init_exp_; i++)
        experiments[i].x0 = x0;

    its_ = 0;
    time_steps_ = 0;

    if (read_ != "")
        read(read_, experiments);

    int converged = 0;
    double tmax = 100 * tmax_;
    int time_steps_previous_write = 0;

    for (int i = 0; i < num_init_exp_; i++)
    {
        if (experiments[i].initialized)
            continue;

        transient_start(x0, dt_, tmax, experiments[i]);

        if (experiments[i].xlist.size() == 0)
        {
            std::cerr << "Initialization failed" << std::endl;
            exit(-1);
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

        std::cout << "Initialization: " << i+1 << " / " << num_init_exp_ << ", "
                  << converged << " / " << num_init_exp_
                  << " converged with t="
                  << experiments[i].initial_time + experiments[i].time << std::endl;

        write_helper(experiments, i+1, time_steps_previous_write);
    }
    std::cout << std::endl;

    converged = 0;

    std::vector<AMSExperiment<T> *> unconverged_experiments;
    std::vector<AMSExperiment<T> *> unused_experiments;
    std::vector<AMSExperiment<T> *> reactive_experiments;

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

    for (int i = 0; i < maxit_; i++)
    {
        AMSExperiment<T> *exp = NULL;
        int num_unused_exp = 0;


        if (unconverged_experiments.size() > 0 && unused_experiments.size() > 0)
        {
            exp = unconverged_experiments.back();
            unconverged_experiments.pop_back();
            unused_experiments.erase(
                std::find(unused_experiments.begin(),
                          unused_experiments.end(), exp));
        }
        num_unused_exp = unused_experiments.size();

        if (exp == NULL || num_unused_exp == 0)
            continue;

        double max_distance = exp->max_distance;


        int rnd_idx = randint(0, unused_experiments.size()-1);
        while (unused_experiments[rnd_idx]->max_distance <= exp->max_distance)
            rnd_idx = randint(0, unused_experiments.size()-1);

        AMSExperiment<T> *rnd_exp = unused_experiments[rnd_idx];
        if (rnd_exp->dlist.size() == 0)
        {
            std::cerr << "Experiment " << rnd_idx << " has size 0." << std::endl;
            exit(-1);
        }

        int same_distance_idx = -1;
        while (++same_distance_idx < (int)rnd_exp->dlist.size() &&
               rnd_exp->dlist[same_distance_idx] < exp->max_distance);

        if (same_distance_idx == (int)rnd_exp->dlist.size())
        {
            std::cerr << "Distance larger than " << exp->max_distance
                      << " not found in experiment with max distance "
                      << rnd_exp->max_distance << "." << std::endl;
            exit(-1);
        }

        exp->xlist = std::vector<T>(
            rnd_exp->xlist.begin(), rnd_exp->xlist.begin() + same_distance_idx + 1);
        exp->dlist = std::vector<double>(
            rnd_exp->dlist.begin(), rnd_exp->dlist.begin() + same_distance_idx + 1);
        exp->tlist = std::vector<double>(
            rnd_exp->tlist.begin(), rnd_exp->tlist.begin() + same_distance_idx + 1);

        transient_ams(dt_, tmax, *exp);

        its_++;
        if (exp->converged)
            converged++;
        else
            unconverged_experiments.push_back(exp);

        unused_experiments.push_back(exp);

        std::sort(unconverged_experiments.begin(),
                  unconverged_experiments.end(), AMSExperiment<T>::sort);

        std::cout << "AMS: " << its_ << " / " << maxit_ << ", "
                  << converged << " / " << num_exp_
                  << " converged with max distance "
                  << max_distance << " -> "
                  << exp->max_distance << " and t="
                  << exp->initial_time + exp->time << std::endl;

        // Cleanup
        double min_max_distance = 1.0;
        for (auto exp: reactive_experiments)
            min_max_distance = std::min(min_max_distance, exp->max_distance);

        if (its_ % 10 == 0)
        {
            std::cout << "Starting cleanup" << std::endl;
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
            std::cout << "Finished cleanup" << std::endl;
        }

        write_helper(experiments, its_, time_steps_previous_write);
    }
    std::cout << std::endl;

    if (write_final_)
        write(write_, experiments);

    double total_tr = 0;
    double total_t1 = 0;
    double total_t2 = 0;
    int num_t1 = num_init_exp_;
    int num_t2 = 0;
    for (auto exp: reactive_experiments)
        total_tr += exp->time;

    for (auto &exp: experiments)
    {
        total_t1 += exp.initial_time;
        total_t2 += exp.return_time;
        if (exp.return_time > dt_ / 2.0)
            num_t2++;
    }

    double alpha = (double)converged / (double)num_exp_ *
        pow((1.0 - 1.0 / (double)num_exp_), its_);

    double meann = 1.0 / alpha - 1.0;
    double fpt = meann * (total_t1 / (double)num_t1 + total_t2 / (double)num_t2) +
        total_t1 / (double)num_t1 + total_tr / (double)converged;
    std::cout << "Alpha: " << alpha << std::endl;
    std::cout << "Mean first passage time: " << fpt << std::endl;

    double tp = 1.0 -  exp(-1.0 / fpt * tmax_);
    std::cout << "Transition probability T=" << tmax_ << ": " << tp << std::endl;
}

template<class T>
void TimeStepper<T>::tams(T const &x0) const
{
    std::cout << "Using TAMS with an estimated memory usage of: "
              << mem2string((long long)(tmax_ / dt_ * num_exp_ *
                                        vector_length_ * sizeof(double)))
              << std::endl;

    std::vector<AMSExperiment<T>> experiments(num_exp_);
    for (int i = 0; i < num_exp_; i++)
        experiments[i].x0 = x0;

    its_ = 0;
    time_steps_ = 0;

    if (read_ != "")
        read(read_, experiments);

    int converged = 0;
    int time_steps_previous_write = 0;

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

        std::cout << "Initialization: " << i+1 << " / " << num_exp_ << ", "
                  << converged << " / " << num_exp_
                  << " converged with t="
                  << experiments[i].time << std::endl;

        write_helper(experiments, i+1, time_steps_previous_write);
    }
    std::cout << std::endl;
    converged = 0;

    std::vector<AMSExperiment<T> *> unconverged_experiments;
    std::vector<AMSExperiment<T> *> unused_experiments;
    std::vector<AMSExperiment<T> *> reactive_experiments;

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

    for (int i = 0; i < maxit_; i++)
    {
        AMSExperiment<T> *exp = NULL;
        int num_unused_exp = 0;

        if (unconverged_experiments.size() > 0 && unused_experiments.size() > 0)
        {
            exp = unconverged_experiments.back();
            unconverged_experiments.pop_back();
            unused_experiments.erase(
                std::find(unused_experiments.begin(),
                          unused_experiments.end(), exp));
        }
        num_unused_exp = unused_experiments.size();

        if (exp == NULL || num_unused_exp == 0)
            continue;

        double max_distance = exp->max_distance;

        int rnd_idx = randint(0, unused_experiments.size()-1);
        while (unused_experiments[rnd_idx]->max_distance <= exp->max_distance)
            rnd_idx = randint(0, unused_experiments.size()-1);

        AMSExperiment<T> *rnd_exp = unused_experiments[rnd_idx];
        if (rnd_exp->dlist.size() == 0)
        {
            std::cerr << "Experiment " << rnd_idx << " has size 0." << std::endl;
            exit(-1);
        }

        int same_distance_idx = -1;
        while (++same_distance_idx < (int)rnd_exp->dlist.size() &&
               rnd_exp->dlist[same_distance_idx] < exp->max_distance);

        if (same_distance_idx == (int)rnd_exp->dlist.size())
        {
            std::cerr << "Distance larger than " << exp->max_distance
                      << " not found in experiment with max distance "
                      << rnd_exp->max_distance << "." << std::endl;
            exit(-1);
        }

        exp->xlist = std::vector<T>(
            rnd_exp->xlist.begin(), rnd_exp->xlist.begin() + same_distance_idx + 1);
        exp->dlist = std::vector<double>(
            rnd_exp->dlist.begin(), rnd_exp->dlist.begin() + same_distance_idx + 1);
        exp->tlist = std::vector<double>(
            rnd_exp->tlist.begin(), rnd_exp->tlist.begin() + same_distance_idx + 1);

        transient_tams(dt_, tmax_, *exp);

        its_++;
        if (exp->converged)
            converged++;
        else
            unconverged_experiments.push_back(exp);

        unused_experiments.push_back(exp);

        std::sort(unconverged_experiments.begin(),
                  unconverged_experiments.end(), AMSExperiment<T>::sort);

        std::cout << "TAMS: " << its_ << " / " << maxit_ << ", "
                  << converged << " / " << num_exp_
                  << " converged with max distance "
                  << max_distance << " -> "
                  << exp->max_distance << " and t="
                  << exp->time << std::endl;

        // Cleanup
        double min_max_distance = 1.0;
        for (auto exp: reactive_experiments)
            min_max_distance = std::min(min_max_distance, exp->max_distance);

        if (its_ % 10 == 0)
        {
            std::cout << "Starting cleanup" << std::endl;
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
            std::cout << "Finished cleanup" << std::endl;
        }

        write_helper(experiments, its_, time_steps_previous_write);
    }
    std::cout << std::endl;

    if (write_final_)
        write(write_, experiments);

    double W = num_exp_ * pow(1.0 - 1.0 / (double)num_exp_, its_);
    for (int i = 1; i < its_; i++)
        W += pow(1.0 - 1.0 / (double)num_exp_, i);

    double tp = (double)converged * pow(1.0 - 1.0 / (double)num_exp_, its_) / W;
    std::cout << "Transition probability T=" << tmax_ << ": " << tp << std::endl;
}

template<class T>
void TimeStepper<T>::gpa(T const &x0) const
{
    std::cout << "Using GPA with an estimated memory usage of: "
              << mem2string((long long)(num_exp_ * vector_length_ * sizeof(double)))
              << std::endl;

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
                    std::cerr << "Particle not found with sum " << sum
                              << " and random value " << val << std::endl;
                    exit(-1);
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

        std::cout << "GPA: " << converged << " / " << num_exp_
                  << " converged with t=" << t << " and eta=" << eta << std::endl;
    }
    std::cout << std::endl;

    double tp = 0.0;
    for (int i = 0; i < num_exp_; i++)
        if (experiments[i].converged)
            tp += experiments[i].probability;
    tp /= (double)num_exp_;

    std::cout << "Transition probability T=" << tmax_ << ": " << tp << std::endl;
}

template<class T>
void TimeStepper<T>::set_random_engine(unsigned int seed)
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
int TimeStepper<T>::randint(int a, int b) const
{
    if (engine_initialized_)
        return randint_(a, b);

    static bool first = true;
    if (first)
    {
        first = false;
        std::cout << "WARNING: Random engine not initialized." << std::endl;
    }

    static thread_local std::random_device rd;
    static thread_local std::seed_seq seeder{rd()};
    static thread_local std::mt19937_64 engine(seeder);
    std::uniform_int_distribution<int> int_distribution(a, b);
    return int_distribution(engine);
}

template<class T>
int TimeStepper<T>::randreal(double a, double b) const
{
    if (engine_initialized_)
        return randreal_(a, b);

    static bool first = true;
    if (first)
    {
        first = false;
        std::cout << "WARNING: Random engine not initialized." << std::endl;
    }

    static thread_local std::random_device rd;
    static thread_local std::seed_seq seeder{rd()};
    static thread_local std::mt19937_64 engine(seeder);
    std::uniform_real_distribution<double> real_distribution(a, b);
    return real_distribution(engine);
}

template<class T>
T TimeStepper<T>::time_step_helper(T const &x, double dt) const
{
    time_steps_++;
    return std::move(time_step_(x, dt));
}

template<class T>
void TimeStepper<T>::write_helper(std::vector<AMSExperiment<T> > const &experiments,
                                  int its, int &time_steps_previous_write) const
{
    if (write_ == "")
        return;

    if (write_steps_ > 0 && its % write_steps_ == 0)
    {
        time_steps_previous_write = time_steps_;
        write(write_, experiments);
        return;
    }

    if (write_time_steps_ > 0 && time_steps_ - time_steps_previous_write >= write_time_steps_)
    {
        time_steps_previous_write = time_steps_;
        write(write_, experiments);
        return;
    }
}

#endif
