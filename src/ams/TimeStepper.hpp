#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <iostream>

template<class T>
struct AMSExperiment {
    std::vector<T> xlist;
    std::vector<double> dlist;
    std::vector<double> tlist;
    double max_distance;
    double time;
    double initial_time;
    double return_time;
    bool converged;

    static bool sort(AMSExperiment<T> const *e1, AMSExperiment<T> const *e2)
        {
            return e1->max_distance > e2->max_distance;
        };
};

template<class T>
class TimeStepper
{
    std::function<T(T const &, double)> time_step_;
    std::function<double(T const &)> dist_fun_;

    double adist_;
    double bdist_;
    double cdist_;

    // RNG methods
    bool engine_initialized_;
    std::function<int(int, int)> randint_;

public:
    TimeStepper(std::function<T(T const &, double)> time_step);
    TimeStepper(std::function<T(T const &, double)> time_step,
                std::function<double(T const &)> dist_fun,
                double rho);
    TimeStepper(std::function<T(T const &, double)> time_step,
                std::function<double(T const &)> dist_fun,
                double adist,
                double bdist,
                double cdist);

    T transient(T x, double dt, double tmax) const;
    double transient_max_distance(
        T x, double dt, double tmax, double max_distance) const;

    void transient_start(
        T const &x0, double dt, double tmax,
        AMSExperiment<T> &experiment) const;

    void transient_ams(
        double dt, double tmax,
        AMSExperiment<T> &experiment) const;

    void ams(int num_exp, int num_init_exp,
             T const &x0, double dt, double tmax) const;

    template<class URNG>
    void set_random_engine(URNG &engine);

protected:
    int randint(int a, int b) const;
};

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &, double)> time_step)
    :
    time_step_(time_step),
    engine_initialized_(false)
{}

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &, double)> time_step,
                            std::function<double(T const &)> dist_fun,
                            double rho)
    :
    TimeStepper(time_step, dist_fun, rho, rho, 2 * rho)
{}

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &, double)> time_step,
                            std::function<double(T const &)> dist_fun,
                            double adist,
                            double bdist,
                            double cdist)
    :
    time_step_(time_step),
    dist_fun_(dist_fun),
    adist_(adist),
    bdist_(bdist),
    cdist_(cdist),
    engine_initialized_(false)
{
}

template<class T>
T TimeStepper<T>::transient(T x, double dt, double tmax) const
{
    for (double t = 0; t < tmax; t += dt)
        x = std::move(time_step_(x, dt));
    return x;
}

template<class T>
double TimeStepper<T>::transient_max_distance(
    T x, double dt, double tmax, double max_distance) const
{
    const double lim = max_distance - bdist_;
    for (double t = 0; t < tmax; t += dt)
    {
        x = std::move(time_step_(x, dt));
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

    for (double t = 0; t < tmax; t += dt)
    {
        x = std::move(time_step_(x, dt));

        double dist = dist_fun_(x);
        if (dist > cdist_)
        {
            experiment.xlist.push_back(x);
            experiment.dlist.push_back(dist);
            experiment.tlist.push_back(0);
            experiment.max_distance = dist;
            experiment.initial_time = t;
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
    double t = experiment.tlist.back();
    double tend = t + tmax;

    double max_distance = experiment.max_distance;

    for (; t < tend; t += dt)
    {
        x = std::move(time_step_(x, dt));

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
        if (dist > max_distance)
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
void TimeStepper<T>::ams(int num_exp, int num_init_exp,
                         T const &x0, double dt, double tmax) const
{

    int its = 0;
    int converged = 0;
    std::vector<AMSExperiment<T>> experiments(num_init_exp);

#pragma omp parallel for default(none),                                 \
    firstprivate(num_exp, num_init_exp, dt, tmax),                      \
    shared(std::cout, its, experiments, converged, x0), schedule(dynamic)
    for (int i = 0; i < num_init_exp; i++)
    {
        transient_start(x0, dt, tmax,
                        experiments[i]);

        transient_ams(
            dt, tmax, experiments[i]);

        // Erase data that we do not need for later experiments
        if (i >= num_exp)
        {
            experiments[i].xlist = std::vector<T>();
            experiments[i].dlist = std::vector<double>();
            experiments[i].tlist = std::vector<double>();
        }

        if (experiments[i].converged)
#pragma omp atomic
            converged++;

#pragma omp atomic
        its++;

        std::cout << "Initialization: " << its << " / " << num_init_exp << ", "
                  << converged << " / " << num_init_exp
                  << " converged with t="
                  << experiments[i].initial_time + experiments[i].time << std::endl;
    }
    std::cout << std::endl;

    int maxit = 10000000;
    its = 0;
    converged = 0;

    std::vector<AMSExperiment<T> *> unconverged_experiments;
    std::vector<AMSExperiment<T> *> unused_experiments;
    std::vector<AMSExperiment<T> *> reactive_experiments;

    for (int i = 0; i < num_exp; i++)
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

#pragma omp parallel for default(none),                                 \
    firstprivate(num_exp, maxit, dt, tmax),                             \
    shared(std::cout, std::cerr, its, converged,                        \
           experiments, unused_experiments, unconverged_experiments,    \
           reactive_experiments),                                       \
    schedule(static)
    for (int i = 0; i < maxit; i++)
    {
        AMSExperiment<T> *exp = NULL;

#pragma omp critical (experiments)
        {
#pragma omp flush (unconverged_experiments, unused_experiments)
            if (unconverged_experiments.size() > 0 && unused_experiments.size() > 0)
            {
                exp = unconverged_experiments.back();
                unconverged_experiments.pop_back();
                unused_experiments.erase(
                    std::find(unused_experiments.begin(),
                              unused_experiments.end(), exp));
            }
        }

        if (exp == NULL)
            continue;

        double max_distance = exp->max_distance;

#pragma omp critical (experiments)
        {
#pragma omp flush (unused_experiments)
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
            while (++same_distance_idx < rnd_exp->dlist.size() &&
                   rnd_exp->dlist[same_distance_idx] < exp->max_distance);

            if (same_distance_idx == rnd_exp->dlist.size())
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
        }

        transient_ams(dt, tmax, *exp);

#pragma omp critical (experiments)
        {
#pragma omp flush (experiments, unconverged_experiments, unused_experiments, \
                   its, converged)
            its++;
            if (exp->converged)
                converged++;
            else
                unconverged_experiments.push_back(exp);

            unused_experiments.push_back(exp);

            std::sort(unconverged_experiments.begin(),
                      unconverged_experiments.end(), AMSExperiment<T>::sort);

            std::cout << "AMS: " << its << " / " << maxit << ", "
                      << converged << " / " << num_exp
                      << " converged with max distance "
                      << max_distance << " -> "
                      << exp->max_distance << " and t="
                      << exp->initial_time + exp->time << std::endl;

            // Cleanup
            double min_max_distance = 1.0;
            for (auto exp: reactive_experiments)
                min_max_distance = std::min(min_max_distance, exp->max_distance);

            if (its % 10 == 0)
            {
                std::cout << "Starting cleanup" << std::endl;
                for (auto exp: unused_experiments)
                {
                    int min_max_idx = -1;
                    while (++min_max_idx < exp->dlist.size() &&
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
        }
    }
    std::cout << std::endl;

    double total_tr = 0;
    double total_t1 = 0;
    double total_t2 = 0;
    int num_t1 = num_init_exp;
    int num_t2 = 0;
    for (auto exp: reactive_experiments)
        total_tr += exp->time;

    for (auto &exp: experiments)
    {
        total_t1 += exp.initial_time;
        total_t2 += exp.return_time;
        if (exp.return_time > dt / 2.0)
            num_t2++;
    }

    double alpha = (double)converged / (double)num_exp * pow((1.0 - 1.0 / (double)num_exp), its);

    double meann = 1.0 / alpha - 1.0;
    double fpt = meann * (total_t1 / (double)num_t1 + total_t2 / (double)num_t2) +
        total_t1 / (double)num_t1 + total_tr / (double)converged;
    std::cout << "Alpha: " << alpha << std::endl;
    std::cout << "Mean first passage time: " << fpt << std::endl;

    tmax = 10;
    double tp = 1.0 -  exp(-1.0 / fpt * tmax);
    std::cout << "Transition probability T=" << tmax << ": " << tp << std::endl;
}

template<class T>
template<class URNG>
void TimeStepper<T>::set_random_engine(URNG &engine)
{
    randint_ = [&engine](int a, int b) {
        static thread_local std::uniform_int_distribution<int> int_distribution(a, b);
        return int_distribution(engine);
    };
    engine_initialized_ = true;
}

template<class T>
int TimeStepper<T>::randint(int a, int b) const
{
    if (engine_initialized_)
        return randint_(a, b);

    static thread_local std::random_device rd;
    static thread_local std::default_random_engine engine(rd());
    static thread_local std::uniform_int_distribution<int> int_distribution(a, b);
    return int_distribution(engine);
}

#endif
