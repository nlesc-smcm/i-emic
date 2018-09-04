#ifndef TIMESTEPPERDECL_HPP
#define TIMESTEPPERDECL_HPP

#include <random>
#include <functional>

template<class T>
class AMSExperiment;

template<class T>
class GPAExperiment;

template<class T>
class TimeStepper
{
    std::function<T(T const &, double)> time_step_;
    std::function<double(T const &)> dist_fun_;

    int vector_length_;

    double dt_;
    double tstep_;
    double tmax_;
    double beta_;

    double adist_;
    double bdist_;
    double cdist_;
    double dist_tol_;

    int num_exp_;
    int num_init_exp_;

    int maxit_;
    mutable int its_;
    mutable int time_steps_;

    std::string read_;
    std::string write_;
    bool write_final_;
    int write_steps_;
    int write_time_steps_;

    // RNG methods
    bool engine_initialized_;
    std::function<int(int, int)> randint_;
    std::function<double(double, double)> randreal_;

    // Random engine. FIXME: Does not work with omp!!!
    std::mt19937_64 *engine_;

public:
    TimeStepper(std::function<T(T const &, double)> time_step);
    TimeStepper(std::function<T(T const &, double)> time_step,
                std::function<double(T const &)> dist_fun,
                int vector_length);

   virtual ~TimeStepper();

    template<class ParameterList>
    void set_parameters(ParameterList &params);

    T transient(T x, double dt, double tmax) const;

    double transient_max_distance(
        T x, double dt, double tmax, double max_distance) const;

    void transient_start(
        T const &x0, double dt, double tmax,
        AMSExperiment<T> &experiment) const;

    void transient_ams(
        double dt, double tmax,
        AMSExperiment<T> &experiment) const;

    void transient_tams(
        double dt, double tmax,
        AMSExperiment<T> &experiment) const;

    void transient_gpa(
        double dt, double tmax,
        GPAExperiment<T> &experiment) const;

    void naive(T const &x0) const;

    void ams(T const &x0) const;

    void tams(T const &x0) const;

    void gpa(T const &x0) const;

    void read(
        std::string const &name,
        std::vector<AMSExperiment<T> > &experiments) const;

    void write(
        std::string const &name,
        std::vector<AMSExperiment<T> > const &experiments) const;

    void set_random_engine(unsigned int seed);

protected:
    int randint(int a, int b) const;
    int randreal(double a, double b) const;

    T time_step_helper(T const &x, double dt) const;

    void write_helper(std::vector<AMSExperiment<T> > const &experiments,
                      int its, int &time_steps_previous_write) const;
};

#endif
