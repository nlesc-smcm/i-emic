#ifndef AMS_H
#define AMS_H

//======================================================================
#include "GlobalDefinitions.H"
#include "Utils.H"
#include "TimeStepper.hpp"

#include "EpetraExt_CrsMatrixIn.h"

template<typename Model>
class AMS
{
    using VectorPtr  = typename Model::element_type::VectorPtr;
    using Vector     = typename VectorPtr::element_type;

    Model model_;
    double dt_;
    double tmax_;
    double adist_;
    double bdist_;
    double cdist_;
    double sigma_;

    int num_exp_;
    int num_init_exp_;

    VectorPtr sol1_;
    VectorPtr sol2_;

    Teuchos::RCP<TimeStepper<VectorPtr> > timestepper_;


public:
    template<typename ParameterList>
    AMS(Model model, ParameterList pars,
        VectorPtr sol1, VectorPtr sol2);
    ~AMS();

    void run();
};

template<typename VectorPtr>
VectorPtr newton(
    std::function<VectorPtr(VectorPtr const &)> F,
    std::function<VectorPtr(VectorPtr const &, VectorPtr const &)> Jsol,
    VectorPtr const &x0,
    double tol)
{
    double nrm = -1;
    VectorPtr x = Teuchos::rcp(new Epetra_Vector(*x0));
    VectorPtr Fx = F(x);
    for (int i = 0; i < 20; i++)
    {
        VectorPtr dx = Jsol(x, Fx);
        x->Update(-1.0, *dx, 1.0);
        Fx = F(x);
        nrm = Utils::norm(Fx);
        if (nrm < tol)
            return x;
    }
    std::cerr << "Newton unconverged with norm " << nrm << std::endl;
    return x;
}

//======================================================================
//Constructor
template<typename Model>
template<typename ParameterList>
AMS<Model>::
AMS(Model model, ParameterList pars,
    VectorPtr sol1, VectorPtr sol2)
    :
    model_(model),
    dt_(pars->get("time step", 0.01)),
    tmax_(pars->get("maximum time", 1000.0)),
    adist_(pars->get("A distance", 0.05)),
    bdist_(pars->get("B distance", adist_)),
    cdist_(pars->get("C distance", 2 * adist_)),
    sigma_(pars->get("sigma", 1.0)),
    num_exp_(pars->get("number of experiments", 1000)),
    num_init_exp_(pars->get("number of initial experiments", num_exp_)),
    sol1_(sol1),
    sol2_(sol2)
{
    Epetra_CrsMatrix* Bptr;
    CHECK_ZERO(EpetraExt::MatrixMarketFileToCrsMatrix(
                   "B.mtx", *model_->getDomain()->GetSolveMap(), Bptr));
    Teuchos::RCP<Epetra_CrsMatrix> B = Teuchos::rcp(Bptr);

    INFO("Time step: " << dt_);

    // Function to perform one stochastic time step
    auto time_step = [this, B](VectorPtr const &x, double dt) {
        double theta = 1;
        double sig = 1.0 / (dt * theta);
        *model_->getState('V') = *x;

        // FIXME!!!
        // VectorPtr G = B(x);

        // Compute noise for forcing which is independent per processor
        static thread_local std::random_device rd;
        static thread_local std::default_random_engine engine(rd());
        static thread_local std::normal_distribution<double> distribution(0.0, 1.0);
        static thread_local auto generator = std::bind(distribution, std::ref(engine));

        int m = B->NumMyCols();
        if (!B->ColMap().UniqueGIDs())
        {
            ERROR("The values of B are distributed", __FILE__, __LINE__);
        }

        Epetra_Vector pert(B->ColMap());
        std::generate(pert.Values(), pert.Values() + m, generator);
        VectorPtr G = model_->getState('C');
        B->Apply(pert, *G);
        G->Scale(sqrt(dt_) * sigma_);

        // Deterministic theta stepper:
        // M * u_n + dt * theta * F(u_(n+1)) + (1-theta) * F(u_n) - M * u_(n+1) = 0
        // Noise is added in an explicit manner
        model_->computeRHS();
        VectorPtr Fx = model_->getRHS('C');
        INFO("||F|| = " << Utils::norm(Fx));
        auto thetaF = [this, G, Fx, x, dt, theta](VectorPtr const &xnew) {
            *model_->getState('V') = *xnew;

            VectorPtr M = model_->getMassMat('V');
            VectorPtr ret = Teuchos::rcp(new Vector(*x));
            ret->Update(-1.0, *xnew, 1.0);
            ret->Multiply(1.0, *ret, *M, 0.0);

            model_->computeRHS();
            VectorPtr Fxnew = model_->getRHS('V');
            ret->Update(dt * theta, *Fxnew, 1.0);
            ret->Update(dt * (1.0 - theta), *Fx, 1.0);
            ret->Update(1.0, *G, 1.0);
            return ret;
        };

        // J2 = theta * dt * J - M, so J2*x = b
        // We write this as
        // J2 = J - 1/(theta*dt) * M, J2 * x = 1/(theta*dt) * b
        auto thetaJ = [this, sig](VectorPtr const &xnew, VectorPtr const &b) {
            *model_->getState('V') = *xnew;
            VectorPtr y = Teuchos::rcp(new Vector(*b));
            y->Scale(sig);
            model_->setShift(-sig);
            model_->computeJacobian();
            model_->solve(y);
            return model_->getSolution('V');
        };

        return newton<VectorPtr>(thetaF, thetaJ, x, 1e-8);
    };

    auto var_norm = [](VectorPtr const &x, int var) {
        double tmp = 0;
        int n = x->MyLength();
        auto &map = x->Map();
        Vector const &xref = *x;
        for (int i = 0; i < n; i++)
            if (map.GID(i) % 6 == var)
                tmp += xref[i] * xref[i];
        double out = 0;
        if (x->DistributedGlobal())
            x->Comm().SumAll(&tmp, &out, 1);
        else
            out = tmp;
        return sqrt(out);
    };

    double nrm[6];
    for (int var = 0; var < 6; var++)
    {
        VectorPtr diff = Teuchos::rcp(new Vector(*sol1));
        diff->Update(-1.0, *sol2, 1.0);
        nrm[var] = var_norm(diff, var);
    }

    auto dist_fun_var = [nrm, var_norm, sol1, sol2](VectorPtr const &x, int var) {
        VectorPtr d1v = Teuchos::rcp(new Vector(*x));
        d1v->Update(-1.0, *sol1, 1.0);
        double d1 =  var_norm(d1v, var) / nrm[var];
        VectorPtr d2v = Teuchos::rcp(new Vector(*x));
        d2v->Update(-1.0, *sol2, 1.0);
        double d2 =  var_norm(d2v, var) / nrm[var];
        double dist = 0.25 - 0.25 * exp(-0.5 * pow(d1 / 0.25, 2.))
        + 0.75 * exp(-0.5 * pow(d2 / 0.25, 2.));
        return dist;
    };

    auto dist_fun = [this, dist_fun_var](VectorPtr const &x) {
        // BEGIN DEBUG
        for (int var = 0; var < 6; var++)
        {
            double dist = dist_fun_var(x, var);
            INFO("distance " << var << " = " << dist);
        }

        INFO(model_->writeData(true));
        INFO(model_->writeData(false));
        // END DEBUG

        double dist = dist_fun_var(x, 1);
        INFO("distance = " << dist);
        return dist;
    };

    timestepper_ = Teuchos::rcp(new TimeStepper<VectorPtr>(
                                    time_step, dist_fun, adist_, bdist_, cdist_));
}

//======================================================================
//Destructor
template<typename Model>
AMS<Model>::~AMS()
{
    INFO("AMS destructor");
}

template<typename Model>
void AMS<Model>::run()
{
    timestepper_->ams(num_exp_, num_init_exp_, sol1_, dt_, tmax_);
}

#endif
