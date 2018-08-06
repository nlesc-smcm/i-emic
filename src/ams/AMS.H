#ifndef AMS_H
#define AMS_H

//======================================================================
#include "GlobalDefinitions.H"
#include "Utils.H"
#include "TimeStepper.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"

template<typename Model>
class AMS
{
    using VectorPtr  = typename Model::element_type::VectorPtr;
    using Vector     = typename VectorPtr::element_type;
    using TimeStep   = typename std::function<VectorPtr(VectorPtr const &, double)>;
    using DistFun    = typename std::function<double(VectorPtr const &)>;

    Model model_;
    double dt_;
    double tstep_;
    double tmax_;
    double beta_;
    double theta_;

    double adist_;
    double bdist_;
    double cdist_;
    double dist_tol_;
    double sigma_;

    std::string method_;

    int num_exp_;
    int num_init_exp_;
    int maxit_;

    VectorPtr sol1_;
    VectorPtr sol2_;

    Teuchos::RCP<TimeStepper<VectorPtr> > timestepper_;

    Teuchos::RCP<Epetra_CrsMatrix> B_;
    Teuchos::RCP<Epetra_MultiVector> V_;
    Teuchos::RCP<Epetra_MultiVector> BV_;

    Teuchos::RCP<std::mt19937_64> engine_;
public:
    template<typename ParameterList>
    AMS(Model model, ParameterList pars,
        VectorPtr sol1, VectorPtr sol2);
    ~AMS();

    auto get_time_step();
    auto get_projected_time_step();

    void run();

protected:
    Teuchos::RCP<Epetra_Vector> restrict(Epetra_MultiVector const &x) const;
    Teuchos::RCP<Epetra_Vector> prolongate(Epetra_MultiVector const &x) const;

    unsigned int get_seed();
};

template<typename Vector>
Teuchos::RCP<Vector> newton(
    std::function<Teuchos::RCP<Vector>(Teuchos::RCP<Vector> const &)> F,
    std::function<Teuchos::RCP<Vector>(Teuchos::RCP<Vector> const &,
                                       Teuchos::RCP<Vector> const &)> Jsol,
    Teuchos::RCP<Vector> const &x0,
    double tol)
{
    TIMER_SCOPE("AMS: Newton");
    double nrm = -1;
    Teuchos::RCP<Vector> x = Teuchos::rcp(new Vector(*x0));
    Teuchos::RCP<Vector> Fx = F(x);
    for (int i = 0; i < 20; i++)
    {
        Teuchos::RCP<Vector> dx = Jsol(x, Fx);
        CHECK_ZERO(x->Update(-1.0, *dx, 1.0));
        Fx = F(x);
        nrm = Utils::norm(Fx);
        if (nrm < tol)
            return x;
    }
    std::cerr << "Newton unconverged with norm " << nrm << std::endl;
    return x;
}

// =====================================================================
// Helper methods

Teuchos::RCP<Epetra_MultiVector> dot(
    Epetra_MultiVector const &x, Epetra_MultiVector const &y)
{
    TIMER_SCOPE("AMS: dot");
    int m = x.NumVectors();
    int n = y.NumVectors();

    Epetra_LocalMap map(m, 0, x.Comm());
    Teuchos::RCP<Epetra_MultiVector> out = Teuchos::rcp(new Epetra_MultiVector(map, n));

    CHECK_ZERO(out->Multiply('T', 'N', 1.0, x, y, 0.0));
    return out;
}

template<typename Model>
Teuchos::RCP<Epetra_Vector> AMS<Model>::restrict(Epetra_MultiVector const &x) const
{
    TIMER_SCOPE("AMS: restrict");
    int m = V_->NumVectors();

    Epetra_LocalMap map(m, 0, x.Comm());
    Teuchos::RCP<Epetra_Vector> out = Teuchos::rcp(new Epetra_Vector(map));

    CHECK_ZERO(out->Multiply('T', 'N', 1.0, *V_, x, 0.0));
    return out;
}

template<typename Model>
Teuchos::RCP<Epetra_Vector> AMS<Model>::prolongate(Epetra_MultiVector const &x) const
{
    TIMER_SCOPE("AMS: prolongate");
    VectorPtr out = Teuchos::rcp(new Vector(V_->Map()));
    CHECK_ZERO(out->Multiply('N', 'N', 1.0, *V_, x, 0.0));
    return out;
}

void write_seed(Epetra_Comm const &comm, unsigned int seed, std::string const &label)
{
    unsigned int *seeds = new unsigned int[comm.NumProc()];

    int *seedptr = reinterpret_cast<int *>(&seed);
    int *seedsptr = reinterpret_cast<int *>(seeds);

    CHECK_ZERO(comm.GatherAll(seedptr, seedsptr, 1));
    for (int i = 0; i < comm.NumProc(); i++)
        INFO(label << ": " << seeds[i]);
    delete[] seeds;
}

// =====================================================================
// Class implementation

template<typename Model>
template<typename ParameterList>
AMS<Model>::
AMS(Model model, ParameterList pars,
    VectorPtr sol1, VectorPtr sol2)
    :
    model_(model),
    dt_(pars->get("time step", 0.01)),
    tstep_(pars->get("GPA time step", 1.0)),
    tmax_(pars->get("maximum time", 1000.0)),
    beta_(pars->get("beta", 1.0)),
    theta_(pars->get("theta", 1.0)),
    adist_(pars->get("A distance", 0.05)),
    bdist_(pars->get("B distance", adist_)),
    cdist_(pars->get("C distance", 2 * adist_)),
    dist_tol_(pars->get("distance tolerance", 0.0005)),
    sigma_(pars->get("sigma", 1.0)),
    method_(pars->get("method", "AMS")),
    num_exp_(pars->get("number of experiments", 1000)),
    num_init_exp_(pars->get("number of initial experiments", num_exp_)),
    maxit_(pars->get("maximum iterations", num_exp_ * 10)),
    sol1_(sol1),
    sol2_(sol2),
    V_(Teuchos::null)
{
    std::random_device rd;
    std::random_device::result_type seed = rd();
    std::seed_seq seeder{seed};
    engine_ = Teuchos::rcp(new std::mt19937_64(seeder));

    write_seed(*model_->Comm(), seed, "Time step seed");

    model_->computeForcing();
    B_ = model_->getForcing();

    std::string space = pars->get("space", "");
    TimeStep time_step;
    DistFun dist_fun;
    if (space != "")
    {
        // Use a map that is constructed with the most basic
        // constructor for loading V since
        // MatrixMarketFileToMultiVector does not allow for anything
        // else...
        Epetra_Map map(B_->NumGlobalRows(), 0, *model_->Comm());
        Epetra_MultiVector* Vptr;
        CHECK_ZERO(EpetraExt::MatrixMarketFileToMultiVector(
                       space.c_str(), map, Vptr));

        Epetra_Import import(B_->RangeMap(), map);
        V_ = Teuchos::rcp(new Epetra_MultiVector(B_->RangeMap(), Vptr->NumVectors()));
        CHECK_ZERO(V_->Import(*Vptr, import, Insert));

        delete Vptr;

        BV_ = Teuchos::rcp(new Epetra_MultiVector(B_->DomainMap(), V_->NumVectors()));
        CHECK_ZERO(B_->Multiply('T', *V_, *BV_));

        time_step = get_projected_time_step();

        auto vvec = Teuchos::rcp(new Vector(B_->RangeMap()));
        CHECK_ZERO(vvec->PutScalar(0.0));
        for (int i = 0; i < vvec->MyLength(); i++)
            if (vvec->Map().GID(i) % 6 == 1)
                (*vvec)[i] = 1.0;

        auto Vvvec = Teuchos::rcp(new Epetra_MultiVector(*V_));
        CHECK_ZERO(Vvvec->Multiply(1.0, *vvec, *V_, 0.0));

        auto VvvV = dot(*Vvvec, *Vvvec);

        auto projected_v_norm = [VvvV](VectorPtr const &x) {
            auto tmp = Teuchos::rcp(new Vector(*x));
            CHECK_ZERO(tmp->Multiply('N', 'N', 1.0, *VvvV, *x, 0.0));
            return sqrt((*dot(*x, *tmp))[0][0]);
        };

        sol1_ = restrict(*sol1_);
        sol2_ = restrict(*sol2_);

        auto diff = Teuchos::rcp(new Vector(*sol1_));
        CHECK_ZERO(diff->Update(-1.0, *sol2_, 1.0));
        double nrm = projected_v_norm(diff);

        dist_fun = [nrm, projected_v_norm, this](VectorPtr const &x) {
            VectorPtr d1v = Teuchos::rcp(new Vector(*x));
            CHECK_ZERO(d1v->Update(-1.0, *sol1_, 1.0));
            double d1 =  projected_v_norm(d1v) / nrm;
            VectorPtr d2v = Teuchos::rcp(new Vector(*x));
            CHECK_ZERO(d2v->Update(-1.0, *sol2_, 1.0));
            double d2 =  projected_v_norm(d2v) / nrm;
            double dist = 0.5 - 0.5 * exp(-0.5 * pow(d1 / 0.25, 2.))
            + 0.5 * exp(-0.5 * pow(d2 / 0.25, 2.));
            INFO("distance = " << dist);
            return dist;
        };
    }
    else
    {
        time_step = get_time_step();

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
            {
                CHECK_ZERO(x->Comm().SumAll(&tmp, &out, 1));
            }
            else
            {
                out = tmp;
            }
            return sqrt(out);
        };

        double nrm[6];
        for (int var = 0; var < 6; var++)
        {
            VectorPtr diff = Teuchos::rcp(new Vector(*sol1_));
            CHECK_ZERO(diff->Update(-1.0, *sol2_, 1.0));
            nrm[var] = var_norm(diff, var);
        }

        auto dist_fun_var = [nrm, var_norm, this](VectorPtr const &x, int var) {
            VectorPtr d1v = Teuchos::rcp(new Vector(*x));
            CHECK_ZERO(d1v->Update(-1.0, *sol1_, 1.0));
            double d1 =  var_norm(d1v, var) / nrm[var];
            VectorPtr d2v = Teuchos::rcp(new Vector(*x));
            CHECK_ZERO(d2v->Update(-1.0, *sol2_, 1.0));
            double d2 =  var_norm(d2v, var) / nrm[var];
            double dist = 0.5 - 0.5 * exp(-0.5 * pow(d1 / 0.25, 2.))
            + 0.5 * exp(-0.5 * pow(d2 / 0.25, 2.));
            return dist;
        };

        dist_fun = [this, dist_fun_var](VectorPtr const &x) {
            double dist = dist_fun_var(x, 1);
            INFO("distance = " << dist);
            return dist;
        };
    }

    INFO("Time step: " << dt_);

    timestepper_ = Teuchos::rcp(new TimeStepper<VectorPtr>(
                                    time_step, dist_fun, adist_, bdist_, cdist_, dist_tol_));

    timestepper_->set_random_engine(get_seed());
}

template<typename Model>
auto AMS<Model>::get_time_step()
{
    // Function to perform one stochastic time step
    return [this](VectorPtr const &x, double dt) {
        TIMER_SCOPE("AMS: Time step");
        double sig = 1.0 / (dt * theta_);

        *model_->getState('V') = *x;
        // model_->saveState();

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

        VectorPtr G = model_->getState('C');
        CHECK_ZERO(B_->Apply(pert, *G));
        CHECK_ZERO(G->Scale(sqrt(dt_) * sigma_));

        // Deterministic theta stepper:
        // M * u_n + dt * theta * F(u_(n+1)) + (1-theta) * F(u_n) - M * u_(n+1) = 0
        // Noise is added in an explicit manner
        model_->computeRHS();
        VectorPtr Fx = model_->getRHS('C');
        // INFO("||F|| = " << Utils::norm(Fx));
        auto thetaF = [this, G, Fx, x, dt](VectorPtr const &xnew) {
            TIMER_SCOPE("AMS: F");
            *model_->getState('V') = *xnew;

            VectorPtr M = model_->getMassMat('V');
            VectorPtr ret = Teuchos::rcp(new Vector(*x));
            CHECK_ZERO(ret->Update(-1.0, *xnew, 1.0));
            CHECK_ZERO(ret->Multiply(1.0, *ret, *M, 0.0));

            model_->computeRHS();
            VectorPtr Fxnew = model_->getRHS('V');
            CHECK_ZERO(ret->Update(dt * theta_, *Fxnew, 1.0));
            CHECK_ZERO(ret->Update(dt * (1.0 - theta_), *Fx, 1.0));
            CHECK_ZERO(ret->Update(1.0, *G, 1.0));
            return ret;
        };

        // J2 = theta * dt * J - M, so J2*x = b
        // We write this as
        // J2 = J - 1/(theta*dt) * M, J2 * x = 1/(theta*dt) * b
        auto thetaJ = [this, sig](VectorPtr const &xnew, VectorPtr const &b) {
            TIMER_SCOPE("AMS: Jacobian solve");
            *model_->getState('V') = *xnew;
            VectorPtr y = Teuchos::rcp(new Vector(*b));
            CHECK_ZERO(y->Scale(sig));
            model_->setShift(-sig);
            model_->computeJacobian();
            model_->solve(y);
            return model_->getSolution('V');
        };

        return newton<Vector>(thetaF, thetaJ, x, 1e-8);
    };
}

template<typename Model>
auto AMS<Model>::get_projected_time_step()
{
    // Function to perform one stochastic time step
    return [this](VectorPtr const &x, double dt) {
        TIMER_SCOPE("AMS: Projected time step");
        double sig = 1.0 / (dt * theta_);

        *model_->getState('V') = *prolongate(*x);
        // model_->saveState();

        TIMER_START("AMS: Compute noise");
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

        Epetra_LocalMap Gmap(BV_->NumVectors(), 0, *model_->Comm());
        auto G = Teuchos::rcp(new Epetra_MultiVector(Gmap, 1));
        CHECK_ZERO(G->Multiply('T', 'N', 1.0, *BV_, pert, 0.0));
        CHECK_ZERO(G->Scale(sqrt(dt_) * sigma_));
        TIMER_STOP("AMS: Compute noise");

        // Deterministic theta stepper:
        // M * u_n + dt * theta * F(u_(n+1)) + (1-theta) * F(u_n) - M * u_(n+1) = 0
        // Noise is added in an explicit manner
        model_->computeRHS();
        auto Fx = restrict(*model_->getRHS('V'));
        // INFO("||F|| = " << Utils::norm(Fx));

        auto tmp = Teuchos::rcp(new Epetra_MultiVector(*V_));

        TIMER_START("AMS: Compute VMV");
        auto M = model_->getMassMat('V');
        CHECK_ZERO(tmp->Multiply(1.0, *M, *V_, 0.0));
        auto VMV = dot(*V_, *tmp);
        TIMER_STOP("AMS: Compute VMV");

        auto thetaF = [this, &G, &Fx, &VMV, &x, dt](VectorPtr const &xnew) {
            TIMER_SCOPE("AMS: Projected F");
            *model_->getState('V') = *prolongate(*xnew);

            auto ret = Teuchos::rcp(new Vector(*x));
            CHECK_ZERO(ret->Update(-1.0, *xnew, 1.0));
            auto tmp = Teuchos::rcp(new Vector(*ret));
            CHECK_ZERO(ret->Multiply('N', 'N', 1.0, *VMV, *tmp, 0.0));

            model_->computeRHS();
            auto Fxnew = restrict(*model_->getRHS('V'));
            CHECK_ZERO(ret->Update(dt * theta_, *Fxnew, 1.0));
            CHECK_ZERO(ret->Update(dt * (1.0 - theta_), *Fx, 1.0));
            CHECK_ZERO(ret->Update(1.0, *G, 1.0));
            return ret;
        };

        // J2 = theta * dt * J - M, so J2*x = b
        // We write this as
        // J2 = J - 1/(theta*dt) * M, J2 * x = 1/(theta*dt) * b
        TIMER_START("AMS: Compute VAV");
        model_->setShift(-sig);
        model_->computeJacobian();
        model_->applyMatrix(*V_, *tmp);
        auto VAV = dot(*V_, *tmp);
        auto VAVmat = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                       Copy, VAV->Values(), VAV->Stride(),
                                       VAV->MyLength(), VAV->NumVectors()));
        auto VAVsolver = Teuchos::rcp(new Epetra_SerialDenseSolver());
        CHECK_ZERO(VAVsolver->SetMatrix(*VAVmat));
        CHECK_ZERO(VAVsolver->Factor());
        TIMER_STOP("AMS: Compute VAV");

        auto thetaJ = [this, sig, &VAVsolver](VectorPtr const &xnew,
                                              VectorPtr const &b) {
            TIMER_SCOPE("AMS: Projected jacobian solve");
            auto y = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                      Copy, b->Values(), b->Stride(),
                                      b->MyLength(), b->NumVectors()));
            auto z = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                      Copy, b->Values(), b->Stride(),
                                      b->MyLength(), b->NumVectors()));
            CHECK_ZERO(z->Scale(sig));
            CHECK_ZERO(VAVsolver->SetVectors(*y, *z));
            CHECK_NONNEG(VAVsolver->Solve());
            return Teuchos::rcp(new Vector(
                                    Copy, Epetra_LocalMap(y->M(), 0, *model_->Comm()),
                                    y->A()));
        };

        return newton<Vector>(thetaF, thetaJ, x, 1e-8);
    };
}

template<typename Model>
AMS<Model>::~AMS()
{
    INFO("AMS destructor");
}

template<typename Model>
void AMS<Model>::run()
{
    if (method_ == "AMS")
        timestepper_->ams(num_exp_, num_init_exp_, sol1_, dt_, tmax_);
    else if (method_ == "TAMS")
        timestepper_->tams(num_exp_, maxit_, sol1_, dt_, tmax_);
    else if (method_ == "GPA")
        timestepper_->gpa(num_exp_, maxit_, beta_, sol1_, dt_, tstep_, tmax_);
    else if (method_ == "Naive")
        timestepper_->naive(num_exp_, sol1_, dt_, tmax_);
    else if (method_ == "Brownian motion")
        timestepper_->brownian_motion(num_exp_, dt_, tmax_);
    else
        std::cerr << "Method " << method_ << " does not exist." << std::endl;
}

template<typename Model>
unsigned int AMS<Model>::get_seed()
{
    static thread_local std::random_device rd;
    unsigned int seed = rd();

    int *seed_ptr = reinterpret_cast<int *>(&seed);
    CHECK_ZERO(model_->Comm()->Broadcast(seed_ptr, 1, 0));

    write_seed(*model_->Comm(), seed, "Global seed");

    return seed;
}

#endif
