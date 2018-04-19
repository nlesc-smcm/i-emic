#ifndef AMS_H
#define AMS_H

//======================================================================
#include "GlobalDefinitions.H"
#include "Utils.H"
#include "TimeStepper.hpp"

#include "Epetra_SerialComm.h"
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

    Teuchos::RCP<Epetra_CrsMatrix> B_;
    Teuchos::RCP<Epetra_MultiVector> V_;
    Teuchos::RCP<Epetra_MultiVector> BV_;
public:
    template<typename ParameterList>
    AMS(Model model, ParameterList pars,
        VectorPtr sol1, VectorPtr sol2);
    ~AMS();

    auto get_time_step();
    auto get_projected_time_step();

    void run();

protected:
    Teuchos::RCP<Epetra_MultiVector> restrict(Epetra_MultiVector const &x) const;
    Teuchos::RCP<Epetra_Vector> prolongate(Epetra_MultiVector const &x) const;
};

template<typename Vector>
Teuchos::RCP<Vector> newton(
    std::function<Teuchos::RCP<Vector>(Teuchos::RCP<Vector> const &)> F,
    std::function<Teuchos::RCP<Vector>(Teuchos::RCP<Vector> const &,
                                       Teuchos::RCP<Vector> const &)> Jsol,
    Teuchos::RCP<Vector> const &x0,
    double tol)
{
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

Teuchos::RCP<Epetra_MultiVector> MakeLocalVector(int m, int n)
{
    Epetra_SerialComm comm;
    Epetra_LocalMap map(m, 0, comm);
    return Teuchos::rcp(new Epetra_MultiVector(map, n));
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
    tmax_(pars->get("maximum time", 1000.0)),
    adist_(pars->get("A distance", 0.05)),
    bdist_(pars->get("B distance", adist_)),
    cdist_(pars->get("C distance", 2 * adist_)),
    sigma_(pars->get("sigma", 1.0)),
    num_exp_(pars->get("number of experiments", 1000)),
    num_init_exp_(pars->get("number of initial experiments", num_exp_)),
    sol1_(sol1),
    sol2_(sol2),
    V_(Teuchos::null)
{
    Epetra_CrsMatrix* Bptr;
    CHECK_ZERO(EpetraExt::MatrixMarketFileToCrsMatrix(
                   "B.mtx", *model_->getDomain()->GetSolveMap(), Bptr));
    B_ = Teuchos::rcp(Bptr);

    std::string space = pars->get("space", "");
    TimeStep time_step;
    if (space != "")
    {
        Epetra_MultiVector* Vptr;
        CHECK_ZERO(EpetraExt::MatrixMarketFileToMultiVector(
                       space.c_str(), *model_->getDomain()->GetSolveMap(), Vptr));
        V_ = Teuchos::rcp(Vptr);
        
        BV_ = MakeLocalVector(B_->NumGlobalCols(), V_->NumVectors());
        CHECK_ZERO(B_->Multiply('T', *V_, *BV_));

        time_step = get_projected_time_step();
    }
    else
        time_step = get_time_step();

    INFO("Time step: " << dt_);

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

    static thread_local std::random_device rd;
    unsigned int seed = rd();

    int *seed_ptr = reinterpret_cast<int *>(&seed);
    CHECK_ZERO(model->Comm()->Broadcast(seed_ptr, 1, 0));

    static thread_local std::default_random_engine engine(seed);
    timestepper_->set_random_engine(engine);
}

template<typename Model>
auto AMS<Model>::get_time_step()
{
    // Function to perform one stochastic time step
    return [this](VectorPtr const &x, double dt) {
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
        INFO("||F|| = " << Utils::norm(Fx));
        auto thetaF = [this, G, Fx, x, dt, theta](VectorPtr const &xnew) {
            *model_->getState('V') = *xnew;

            VectorPtr M = model_->getMassMat('V');
            VectorPtr ret = Teuchos::rcp(new Vector(*x));
            CHECK_ZERO(ret->Update(-1.0, *xnew, 1.0));
            CHECK_ZERO(ret->Multiply(1.0, *ret, *M, 0.0));

            model_->computeRHS();
            VectorPtr Fxnew = model_->getRHS('V');
            CHECK_ZERO(ret->Update(dt * theta, *Fxnew, 1.0));
            CHECK_ZERO(ret->Update(dt * (1.0 - theta), *Fx, 1.0));
            CHECK_ZERO(ret->Update(1.0, *G, 1.0));
            return ret;
        };

        // J2 = theta * dt * J - M, so J2*x = b
        // We write this as
        // J2 = J - 1/(theta*dt) * M, J2 * x = 1/(theta*dt) * b
        auto thetaJ = [this, sig](VectorPtr const &xnew, VectorPtr const &b) {
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

Teuchos::RCP<Epetra_MultiVector> dot(
    Epetra_MultiVector const &x, Epetra_MultiVector const &y)
{
    int m = x.NumVectors();
    int n = y.NumVectors();
    Teuchos::RCP<Epetra_MultiVector> out = MakeLocalVector(m, n);

    CHECK_ZERO(out->Multiply('T', 'N', 1.0, x, y, 0.0));
    return out;
}

template<typename Model>
Teuchos::RCP<Epetra_MultiVector> AMS<Model>::restrict(Epetra_MultiVector const &x) const
{
    return dot(*V_, x);
}

template<typename Model>
Teuchos::RCP<Epetra_Vector> AMS<Model>::prolongate(Epetra_MultiVector const &x) const
{
    VectorPtr out = Teuchos::rcp(new Vector(V_->Map()));
    CHECK_ZERO(out->Multiply('N', 'N', 1.0, *V_, x, 0.0));
    return out;
}

template<typename Model>
auto AMS<Model>::get_projected_time_step()
{
    // Function to perform one stochastic time step
    return [this](VectorPtr const &x, double dt) {
        double theta = 1;
        double sig = 1.0 / (dt * theta);

        auto xsmall = restrict(*x);
        *model_->getState('V') = *x;

        using SmallVector = Epetra_MultiVector;
        using SmallVectorPtr = Teuchos::RCP<Epetra_MultiVector>;

        // FIXME!!!
        // VectorPtr G = B(x);

        // Compute noise for forcing which is independent per processor
        static thread_local std::random_device rd;
        static thread_local std::default_random_engine engine(rd());
        static thread_local std::normal_distribution<double> distribution(0.0, 1.0);
        static thread_local auto generator = std::bind(distribution, std::ref(engine));

        int m = BV_->MyLength();
        auto pert = MakeLocalVector(m, 1);
        std::generate(pert->Values(), pert->Values() + m, generator);
        auto G = MakeLocalVector(BV_->NumVectors(), 1);
        CHECK_ZERO(G->Multiply('T', 'N', 1.0, *BV_, *pert, 0.0));
        CHECK_ZERO(G->Scale(sqrt(dt_) * sigma_));

        // Deterministic theta stepper:
        // M * u_n + dt * theta * F(u_(n+1)) + (1-theta) * F(u_n) - M * u_(n+1) = 0
        // Noise is added in an explicit manner
        model_->computeRHS();
        auto Fx = restrict(*model_->getRHS('V'));
        INFO("||F|| = " << Utils::norm(Fx));

        auto tmp = Teuchos::rcp(new Epetra_MultiVector(*V_));

        auto M = model_->getDiagB('V');
        CHECK_ZERO(tmp->Multiply(1.0, *M, *V_, 0.0));
        auto VMV = dot(*V_, *tmp);

        auto thetaF = [this, &G, &Fx, &VMV, &xsmall, dt, theta](SmallVectorPtr const &xnew) {
            *model_->getState('V') = *prolongate(*xnew);

            auto ret = Teuchos::rcp(new SmallVector(*xsmall));
            CHECK_ZERO(ret->Update(-1.0, *xnew, 1.0));
            CHECK_ZERO(ret->Multiply('N', 'N', 1.0, *VMV, *ret, 0.0));

            model_->computeRHS();
            auto Fxnew = restrict(*model_->getRHS('V'));
            CHECK_ZERO(ret->Update(dt * theta, *Fxnew, 1.0));
            CHECK_ZERO(ret->Update(dt * (1.0 - theta), *Fx, 1.0));
            CHECK_ZERO(ret->Update(1.0, *G, 1.0));
            return ret;
        };

        // J2 = theta * dt * J - M, so J2*x = b
        // We write this as
        // J2 = J - 1/(theta*dt) * M, J2 * x = 1/(theta*dt) * b
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

        auto thetaJ = [this, sig, &VAVsolver](SmallVectorPtr const &xnew,
                                              SmallVectorPtr const &b) {
            auto y = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                      Copy, b->Values(), b->Stride(),
                                      b->MyLength(), b->NumVectors()));
            auto z = Teuchos::rcp(new Epetra_SerialDenseMatrix(
                                      Copy, b->Values(), b->Stride(),
                                      b->MyLength(), b->NumVectors()));
            CHECK_ZERO(z->Scale(sig));
            CHECK_ZERO(VAVsolver->SetVectors(*y, *z));
            CHECK_NONNEG(VAVsolver->Solve());
            return Teuchos::rcp(new Epetra_MultiVector(
                                    Copy, Epetra_LocalMap(y->M(), 0, Epetra_SerialComm()),
                                    y->A(), y->LDA(), y->N()));
        };

        return prolongate(*newton<SmallVector>(thetaF, thetaJ, xsmall, 1e-8));
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
    timestepper_->ams(num_exp_, num_init_exp_, sol1_, dt_, tmax_);
}

#endif
