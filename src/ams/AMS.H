#ifndef AMS_H
#define AMS_H

//======================================================================
#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Transient.hpp"
#include "StochasticThetaModel.H"
#include "StochasticProjectedThetaModel.H"
#include "ScoreFunctions.H"

#include "Epetra_MultiVector.h"
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
    Teuchos::RCP<ThetaModel<typename Model::element_type> > thetaModel_;

    int dof_;
    int var_;

    unsigned int ams_seed_;

    VectorPtr sol1_;
    VectorPtr sol2_;
    VectorPtr sol3_;

    Teuchos::RCP<Transient<VectorPtr> > timestepper_;

    Teuchos::RCP<Epetra_MultiVector> V_;

public:
    template<typename ParameterList>
    AMS(Model model, ParameterList pars,
        VectorPtr sol1, VectorPtr sol2,
        VectorPtr sol3 = Teuchos::null);
    ~AMS();

    template<typename ParameterList>
    void initialize(ParameterList pars);
    void set_space(Teuchos::RCP<Epetra_MultiVector> const &V);

    auto get_time_step();

    void run();

    double get_probability();
    double get_mfpt();

protected:

    unsigned int get_ams_seed();
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
    VectorPtr sol1, VectorPtr sol2,
    VectorPtr sol3)
    :
    model_(model),
    dof_(pars->get("dof", 6)),
    var_(pars->get("var", 1)),
    ams_seed_(pars->get("ams seed", 0)),
    sol1_(sol1),
    sol2_(sol2),
    sol3_(sol3),
    V_(Teuchos::null)
{
    initialize(pars);
}

template<typename Model>
template<typename ParameterList>
void AMS<Model>::initialize(ParameterList pars)
{
    std::string space = pars->get("space", "");
    if (space != "")
    {
        // Use a map that is constructed with the most basic
        // constructor for loading V since
        // MatrixMarketFileToMultiVector does not allow for anything
        // else...
        Epetra_BlockMap const &solveMap = model_->getState('V')->Map();
        Epetra_Map map(solveMap.NumGlobalElements(), 0, *model_->Comm());
        Epetra_MultiVector* Vptr;
        CHECK_ZERO(EpetraExt::MatrixMarketFileToMultiVector(
                       space.c_str(), map, Vptr));

        Epetra_Import import(solveMap, map);
        V_ = Teuchos::rcp(new Epetra_MultiVector(solveMap, Vptr->NumVectors()));
        CHECK_ZERO(V_->Import(*Vptr, import, Insert));

        delete Vptr;
    }

    TimeStep time_step;
    DistFun dist_fun;

    if (V_ != Teuchos::null)
    {
        time_step = get_time_step();

        Teuchos::RCP<StochasticProjectedThetaModel<typename Model::element_type> >
            projectedThetaModel = Teuchos::rcp(
                new StochasticProjectedThetaModel<typename Model::element_type>(
                    *model_, pars, V_));

        sol1_ = projectedThetaModel->restrict(*sol1_);
        sol2_ = projectedThetaModel->restrict(*sol2_);
        sol3_ = projectedThetaModel->restrict(*sol3_);

        // if (std::is_same<Model, Ocean>::value)
        if (dof_ == 6)
            dist_fun = get_projected_ocean_score_function(sol1_, sol2_, sol3_, V_);
        else
            dist_fun = get_projected_default_score_function(sol1_, sol2_, sol3_, V_);

        thetaModel_ = projectedThetaModel;
    }
    else
    {
        time_step = get_time_step();

        // if (std::is_same<Model, Ocean>::value)
        if (dof_ == 6)
            dist_fun = get_ocean_score_function(sol1_, sol2_, sol3_);
        else
            dist_fun = get_default_score_function(sol1_, sol2_, sol3_);

        thetaModel_ = Teuchos::rcp(new StochasticThetaModel<typename Model::element_type>(
                                       *model_, pars));
    }

    timestepper_ = Teuchos::rcp(new Transient<VectorPtr>(
                                    time_step, dist_fun, sol1_->GlobalLength()));

    timestepper_->set_parameters(*pars);
    timestepper_->set_random_engine(get_ams_seed());
}

template<typename Model>
void AMS<Model>::set_space(Teuchos::RCP<Epetra_MultiVector> const &V)
{
    V_ = V;
}

template<typename Model>
auto AMS<Model>::get_time_step()
{
    // Function to perform one stochastic time step
    return [this](VectorPtr const &x, double dt) {
        TIMER_SCOPE("AMS: Time step");

        thetaModel_->setState(x);
        thetaModel_->initStep(dt);

        // Deterministic theta stepper:
        // M * u_n + dt * theta * F(u_(n+1)) + dt * (1-theta) * F(u_n) - M * u_(n+1) = 0
        // Noise is added in an explicit manner
        auto thetaF = [this](VectorPtr const &xnew) {
            TIMER_SCOPE("AMS: F");

            thetaModel_->setState(xnew);
            thetaModel_->computeRHS();
            return thetaModel_->getRHS('V');
        };

        // J2 = theta * dt * J - M, so J2*x = b
        // We write this as
        // J2 = J - 1/(theta*dt) * M, J2 * x = 1/(theta*dt) * b
        auto thetaJ = [this](VectorPtr const &xnew, VectorPtr const &b) {
            TIMER_SCOPE("AMS: Jacobian solve");
            thetaModel_->setState(xnew);
            thetaModel_->computeJacobian();
            thetaModel_->solve(b);
            return thetaModel_->getSolution('V');
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
    timestepper_->run(sol1_);
}

template<typename Model>
double AMS<Model>::get_probability()
{
    return timestepper_->get_probability();
}

template<typename Model>
double AMS<Model>::get_mfpt()
{
    return timestepper_->get_mfpt();
}

template<typename Model>
unsigned int AMS<Model>::get_ams_seed()
{
    unsigned int seed = ams_seed_;
    if (seed == 0)
    {
        static thread_local std::random_device rd;
        seed = rd();
    }

    int *seed_ptr = reinterpret_cast<int *>(&seed);
    CHECK_ZERO(sol1_->Map().Comm().Broadcast(seed_ptr, 1, 0));

    write_seed(sol1_->Map().Comm(), seed, "Global seed");

    return seed;
}

#endif
