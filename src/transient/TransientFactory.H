#ifndef TRANSIENTFACTORY_H
#define TRANSIENTFACTORY_H

//======================================================================
#include "GlobalDefinitions.H"
#include "Utils.H"

#include "Newton.H"
#include "Transient.hpp"
#include "AdaptiveTransient.H"
#include "StochasticThetaModel.H"
#include "StochasticProjectedThetaModel.H"
#include "ScoreFunctions.H"

#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"

//! This file implements factory functions that construct a Transient.
//! This is mainly to obfuscate the complicated dependency hierarchy that
//! is described here.
//!
//! Model is a model such as Ocean.
//!
//! ThetaModel is a wrapper which overrides the RHS and Jacobian such
//! that they are the RHS and Jacobian of the Theta iteration.
//!
//! ((Stochastic)?(Projected)?)ThetaModel are the same as the previous
//! but then their stochastic and projected (time stepping) variants.
//!
//! Transient is a class that implements all time steppers and rare
//! event methods, which expects a time step. Note that a ThetaModel
//! itself does not implement a time step.
//!
//! Newton implements the Newton method which can be used in an implicit
//! time step.
//!
//! get_time_step uses Newton and a ThetaModel to apply one time step.
//! This can then be passed to Transient to implement a timestepper.
//!
//! AdaptiveTransient is a specialization of Transient which also allows
//! for adaptive time steps based on the convergence of Newton. This is
//! useful when computing a bifurcation diagram using time stepping.
//!
//! Using this terminology, we have a ThetaModel<Model> which is passed
//! to get_time_step, where it is combined with Newton to construct a
//! time step. This time step is then passed to Transient to make an
//! actual transient or rare event method.

// Most methods in this file can be used as
// Teuchos::RCP<Transient<Teuchos::RCP<const Epetra_Vector> > >

//! Obtain a time step from a time step model, for instance a ThetaModel.
//! This should only be used internally in the TransientFactory methods.
template<typename Model, typename ParameterList>
auto get_time_step(Model const &model, ParameterList pars)
{
    Teuchos::RCP<Newton<Model> > newton = Teuchos::rcp(new Newton<Model>(model, pars));
    // Function to perform one stochastic time step
    return [newton, model](Teuchos::RCP<const Epetra_Vector> const &x, double dt) {
        TIMER_SCOPE("TransientFactory: Time step");

        model->setState(x);
        model->initStep(dt);
        return newton->run(x);
    };
}

//! Factory function for a standard time stepper, which starts from the current
//! state of the model.
template<typename Model, typename ParameterList>
auto TransientFactory(
    Model model, ParameterList pars)
{
    auto theta_model = Teuchos::rcp(new ThetaModel<typename Model::element_type>(
                                        *model, pars));
    auto timestepper = Teuchos::rcp(
        new AdaptiveTransient<decltype(theta_model)>(theta_model, pars));

    timestepper->set_parameters(*pars);
    return timestepper;
}

//! Factory function for a standard time stepper, which starts from the initial
//! guess x0.
template<typename Model, typename ParameterList>
auto TransientFactory(
    Model model, ParameterList pars,
    Teuchos::RCP<const Epetra_Vector> x0)
{
    auto theta_model = Teuchos::rcp(new ThetaModel<typename Model::element_type>(
                                        *model, pars));
    auto timestepper = Teuchos::rcp(
        new AdaptiveTransient<decltype(theta_model)>(theta_model, pars, x0));

    timestepper->set_parameters(*pars);
    return timestepper;
}

//! Factory function for a rare event method, which computes transitions between
//! sol1 and sol2, with a possibly Teuchos::null unstable steady state sol3 in
//! between. V is the space which can be used for a projected time step. This
//! should be Teuchos::null in case not projected time step is desired.
template<typename Model, typename ParameterList>
auto TransientFactory(
    Model model, ParameterList pars,
    Teuchos::RCP<const Epetra_Vector> sol1,
    Teuchos::RCP<const Epetra_Vector> sol2,
    Teuchos::RCP<const Epetra_Vector> sol3,
    Teuchos::RCP<const Epetra_MultiVector> V)
{
    std::function<double(Teuchos::RCP<const Epetra_Vector> const &)> score_fun;
    Teuchos::RCP<ThetaModel<typename Model::element_type> > theta_model;

    if (V != Teuchos::null)
    {
        Teuchos::RCP<StochasticProjectedThetaModel<typename Model::element_type> >
            projected_theta_model = Teuchos::rcp(
                new StochasticProjectedThetaModel<typename Model::element_type>(
                    *model, pars, V));

        sol1 = projected_theta_model->restrict(*sol1);
        sol2 = projected_theta_model->restrict(*sol2);
        sol3 = projected_theta_model->restrict(*sol3);

        // if (std::is_same<Model, Ocean>::value)
        if (pars->get("dof", 1) == 6)
            score_fun = get_projected_ocean_score_function(sol1, sol2, sol3, V);
        else
            score_fun = get_projected_default_score_function(sol1, sol2, sol3, V);

        theta_model = projected_theta_model;
    }
    else
    {
        // if (std::is_same<Model, Ocean>::value)
        if (pars->get("dof", 1) == 6)
            score_fun = get_ocean_score_function(sol1, sol2, sol3);
        else
            score_fun = get_default_score_function(sol1, sol2, sol3);

        theta_model = Teuchos::rcp(new StochasticThetaModel<typename Model::element_type>(
                                       *model, pars));
    }

    auto time_step = get_time_step(theta_model, pars);
    auto timestepper = Teuchos::rcp(
        new Transient<Teuchos::RCP<const Epetra_Vector> >(
            time_step, score_fun, sol1, sol1->GlobalLength()));

    timestepper->set_parameters(*pars);

    unsigned int seed = pars->get("ams seed", 0);
    if (seed == 0)
    {
        static thread_local std::random_device rd;
        seed = rd();
    }

    int *seed_ptr = reinterpret_cast<int *>(&seed);
    CHECK_ZERO(sol1->Map().Comm().Broadcast(seed_ptr, 1, 0));

    StochasticBase::write_seed(sol1->Map().Comm(), seed, "Global seed");
    timestepper->set_random_engine(seed);
    return timestepper;
}

//! Wrapper for the previous factory function where the space V is loaded
//! from an mtx file if the "space" parameter is not empty.
template<typename Model, typename ParameterList>
auto TransientFactory(
    Model model, ParameterList pars,
    Teuchos::RCP<const Epetra_Vector> sol1,
    Teuchos::RCP<const Epetra_Vector> sol2,
    Teuchos::RCP<const Epetra_Vector> sol3)
{
    Teuchos::RCP<Epetra_MultiVector> V = Teuchos::null;
    std::string space = pars->get("space", "");
    if (space != "")
    {
        // Use a map that is constructed with the most basic
        // constructor for loading V since
        // MatrixMarketFileToMultiVector does not allow for anything
        // else...
        Epetra_BlockMap const &solveMap = model->getState('V')->Map();
        Epetra_Map map(solveMap.NumGlobalElements(), 0, *model->Comm());
        Epetra_MultiVector* Vptr;
        CHECK_ZERO(EpetraExt::MatrixMarketFileToMultiVector(
                       space.c_str(), map, Vptr));

        Epetra_Import import(solveMap, map);
        V = Teuchos::rcp(new Epetra_MultiVector(solveMap, Vptr->NumVectors()));
        CHECK_ZERO(V->Import(*Vptr, import, Insert));

        delete Vptr;
    }
    return TransientFactory(model, pars, sol1, sol2, sol3, V);
}

#endif
