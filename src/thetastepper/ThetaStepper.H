//==================================================================
#ifndef THETASTEPPER_H
#define THETASTEPPER_H

//==================================================================
#include "ThetaStepperDecl.H"
#include "GlobalDefinitions.H"

//==================================================================
template<typename ThetaModel, typename ParameterList>
ThetaStepper<ThetaModel, ParameterList>::
ThetaStepper(ThetaModel model, ParameterList params)
	:
    model_    (model), 
    theta_    (params->get("theta", 1.0)),
    time_     (0.0),
    dt_       (params->get("initial time step size", 1.0e-03)),
    mindt_    (params->get("minimum step size", 1.0e-8)),
    maxdt_    (params->get("maximum step size", 1.0)),
    iscale_   (params->get("increase step size", 2.0)),
    dscale_   (params->get("decrease step size", 2.0)),
    inDays_   (params->get("timescale in days", 737.2685)), 
    inYears_  (inDays_ / 365.), 
    tend_     (params->get("end time (in y)", 10.0)),  
    step_     (0),
    nsteps_   (params->get("number of time steps", 10)),
    output_   (params->get("HDF5 output frequency", 1)),
    minK_     (params->get("minimum desired Newton iterations", 3)), 
    maxK_     (params->get("maximum desired Newton iterations", 3)), 
    Ntol_     (params->get("Newton tolerance", 1e-6)), 
    Niters_   (params->get("maximum Newton iterations", 8)),
    initWD_   (true)
{
    F_    = model_->getRHS('V');
    x_    = model_->getState('V');
    dx_   = model_->getSolution('V');

    dx_->PutScalar(0.0);

    // set theta in ThetaModel
    model_->setTheta(theta_);
}

//==================================================================
template<typename ThetaModel, typename ParameterList>
void ThetaStepper<ThetaModel, ParameterList>::run()
{
    
    bool test_step = ( nsteps_ < 0 ) ? true : step_ < nsteps_;
    while ( ( time_ < tend_ ) && ( test_step ) )
    {
        INFO("----------------------------------------------------------");
        INFO("Timestepping:    t = " << time_ << " y");

        model_->preProcess();

        // save current state
        model_->store();

        k_ = 0;
        for (; k_ != Niters_; ++k_)
        {
            model_->setTimestep(dt_);

            // compute time discretization F:
            // F(x) =  -(B d/dt x)/theta + F(x) + (theta-1)/theta * F(x_old)
            model_->computeRHS();

            // create jacobian of time discretization
            model_->computeJacobian();

            // solve for dx
            F_->Scale(-1.0);
            model_->solve(F_);
            normdx_ = Utils::normInf(dx_);

            // update state
            x_->Update(1.0, *dx_, 1.0);

            // compute new rhs
            model_->computeRHS();
            normF_ = Utils::norm(F_);

            INFO("  Newton solver ------------------------------------");
            INFO("                            iter     = " << k_);
            INFO("                           ||F||2    = " << normF_);
            INFO("                           ||dx||inf = " << normdx_);
            INFO("\n");

            if ( (normdx_ < Ntol_ ) && (normF_ < Ntol_) )
                break;
            else if ( normdx_ > 1e2)
            {
                WARNING("Norm exploding! ||dx||inf = "
                        << normdx_ << "\n", __FILE__, __LINE__);
                k_ = Niters_;
                break;
            }
        }

        if (k_ == Niters_)
        {
            WARNING("Newton did not converge! ||F|| = "
                    << normF_ << "\nRestoring model",
                    __FILE__, __LINE__);

            Utils::save(F_, "failed_rhs"); // Print failed residual
            
            INFO("    adjusting time step.. old dt = " << dt_);
            dt_ = std::max(dt_ / dscale_, mindt_);
            INFO("    adjusting time step.. new dt = " << dt_);

            if (dt_ == mindt_)
            {
                INFO("\n ================================ \n  " <<
                     " minimum timestep reached, exiting... " <<
                     "\n ================================ \n");
                return;
            }
            model_->restore();
            continue;
        }

        step_++;
        time_ += dt_ * inYears_;

        INFO("  Newton converged, next time step -----------------");
        INFO("                step = " << step_);
        INFO("                time = " << time_);
        INFO("           ||F||2    = " << normF_);
        INFO("           ||dx||inf = " << normdx_);
        INFO("\n");

        model_->postProcess();

        // additional save
        if ((output_ > 0) && ((step_ % output_) == 0))
        {
            std::stringstream outFile;
            outFile << "transient_" << std::setprecision(8)
                    << time_ << ".h5";

            model_->saveStateToFile(outFile.str());
        }

        writeData();

        // Timestep adjustments
        if (k_ < minK_)
            dt_ = std::min(dt_ * iscale_, maxdt_);
        else if (k_ > maxK_)
            dt_ = std::max(dt_ / dscale_, mindt_);

        test_step = ( nsteps_ < 0 ) ? true : step_ < nsteps_;
    }
}

//==================================================================
template<typename ThetaModel, typename ParameterList>
void ThetaStepper<ThetaModel, ParameterList>::
writeData()
{
    // Data string
    std::ostringstream tdatastring;
    
    if (initWD_) // write description of entries
    {
        tdatastring << std::setw(_FIELDWIDTH_)
                    << "# time_(y)"
                    << std::setw(_FIELDWIDTH_/2.)
                    << "step"
                    << std::setw(_FIELDWIDTH_)
                    << "dt_(y)"
                    << std::setw(_FIELDWIDTH_)
                    << "|x|"
                    << std::setw(_FIELDWIDTH_/2.)
                    << "NR"
                    << model_->writeData(true);

        WRITETDATA(tdatastring.str());
        initWD_ = false;
    }
    tdatastring.str("");
    tdatastring.clear();
    tdatastring.precision(_PRECISION_);

    tdatastring << std::scientific
                << std::setw(_FIELDWIDTH_) <<  time_
                << std::setw(_FIELDWIDTH_/2) <<  step_
                << std::setw(_FIELDWIDTH_) <<  dt_ * inYears_
                << std::setw(_FIELDWIDTH_)
                << Utils::norm(model_->getState('V'))
                << std::setw(_FIELDWIDTH_/2) <<  k_
                << model_->writeData();

    WRITETDATA(tdatastring.str());
}

#endif
