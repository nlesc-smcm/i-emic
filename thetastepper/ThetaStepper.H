//==================================================================
#ifndef THETASTEPPER_H
#define THETASTEPPER_H

//==================================================================
#include "ThetaStepperDecl.H"
#include "Newton.H"
#include "GlobalDefinitions.H"

//==================================================================
template<typename Model, typename Vector>
ThetaStepper<Model, Vector>::ThetaStepper(Model model)
	:
	theta_(1.0),
	time_(0.0),
	timestep_(1.0e-3),
	minTimestep_(1.0e-8),
	maxTimestep_(1.0),
	startTime_(0.0),
	endTime_(500.0),
	timestepScale1_(1.20),
	timestepScale2_(2.00),
	timestepScale3_(1.05),
	inYears_(2.0),
	inDays_(2.0 * 365),
	minAllowedIterations_(4),
	maxAllowedIterations_(5)	
{
	model_ = model;
	state_ = model_->GetState('V');
}

//==================================================================
template<typename Model, typename Vector>
void ThetaStepper<Model, Vector>::Run()
{
	Newton<Model, Vector> newtonSolver(model_);
	int  iters;
	bool convergence;

	Info();
	
	while (time_ < endTime_)
	{
		model_->DumpState();
		model_->Store();
		model_->SetTheta(theta_);
		model_->SetTimestep(timestep_);

		INFO("  time     = " << time_ * inYears_ << " Years");
		INFO("  timestep = " << timestep_ * inDays_ << " Days");
		
		newtonSolver.Run();
		convergence = newtonSolver.Converged();
		
		if (!convergence)
		{
			INFO("No convergence in Newton solve");
			INFO("   Restoring the model and decreasing the timestep");
			model_->Restore();
			AdjustTimestep(iters, convergence); 
			continue;
		}
		
		time_ += timestep_;
		
		iters = newtonSolver.Iterations();

		AdjustTimestep(iters, convergence); 

		// AdjustTheta(iters); // not implemented yet
	}
}

//==================================================================
template<typename Model, typename Vector>
void ThetaStepper<Model, Vector>::AdjustTimestep(int iters,
												 bool convergence)
{
	if (iters < minAllowedIterations_ and convergence)
 	{
		DEBUG("  AdjustTimestep: increasing with a factor "
			  << timestepScale1_);
		timestep_ = MIN(timestep_ * timestepScale1_, maxTimestep_);
	}
	else if (iters > maxAllowedIterations_ and convergence)
	{
		DEBUG("  AdjustTimestep: decreasing with a factor "
			  << timestepScale1_);
		timestep_ = MAX(timestep_ / timestepScale1_, minTimestep_);
	}
	else if (!convergence)
	{
		DEBUG("  AdjustTimestep: RESET ");
		DEBUG("  AdjustTimestep: decreasing with a factor "
			  << timestepScale2_);
		timestep_ = MAX(timestep_ / timestepScale2_, minTimestep_);		
	}	
}

//==================================================================
template<typename Model, typename Vector>
void ThetaStepper<Model, Vector>::Info()
{
	INFO("--------------------------------------------------------");
	INFO("Timestepping run with... ");
	INFO("                  theta: " << theta_);
	INFO("               timestep: " << timestep_);
	INFO("    timestep boundaries: "
		 << "[" << minTimestep_ << ", " << maxTimestep_ << "]");
	INFO("        time boundaries: " 
		 << "[" << startTime_ << ", " << endTime_ << "]");
	INFO("      newton boundaries: " 
		 << "[" << minAllowedIterations_
		 << ", " << maxAllowedIterations_ << "]");
	INFO("      timestep scalings... "
		 << "increase/decrease: " << timestepScale1_
		 << " big decrease: " << timestepScale2_);
	INFO("--------------------------------------------------------");
}

#endif
