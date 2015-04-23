//==================================================================
#ifndef THETASTEPPER_H
#define THETASTEPPER_H

//==================================================================
#include "ThetaStepperDecl.H"
#include "Newton.H"
#include "GlobalDefinitions.H"

//==================================================================
template<typename Model, typename Vector>
ThetaStepper<Model, Vector>::ThetaStepper(Model model, Vector state)
	:
	theta_(1.0),
	time_(0.0),
	timestep_(1.0e-4),
	minTimestep_(1.0e-8),
	maxTimestep_(1.0e-1),
	startTime_(0.0),
	endTime_(50.0),
	timestepScale1_(1.20),
	timestepScale2_(2.00),
	timestepScale3_(1.05),
	minAllowedIterations_(4),
	maxAllowedIterations_(6)	
{
	model_ = model;
	state_ = state;
}

//==================================================================
template<typename Model, typename Vector>
void ThetaStepper<Model, Vector>::Run()
{
	Newton<Model, Vector> newtonSolver(model_, state_);
	int  iters;
	bool convergence;
	
	while (time_ < endTime_)
	{
		model_->DumpState();
		model_->Store();
		model_->SetTheta(theta_);
		model_->SetTimestep(timestep_);

		newtonSolver.Run();
		convergence = newtonSolver.Converged();
		
		if (!convergence)
		{		
			model_->Restore();
			continue;
		}
		
		time_ += timestep_;

		DEBUG("ThetaStepper<Model, Vector>::Run()");
		DEBUG("  time     = " << time_);
		DEBUG("  timestep = " << timestep_);
		
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
	if (iters <= minAllowedIterations_)
 	{
		DEBUG("  AdjustTimestep: increasing with a factor "
			  << timestepScale1_);
		timestep_ = MIN(timestep_ * timestepScale1_, maxTimestep_);
	}
	else if (iters >= maxAllowedIterations_ and convergence)
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
	
#endif
