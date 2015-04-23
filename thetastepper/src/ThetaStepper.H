//==================================================================
#ifndef THETASTEPPER_H
#define THETASTEPPER_H
//==================================================================
#include "ThetaStepperDecl.H"
#include "Newton.H"
//==================================================================
template<typename Model, typename Vector>
ThetaStepper<Model, Vector>::ThetaStepper(Model model, Vector state)
	:
	theta_(1.0),
	time_(0.0),
	timestep_(1.0e-4),
	startTime_(0.0),
	endTime_(50.0)
{
	model_ = model;
	state_ = state;
}
//==================================================================
template<typename Model, typename Vector>
void ThetaStepper<Model, Vector>::Run()
{
	Newton<Model, Vector> newtonSolver(model_, state_);
	int iters;
	while (time_ < endTime_)
	{
		model_->Store();
		model_->SetTheta(theta_);
		model_->SetTimestep(timestep_);
		newtonSolver.Run();

		/*
		if (!newtonSolver.Converged())
		{		
			model->Restore();
			continue;
		}
		*/

		time_ += timestep_;
		DEBUG("ThetaStepper<Model, Vector>::Run()");
		DEBUG("  time = " << time_);
		
		// iters = newtonSolver.Iterations();
		// AdjustTimestep(iters); // not implemented yet
		// AdjustTheta(iters); // not implemented yet
	}
}






#endif
