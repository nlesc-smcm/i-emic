//==================================================================
#ifndef THETASTEPPER_H
#define THETASTEPPER_H
//==================================================================
#include "ThetaStepperDecl.H"
//==================================================================
template<typename Model, typename Vector>
ThetaStepper<Model, Vector>::ThetaStepper(Model model, Vector state)
	:
	theta_(1.0),
	time_(0.0),
	timeStep(1.0e-4),
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
	
	while (time_ < endTime_)
	{
		model_->Store();
		model_->SetTheta(theta_);
		model_->SetTimeStep(timeStep_);
		newtonSolver.Run();		
	}
}






#endif
