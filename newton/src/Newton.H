//======================================================================
#ifndef NEWTON_H
#define NEWTON_H
//======================================================================
#include "NewtonDecl.H"
#include "GlobalDefinitions.H"
//======================================================================
template<typename Model, typename Vector>
Newton<Model, Vector>::Newton(Model model, Vector state)
 	:
 	isInitialized_(false),
	backTracking_(true),
	maxNumIterations_(10),
	toleranceRHS_(1.0e-6),
	normRHS_(1.0),
	numBackTrackingSteps_(10)
{
 	DEBUG("Entering Newton constructor");
 	model_ = model;
	state_ = state;
	Initialize();
 	DEBUG("Leaving Newton constructor");
}
//======================================================================
template<typename Model, typename Vector>
void Newton<Model, Vector>::Initialize()
{
 	DEBUG("Entering Newton::initialize()");
 	isInitialized_ = true;
 	DEBUG("Leaving Newton::initialize()");
}
//======================================================================
template<typename Model, typename Vector>
void Newton<Model, Vector>::Run()
{
 	DEBUG("Entering Newton::run()");

	//
	model_->ComputeRHS();
	normRHS_ = model_->GetNormRHS();
	for (iter_ = 0; iter_ != maxNumIterations_; ++iter_)
	{				
		//
		model_->ComputeJacobian();
		model_->Solve();
		dir_ = model_->GetSolution();
		state_->Update(1.0, *dir_, 1.0);

		//
		model_->ComputeRHS();
		normRHStest_ = model_->GetNormRHS();
		DEBUG("Newton:      iter: " << iter_ );
		DEBUG("Newton:      norm: " << normRHS_ );
		DEBUG("Newton:  new norm: " << normRHStest_ );

		//
		if (normRHStest_ < toleranceRHS_)
		{			
			DEBUG("Success...");
			break;
		}
		
		if (backTracking_ and (normRHS_ < normRHStest_) )
			RunBackTracking();
		
		//
		normRHS_ = normRHStest_;
	}
	if (iter_ == maxNumIterations_)
		DEBUG("Newton: ---> TROUBLE");
	
 	DEBUG("Leaving Newton::run()");
}
//======================================================================
template<typename Model, typename Vector>
void Newton<Model, Vector>::RunBackTracking()
{
	DEBUG("Entering Newton::runBackTracking()");

	//
	double reduction = -1.0 / 2;

	//
	for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_)
	{
		if (normRHStest_ < normRHS_)
		{
			DEBUG("Success...");
			break;
		}
		DEBUG("Newton --> backtracking:       iter: " << iter_ );
		DEBUG("Newton --> backtracking:  backtrack: " << backTrack_);
		DEBUG("Newton --> backtracking:  reduction: " << reduction);
		DEBUG("Newton --> backtracking:       norm: " << normRHStest_);

		//
		state_->Update(reduction, *dir_, 1.0);
		reduction /= 2.0;
		
		//
		model_->ComputeRHS();
		normRHStest_ = model_->GetNormRHS();
	}
	if (backTrack_ == numBackTrackingSteps_)
		DEBUG("Newton: ---> TROUBLE");
			
	DEBUG("Leaving Newton::runBackTracking()");
}

#endif
