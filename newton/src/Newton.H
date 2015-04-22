#ifndef NEWTON_H
#define NEWTON_H

#include "NewtonDecl.H"
#include "GlobalDefinitions.H"

template<typename Model, typename Vector>
Newton<Model, Vector>::Newton(Model model)
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
	state_ = model->getState();
	initialize();
 	DEBUG("Leaving Newton constructor");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::initialize()
{
 	DEBUG("Entering Newton::initialize()");
 	isInitialized_ = true;
 	DEBUG("Leaving Newton::initialize()");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::run()
{
 	DEBUG("Entering Newton::run()");

	//
	model_->computeRHS();
	normRHS_ = model_->getNormRHS();
	for (iter_ = 0; iter_ != maxNumIterations_; ++iter_)
	{				
		//
		model_->computeJacobian();
		model_->solve();
		dir_ = model_->getSolution();
		state_->Update(1.0, *dir_, 1.0);

		//
		model_->computeRHS();
		normRHStest_ = model_->getNormRHS();
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
			runBackTracking();
		
		//
		normRHS_ = normRHStest_;
	}
	if (iter_ == maxNumIterations_)
		DEBUG("Newton: ---> TROUBLE");
	
 	DEBUG("Leaving Newton::run()");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::runBackTracking()
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
		model_->computeRHS();
		normRHStest_ = model_->getNormRHS();
	}
	if (backTrack_ == numBackTrackingSteps_)
		DEBUG("Newton: ---> TROUBLE");
			
	DEBUG("Leaving Newton::runBackTracking()");
}

#endif
