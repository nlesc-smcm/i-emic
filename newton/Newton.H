//======================================================================
#ifndef NEWTON_H
#define NEWTON_H

//======================================================================
#include "NewtonDecl.H"
#include "GlobalDefinitions.H"

//======================================================================
// Constructor
template<typename Model, typename VectorPtr>
Newton<Model, VectorPtr>::Newton(Model model)
 	:
 	isInitialized_(false),
	isConverged_(false),
	backTracking_(true),  
	iter_(0),
	maxNumIterations_(10),
	toleranceRHS_(1.0e-3),
	normRHS_(1.0),
	numBackTrackingSteps_(10)
{
 	model_ = model;

	// Get control of the state in the model
	state_ = model_->GetState('V');
	Initialize();
}

//======================================================================
template<typename Model, typename VectorPtr>
void Newton<Model, VectorPtr>::Initialize()
{
	// This function needs some more stuff or it should be deleted
 	isInitialized_ = true;
}

//======================================================================
template<typename Model, typename VectorPtr>
void Newton<Model, VectorPtr>::Run()
{
	//
	isConverged_ = false;
	
	model_->ComputeRHS();
	normRHS_ = model_->GetNormRHS();
	for (iter_ = 0; iter_ != maxNumIterations_; ++iter_)
	{				
		//
		model_->ComputeJacobian();	
		model_->Solve();
		dir_ = model_->GetSolution('V');
		state_->Update(1.0, *dir_, 1.0);

		//
		model_->ComputeRHS();
		normRHStest_ = model_->GetNormRHS();

 		INFO("Newton:      iter: " << iter_ );
		INFO("Newton:      norm: " << normRHStest_ );

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
	{
		WARNING("Newton: ---> TROUBLE", __FILE__, __LINE__);
	}
	else
	{
		isConverged_ = true;
	}
}

//======================================================================
template<typename Model, typename VectorPtr>
void Newton<Model, VectorPtr>::RunBackTracking()
{
	// Initialize reduction with -1/2
	double reduction = -1.0 / 2;

	//
	for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_)
	{
		if (normRHStest_ < normRHS_)
		{
			DEBUG("Success...");
			break;
		}

		// Apply reduction to the state vector
		state_->Update(reduction, *dir_, 1.0);

		//
		model_->ComputeRHS();
		normRHStest_ = model_->GetNormRHS();

		INFO("Newton: --> backtracking: "
			 <<	" step: "      << backTrack_ 
			 << " reduction: " << reduction
			 << " norm: "      << normRHStest_);

		// Update reduction
		reduction /= 2.0;
	}
	if (backTrack_ == numBackTrackingSteps_)
		WARNING("Newton: --> BACKTRACKING FAILED", __FILE__, __LINE__);
}

#endif
