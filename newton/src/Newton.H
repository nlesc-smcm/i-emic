#ifndef NEWTON_H
#define NEWTON_H
#include "GlobalDefinitions.H"

//! This class finds a root of F(x) = 0. 
//!
//! Model can be a wrapper for THCM and other implicit models.
//! The user should implement the following methods:
//!     Model::evaluateF(Vector state_); // rhs evaluation
//!     Model::evaluateJ(Vector state_); // Jacobian evaluation
//!     Model::setState(Vector state_);  // set the state in the model
//!     Model::getState();               // get the state in the model
//!     Model::solve();                  // solve 
//!
//! Vector should be the type returned by Model::evaluateF().
//! It should have the following methods:
//!     Vector::update(double a, double b, Vector y) // perform  x = a*x+b*y
//!
//! Matrix should be the type returned by Model::evaluateJ().
//! It should have the following methods:


template<typename Model, typename Vector>
class Newton
{
	Model  model_;
	Vector state_;
	Vector dir_;
	
	//!
	int iter_;
	int maxNumIterations_;

	int backTrack_;
	int numBackTrackingSteps_;
	double toleranceRHS_;
	double normRHS_;
	double normRHStest_;

	//!
	bool isInitialized;
	bool backTracking;
	
public:
	Newton(Model model);
	void initialize();
	void run();
	void runBackTracking();
	
};

template<typename Model, typename Vector>
Newton<Model, Vector>::Newton(Model model)
 	:
 	isInitialized(false),
	backTracking(false),
	maxNumIterations_(1),
	toleranceRHS_(1.0e-6),
	normRHS_(1.0),
	numBackTrackingSteps_(10)
{
 	DEBUG("Entering Newton constructor");
 	model_ = model;
	state_ = model->getState();	
 	DEBUG("Leaving Newton constructor");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::initialize()
{
 	DEBUG("Entering Newton::initialize()");
 	isInitialized = true;
 	DEBUG("Leaving Newton::initialize()");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::run()
{
 	DEBUG("Entering Newton::run()");

	//
	model_->computeRHS();
	normRHS_ = model_->getNormRHS();
	double nrms;
	for (iter_ = 0; iter_ != maxNumIterations_; ++iter_)
	{
				
		//
		model_->computeJacobian();

		//
		model_->solve();

		//
		dir_ = model_->getSolution();

		//
 		nrms = model_->getNormState();
		DEBUG("Newton:      state (before): " << nrms );
		state_->Update(1.0, *dir_, 1.0);
		nrms = model_->getNormState();
		DEBUG("Newton:      state (after): " << nrms );

		//
		model_->computeRHS();
		normRHStest_ = model_->getNormRHS();
		DEBUG("Newton:      iter: " << iter_ );
		DEBUG("Newton:      norm: " << normRHS_ );
		DEBUG("Newton:  new norm: " << normRHStest_ );

		//
		if (backTracking)
			runBackTracking();

		//
		normRHS_ = normRHStest_;
	}		
 	DEBUG("Leaving Newton::run()");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::runBackTracking()
{
	DEBUG("Entering Newton::runBackTracking()");

	//
	double reduction = -1.0 / numBackTrackingSteps_;

	//
	for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_)
	{
		DEBUG("Newton --> backtracking:       iter: " << iter_ );
		DEBUG("Newton --> backtracking:  backtrack: " << backTrack_);
		DEBUG("Newton --> backtracking:       norm: " << normRHStest_);

		if (normRHStest_ < normRHS_)
		{
			DEBUG("Success...");
			break;
		}
		
		DEBVAR(backTrack_);
		//
		state_->Update(reduction, *dir_, 1.0);
		
		//
		model_->computeRHS();
		normRHStest_ = model_->getNormRHS();
	}
	if (backTrack_ == numBackTrackingSteps_)
		DEBUG("Newton: ---> TROUBLE");
			
	DEBUG("Leaving Newton::runBackTracking()");
}

#endif
