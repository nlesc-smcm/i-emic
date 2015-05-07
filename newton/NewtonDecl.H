#ifndef NEWTONDECL_H
#define NEWTONDECL_H

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

	int numPrecRecomputes_;
	
	int backTrack_;
	int numBackTrackingSteps_;
	double toleranceRHS_;
	double normRHS_;
	double normRHStest_;

	//!
	bool isInitialized_;
	bool isConverged_;
	bool backTracking_; //perhaps call this enableBacktracking_
	
public:
	Newton(Model model);
	void Initialize();
	void Run();
	void RunBackTracking();

	bool Converged(){ return isConverged_; }
	
	int  Iterations(){ return iter_; }
};

#endif
