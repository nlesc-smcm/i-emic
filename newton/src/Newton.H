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
	Vector dirs_;

	//! 
	int maxNumIterations_;
	double toleranceF_;
	double toleranceX_;

	//!
	bool isInitialized;
	
public:
	Newton(Model model);
	void initialize();
	void run();
	void step();
	
};


template<typename Model, typename Vector>
Newton<Model, Vector>::Newton(Model model)
 	:
 	isInitialized(false)
{
 	INFO("Entering Newton constructor");
 	model_ = model;
 	INFO("Leaving Newton constructor");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::initialize()
{
 	INFO("Entering Newton::initialize()");
 	//
 	isInitialized = true;
 	INFO("Leaving Newton::initialize()");
}
template<typename Model, typename Vector>
void Newton<Model, Vector>::run()
{
 	INFO("Entering Newton::run()");
 	INFO("Leaving Newton::run()");
}
template<typename Model, typename Vector>
void Newton<Model, Vector>::step()
{
 	INFO("Entering Newton::step()");
 	if (!isInitialized) initialize();
 	INFO("Leaving Newton::step()");
}
 
#endif
