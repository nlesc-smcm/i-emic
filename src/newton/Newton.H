#ifndef NEWTON_H
#define NEWTON_H

template<typename Model, typename State,
		 typename Matrix, typename LinearSolver>
class Newton
{
	//! Model can be a wrapper for THCM and other implicit models.
	//! The user should implement the following methods:
	//!     Model::evaluateF(Statetype state_); // rhs evaluation
	//!     Model::evaluateJ(Statetype state_); // Jacobian evaluation
	//!     Model::setState(Statetype state_);  // set the state in the model
	//!     Model::getState();                  // get the state in the model

	Model model_;
	State state_;
	
	LinearSolver linearSolver_;

	//! 
	size_t maxNumIterations_;
	double toleranceF_;
	double toleranceX_;
	
public:
	// constructor
	Newton(Model model, char *parameterFile);
	void step();
	
private:


	
}
#endif
