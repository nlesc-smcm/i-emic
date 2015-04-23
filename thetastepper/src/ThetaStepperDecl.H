#ifndef THETASTEPPERDECL_H
#define THETASTEPPERDECL_H

//! This class performs a timestepping using the theta-method.
//! The templated model is expected to be in the following form:
//! ...
//! The templated model should have the following methods:
//! ...
//! This class depends on the Newton class in this project.

template<typename Model, typename Vector>
class ThetaStepper
{
	Model  model_;
	Vector state_;

	// timestepping specific member
	double theta_;
	double time_;
	double timestep_;
	double startTime_;
	double endTime_;

public:

	ThetaStepper(Model model, Vector state);
	void Run();

};

#endif
