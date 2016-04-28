#ifndef THETASTEPPERDECL_H
#define THETASTEPPERDECL_H

//! This class performs a timestepping using the theta-method.
//! The templated model is expected to be in the following form:
//! ...
//! The templated model should have the following methods:
//! ...
//! This class depends on the Newton class in this project.

template<typename Model, typename VectorPtr>
class ThetaStepper
{
	Model  model_;
	VectorPtr state_;

	// timestepping specific datamembers
	double theta_;
	double time_;
	double timestep_;
	double minTimestep_;
	double maxTimestep_;

	double startTime_;
	double endTime_;

 	double timestepScale1_;
	double timestepScale2_;
	double timestepScale3_;

	double inDays_;
	double inYears_;

	int minAllowedIterations_;
	int maxAllowedIterations_;

	

public:

	ThetaStepper(Model model);
	void Run();
private:
	void AdjustTimestep(int iters, bool convergence);
	void Info();
};

#endif
