#ifndef NEWTONDECL_H
#define NEWTONDECL_H

//! This class finds a root of F(x) = 0. 
//!
//! It should have the following methods:

//! VectorPtr model_->GetState();
//! VectorPtr model_->GetSolution();
//! double    model_->GetNormRHS();

//! model_->ComputeRHS();
//! model_->ComputeJacobian();
//! model_->Solve();

//! state_->Update(1.0, *dir_, 1.0); --> in Epetra_MultiVector

//! In my opinion, a model should have data members:
//!   state
//!   rhs
//!   jacobian


template<typename Model, typename VectorPtr>
class Newton
{
	Model  model_;
	VectorPtr state_;
	VectorPtr dir_;

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
