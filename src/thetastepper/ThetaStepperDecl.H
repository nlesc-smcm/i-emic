#ifndef THETASTEPPERDECL_H
#define THETASTEPPERDECL_H

//! This class performs a time integration using the theta-method.

//! ThetaModel should be an instantiation of the class template
//! Theta<Model> (see Theta.H)

template<typename ThetaModel, typename ParameterList>
class ThetaStepper
{

    using VectorPtr = typename ThetaModel::element_type::VectorPtr;
    using Vector    = typename VectorPtr::element_type;

	ThetaModel model_;

    //! view of theta method RHS
	VectorPtr F_;
    //! view of transient state
    VectorPtr x_;
    //! view of solution from solve with transient problem
    VectorPtr dx_;
    //! our copy of the old transient state
    VectorPtr xold_;

    //! theta
	double theta_;
    //! time
	double time_;
    //! initial time step size
	double dt_;
    //! minimum step size
	double mindt_;
    //! maximum step size
	double maxdt_;
    //! increase step size
 	double iscale_;
    //! decrease step size
	double dscale_;

    //! nondimensional time to days conversion
	double inDays_;
    //! nondimensional time to years conversion
	double inYears_;

    //! end time (in y)
	double tend_;
    //! time step
    int step_;
    //! number of time steps
    int nsteps_;
    //! HDF5 output frequency
    int output_;

    //! Newton step
    int k_;
    //! minimum desired Newton iterations
	int minK_;
    //! maximum desired Newton iterations
	int maxK_;
    //! Newton tolerance
    double Ntol_;
    //! maximum Newton iterations
    int Niters_;
    
    //! norm rhs
    double normF_;
    //! infnorm dx
    double normdx_;

    //! output initialization flag
    bool initWD_;
    

public:
	ThetaStepper(ThetaModel model, ParameterList params);
	void run();
    
private:    
	void writeData();
    
};

#endif
