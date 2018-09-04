#ifndef THETA_H
#define THETA_H

//! Model is assumed to be of RCP type

template<typename Model>
class Theta
{
public:
    using VectorPtr  = typename Model::element_type::VectorPtr;
    using Vector     = typename VectorPtr::element_type;
    using MatrixPtr  = typename Model::element_type::MatrixPtr;
    using Matrix     = typename MatrixPtr::element_type;

private:

    //! rcp to model object
    Model model_;
    
    //! theta method parameter
    double theta_;

    //! time step size
    double timestep_;

    //! view of the model state
    VectorPtr stateView_;

    //! view of the model rhs
    VectorPtr rhsView_;

    //! state to store and restore
    VectorPtr oldState_;

    //! rhs to store and restore
    VectorPtr oldRHS_;
    
    //! state time derivative
    VectorPtr stateDot_;
    
private:
    
    //! constructor
    Theta(Model model)
        :
        model_     (model),
        theta_     (1.0),
        timestep_  (1.0e-3),
        stateView_ (model->getState('V')),
        rhsView_   (model->getRHS('V')),
        oldState_  (model->getState('C')),
        oldRHS_    (model->getRHS('C'))
        stateDot_  (model->getState('C')),
        {}

    //! compute theta method right hand side
    void computeRHS()
        {
            // compute right hand side of discretization
            model_->computeRHS();

            ##//...tododododo
        }
    



        
        
};







#endif
