#ifndef COUPLEDTHETAMODEL_H
#define COUPLEDTHETAMODEL_H

#include "GlobalDefinitions.H"
#include "ThetaModel.H"
#include "CoupledModel.H"

class CoupledThetaModel: public CoupledModel
{
public:
    using ConstVectorPtr = typename CoupledModel::ConstVectorPtr;

    template<typename ParameterList>
    CoupledThetaModel(CoupledModel const &model, ParameterList params)
        :
        CoupledModel(model)
        {}

    virtual ~CoupledThetaModel() {}

    //------------------------------------------------------------------
    void initStep(double dt)
        {
            int num = models_.size();
            for (int i = 0; i != num; ++i)
            {
                auto &model = models_[i];
                if (i == OCEAN)
                    std::static_pointer_cast<ThetaModel<Ocean> >(model)->initStep(dt);
                else if (i == ATMOS)
                    std::static_pointer_cast<ThetaModel<Atmosphere> >(model)->initStep(dt);
                else if (i == SEAICE)
                    std::static_pointer_cast<ThetaModel<SeaIce> >(model)->initStep(dt);
            }
        }

    //------------------------------------------------------------------
    void setState(std::shared_ptr<const Combined_MultiVec> state)
        {
            int num = models_.size();
            for (int i = 0; i != num; ++i)
            {
                auto &model = models_[i];
                if (i == OCEAN)
                    std::static_pointer_cast<ThetaModel<Ocean> >(model)->setState(
                        Teuchos::rcp_dynamic_cast<const Epetra_Vector>((*state)(i)));
                else if (i == ATMOS)
                    std::static_pointer_cast<ThetaModel<Atmosphere> >(model)->setState(
                        Teuchos::rcp_dynamic_cast<const Epetra_Vector>((*state)(i)));
                else if (i == SEAICE)
                    std::static_pointer_cast<ThetaModel<SeaIce> >(model)->setState(
                        Teuchos::rcp_dynamic_cast<const Epetra_Vector>((*state)(i)));
            }
        }
};

#endif
