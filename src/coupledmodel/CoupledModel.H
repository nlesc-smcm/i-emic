#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

//! models
#include "Model.H"
#include "Ocean.H"
#include "AtmospherePar.H"
#include "SeaIce.H"

//! vector and matrix helpers
#include "Combined_MultiVec.H"
#include "CouplingBlock.H"

#include <vector>
#include <memory>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "BelosConfigDefs.hpp"
#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"
#include "BelosTypes.hpp"
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>

/*------------------------------------------------------------------
//! This class assembles multiple sub-models into a CoupledModel.
//! Which models and couplings are used should be supplied in the
//! accompanying xml file (coupledmodel_params.xml).

//! Most routines in CoupledModel simply call their corresponding
//! routines in the sub-models. For instance, a call to
//! computeJacobian() translates to a call to computeJacobian() in the
//! sub-models.

//! The frequency at which the states of the sub-models are
//! synchronized depends on the solving scheme.
------------------------------------------------------------------*/

//! Forward declarations
class Epetra_Comm;
class Epetra_IntVector;
class Epetra_Vector;

template<typename ModelPtr>
class BelosOp;

class CoupledModel
{
public:
    using VectorPtr  = std::shared_ptr<Combined_MultiVec>;
    using LandMask   = Utils::MaskStruct;

private:

    //! submodel identifiers
    int OCEAN, ATMOS, SEAICE;

    //! Trilinos-based parallel implicit ocean model (Trilinos-THCM)
    std::shared_ptr<Ocean>  ocean_;

    //! Parallel implicit atmosphere model
    std::shared_ptr<AtmospherePar> atmos_;

    //! Parallel implicit sea ice model
    std::shared_ptr<SeaIce> seaice_;
    
    //! Combined state vector
    std::shared_ptr<Combined_MultiVec> stateView_;

    //! Combined solution vector
    std::shared_ptr<Combined_MultiVec> solView_;

    //! Combined rhs vectir
    std::shared_ptr<Combined_MultiVec> rhsView_;

    //! Determines the continuation parameter
    //! For a list see THCM::par2int()
    std::string parName_;

    //! Solving scheme
    //!   'D': decoupled     (decoupled, syncs at post-processing)
    //!   'Q': quasi-coupled (no coupling blocks, syncs at every NR step)
    //!   'C': coupled       (fully coupled in FGMRES)
    char solvingScheme_;

    //! Preconditioning
    //!   'S': standard      (1.5 backward block Gauss Seidel)
    //!   'D': diagonal      (do not incorporate coupling blocks in prec)
    char precScheme_;

    //! select whether we should use the atmosphere model in the coupling
    bool useAtmos_;

    //! select whether we should use the sea ice model in the coupling
    bool useSeaIce_;

    //! keep track of syncs
    int syncCtr_;

    //! initialization flag linear solver
    bool solverInitialized_;

    //! Trilinos MPI-like communicator
    Teuchos::RCP<Epetra_Comm> comm_;

    //! Coupling block atmos -> ocean
    CouplingBlock<std::shared_ptr<Ocean>,
                  std::shared_ptr<AtmospherePar> > C12_;
    
    //! Coupling block seaice -> ocean
    CouplingBlock<std::shared_ptr<Ocean>,
                  std::shared_ptr<SeaIce> > C13_;
    
    //! Coupling block ocean -> atmos
    CouplingBlock<std::shared_ptr<AtmospherePar>,
                  std::shared_ptr<Ocean> > C21_;
    
    //! Coupling block seaice -> atmos
    CouplingBlock<std::shared_ptr<AtmospherePar>,
                  std::shared_ptr<SeaIce> > C23_;

    //! Coupling block ocean -> seaice
    CouplingBlock<std::shared_ptr<SeaIce>,
                  std::shared_ptr<Ocean> > C31_;

    //! Coupling block atmosphere -> seaice
    CouplingBlock<std::shared_ptr<SeaIce>,
                  std::shared_ptr<AtmospherePar> > C32_;


    Teuchos::RCP
    <Belos::LinearProblem
     <double, Combined_MultiVec, BelosOp<CoupledModel> > > problem_;

    Teuchos::RCP
    <Belos::BlockGmresSolMgr
     <double, Combined_MultiVec, BelosOp<CoupledModel> > > belosSolver_;

public:
    CoupledModel(std::shared_ptr<Ocean> ocean,
                 std::shared_ptr<AtmospherePar> atmos,
                 Teuchos::RCP<Teuchos::ParameterList> params);

    ~CoupledModel() { INFO("CoupledModel destructor"); }

    //! Compute Jacobian matrix
    void computeJacobian();

    //! Compute RHS
    void computeRHS();

    //! Solve Jx=b
    void solve(std::shared_ptr<Combined_MultiVec> rhs);

    //! Initialize FGMRES (Belos) solver
    void initializeFGMRES();

    //! Apply the Jacobian matrix: out = J*v
    void applyMatrix(Combined_MultiVec const &v, Combined_MultiVec &out);

    void applyMassMat(Combined_MultiVec const &v, Combined_MultiVec &out);

    //! Apply the preconditioning: out = inv(P)*v
    //! In the coupled case: 1.5 step in the block gauss seidel scheme
    //!
    //!        z = inv(P)*x  <=>  P*z = x  -> iterate with
    //!
    //!           [M1 C12] * [z1] = [      ] * [z1] + [x1]
    //!           [    M2]   [z2]   [C21   ]   [z2]   [x2]
    //!
    //! z1^0 = 0 => z2^1 = inv(M2)*x2
    //!          => z1^1 = inv(M1)*(x1-C12*z2^1)
    //!          => z2^1 = inv(M2)*(x2+C21*z1^1)
    void applyPrecon(Combined_MultiVec const &v, Combined_MultiVec &out);

    //! Apply the mass matrix: out=B*v
    void applyMassMatrix(Combined_MultiVec const &v, Combined_MultiVec &out);

    //! Get an RCP to the solution vector
    std::shared_ptr<Combined_MultiVec> getSolution(char mode = 'C');

    //! Get an RCP to the state vector
    std::shared_ptr<Combined_MultiVec> getState(char mode = 'C');

    //! Get an RCP to the rhs vector
    std::shared_ptr<Combined_MultiVec> getRHS(char mode = 'C');

    //! Set the continuation parameter
    void setPar(double value);

    //! Get the continuation parameter
    double getPar();

    //! Pre-processing. A controlling scheme might allow some
    //! pre-processing, for instance at the start of a Newton process.
    void preProcess();

    //! Post-processing. Similarly we can supply some post-processing,
    //! for instance when a Newton process has converged.
    void postProcess();

    //! Gather important continuation data to use in summary file
    std::string const writeData(bool describe = false)
        {
            std::ostringstream datastring;
            if (describe)
            {
                datastring << std::setw(_FIELDWIDTH_/2)
                           << "MV"
                           << std::setw(_FIELDWIDTH_)
                           << "Tol"
                           << ocean_->writeData(describe) << " "
                           << atmos_->writeData(describe) ;            
            }
            else
            {
                datastring.precision(_PRECISION_);
                datastring << std::scientific << std::setw(_FIELDWIDTH_/2)
                           << belosSolver_->getNumIters();
                datastring << std::scientific << std::setw(_FIELDWIDTH_)
                           << belosSolver_->achievedTol();
                
                datastring << ocean_->writeData(describe) << " "
                           << atmos_->writeData(describe) ;
            }

            return datastring.str();
        }

    //! Additional monitor for continuation, return true if
    //! destination reached
    bool monitor() { return false; }

    //! Dump blocks
    void dumpBlocks();

private:

    //! Solve the system using FGMRES
    void FGMRESSolve(std::shared_ptr<Combined_MultiVec> rhs);

    //! For internal use: both modes for the underlying getSolution
    //! methods can be called here.
    std::shared_ptr<Combined_MultiVec> getSolution(char mode1, char mode2);

    //! Compute the residual ||b-Jx||
    double computeResidual(std::shared_ptr<Combined_MultiVec> rhs);

    //! coupled   mode 'C'
    //! decoupled mode 'D'
    double computeResidual(Combined_MultiVec const &rhs,
                           Combined_MultiVec const &x,
                           char mode = 'C');

    //! Synchronize the states between the models that are needed to communicate
    void synchronize();

};

//=============================================================================
// This class wraps a model into an operator suitable for use with Belos.
// We need it for CoupledModel, so we put it here.

template <typename Model>
class BelosOp
{
    Model &model_; // hold a reference to the model
    bool isPrec_;

public:
    BelosOp(Model &model, bool isPrec)
        :
        model_(model),
        isPrec_(isPrec)
        {}

    bool isPrec() const { return isPrec_; }

    void ApplyInverse(Combined_MultiVec const &x, Combined_MultiVec &y) const
        { model_.applyPrecon(x, y); }

    void ApplyMatrix(Combined_MultiVec const &x, Combined_MultiVec &y) const
        { model_.applyMatrix(x, y); }
};

//-----------------------------------------------------------------------------
// Specialize Belos::OperatorTraits using BelosOp
namespace Belos
{
    template <typename Model>
    class OperatorTraits <double, Combined_MultiVec, BelosOp<Model> >
    {
    public:
        static void
        Apply (BelosOp<Model> const &Op,
               Combined_MultiVec const &x,
               Combined_MultiVec &y,
               int trans = 0)
            {
                if ( Op.isPrec() )
                {
                    Op.ApplyInverse(x, y);
                }
                else
                {
                    Op.ApplyMatrix(x, y);
                }
            }
        
        static bool
        HasApplyTranspose (const BelosOp<Model> &Op) { return false; }
    };
}

#endif
