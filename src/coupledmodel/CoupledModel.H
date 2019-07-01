#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

//! models
#include "Model.H"

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
    using VectorPtr      = std::shared_ptr<Combined_MultiVec>;
    using ConstVectorPtr = std::shared_ptr<const Combined_MultiVec>;
    using LandMask       = Utils::MaskStruct;

protected:

    //! submodel identifiers
    int OCEAN, ATMOS, SEAICE;

    //! vector containing pointers to the different submodels
    std::vector<std::shared_ptr<Model> > models_;

    //! Combined state vector
    std::shared_ptr<Combined_MultiVec> stateView_;

    //! Combined solution vector
    std::shared_ptr<Combined_MultiVec> solView_;

    //! Combined rhs vectir
    std::shared_ptr<Combined_MultiVec> rhsView_;

    //! Solving scheme
    //!   'D': decoupled     (decoupled, syncs at post-processing)
    //!   'Q': quasi-coupled (no coupling blocks, syncs at every NR step)
    //!   'C': coupled       (fully coupled in FGMRES)
    char solvingScheme_;

    //! Preconditioning
    //!   'S': standard      (1.5 backward block Gauss Seidel)
    //!   'D': diagonal      (do not incorporate coupling blocks in prec)
    char precScheme_;

    //! select whether we should use the ocean model in the coupling
    bool useOcean_;

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

    //! Coupling matrix containing off-diagonal CouplingBlocks
    std::vector<std::vector<
                    CouplingBlock<std::shared_ptr<Model>,
                                  std::shared_ptr<Model> > > > C_;

    Teuchos::RCP
    <Belos::LinearProblem
     <double, Combined_MultiVec, BelosOp<CoupledModel> > > problem_;

    Teuchos::RCP
    <Belos::BlockGmresSolMgr
     <double, Combined_MultiVec, BelosOp<CoupledModel> > > belosSolver_;

    double effort_;
    int effortCtr_;

    // gid->coord mapping
    std::vector<std::array<int, 5> > gid2coord_;

public:
    //! constructor
    CoupledModel(std::shared_ptr<Model> ocean,
                 std::shared_ptr<Model> atmos,
                 std::shared_ptr<Model> seaice,
                 Teuchos::RCP<Teuchos::ParameterList> params);

    //! constructor
    CoupledModel(std::shared_ptr<Model> ocean,
                 std::shared_ptr<Model> atmos,
                 Teuchos::RCP<Teuchos::ParameterList> params);

    ~CoupledModel() { INFO("CoupledModel destructor"); }


    //! Set xml parameters
    void setParameters(Teuchos::RCP<Teuchos::ParameterList> params);

    //! CoupledModel setup
    void setup();

    //! Compute Jacobian matrix
    void computeJacobian();

    //! Compute RHS
    void computeRHS();

    //! Solve Jx=b
    void solve(std::shared_ptr<const Combined_MultiVec> rhs);

    //! Initialize FGMRES (Belos) solver
    void initializeFGMRES();

    //! Apply the Jacobian matrix: out = J*v
    void applyMatrix(Combined_MultiVec const &v, Combined_MultiVec &out);

    void applyMassMat(Combined_MultiVec const &v, Combined_MultiVec &out);

    //! Apply the preconditioning: out = inv(P)*v
    //! Example: two models, 1.5 step in the block gauss seidel scheme
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

    //! Set a continuation parameter
    void setPar(std::string const &parName, double value);

    //! Get a continuation parameter
    double getPar(std::string const &parName);

    //! initialize state in submodels
    void initializeState();

    //! Pre-processing. A controlling scheme might allow some
    //! pre-processing, for instance at the start of a Newton process.
    void preProcess();

    //! Post-processing. Similarly we can supply some post-processing,
    //! for instance when a Newton process has converged.
    void postProcess();

    //! Save state to file and everything else
    void saveStateToFile(std::string const &filename);

    //! Gather important continuation data to use in summary file
    std::string const writeData(bool describe = false);

    //! Additional monitor for continuation, return true if
    //! destination reached
    bool monitor() { return false; }

    //! Dump blocks
    void dumpBlocks();

    //! Build GID -> coordinate mapping
    void createGID2CoordMap();

    //! Get GID -> coordinate mapping
    void gid2coord(int const &gid, int &modelIdent,
                   int &i, int &j, int &k, int &xx);

private:

    //! Solve the system using FGMRES
    void FGMRESSolve(std::shared_ptr<const Combined_MultiVec> rhs);

    //! Compute the residual ||b-A*x||
    double explicitResNorm(std::shared_ptr<const Combined_MultiVec> rhs);

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
