#include "GlobalDefinitions.H"

#include "Epetra_MultiVectorWrapper.hpp"
#include "Epetra_SerialDenseMatrixWrapper.hpp"
#include "Epetra_OperatorWrapper.hpp"
#include "LyapunovSolver.hpp"
#include "SchurOperator.hpp"
#include "Timer.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialDenseMatrix.h"

#include "AnasaziTypes.hpp"

template<typename Model>
class LyapunovModel: public Model
{
    double trace_;
    std::vector<double> eigenvalues_;
    Teuchos::RCP<Epetra_MultiVector> V_;
public:
    using Model::Model;
    LyapunovModel(Model const &model);

    void postProcess();
    std::string writeData(bool describe = false) const;

    int computeCovarianceMatrix();

    std::vector<double> getEigenvalues();
};


template<typename Model>
LyapunovModel<Model>::LyapunovModel(Model const &model)
    :
    Model(model)
{}

template<typename Model>
int LyapunovModel<Model>::computeCovarianceMatrix()
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("lyapunov_params.xml", params.ptr());

    Model::computeJacobian();
    auto A = Model::getJacobian();

    Model::computeForcing();
    auto B = Model::getForcing();

    CHECK_ZERO(B->Scale(params->get("sigma", 1.0)));

    auto M = Model::getMassMat('V');

    Epetra_Map const &map = A->RowMap();
    Teuchos::RCP<Epetra_MultiVector> nullspace = Teuchos::rcp(new Epetra_MultiVector(map, 2));

    for (int i = 0; i < nullspace->MyLength(); i++)
    {
        int gid = map.GID(i);
        if (gid % 6 == 3)
        {
            int m = gid / 6;
            if ((m % Model::N_ + (m / Model::N_) % Model::M_) % 2 == 0)
                (*nullspace)[0][i] = 1;
            else
                (*nullspace)[1][i] = 1;
        }
    }

    Teuchos::RCP<RAILS::SchurOperator> Schur = Teuchos::rcp(
        new RAILS::SchurOperator(A, M, nullspace));
    Schur->Compute();

    Epetra_Map const &map2 = Schur->OperatorRangeMap();
    Epetra_Import import(map2, map);

    int MaxNumEntriesPerRow = B->MaxNumEntries();
    Teuchos::RCP<Epetra_CrsMatrix> B22 = Teuchos::rcp(
        new Epetra_CrsMatrix(Copy, map2, B->ColMap(), MaxNumEntriesPerRow));
    CHECK_ZERO(B22->Import(*B, import, Insert));
    CHECK_ZERO(B22->FillComplete(B->DomainMap(), map2));

    Teuchos::RCP<Epetra_Operator> Schur_operator = Schur;
    RAILS::Epetra_OperatorWrapper Schur_wrapper = Schur_operator;
    Teuchos::RCP<Epetra_Operator> B22_operator = B22;

    RAILS::Solver<RAILS::Epetra_OperatorWrapper, RAILS::Epetra_MultiVectorWrapper,
                  RAILS::Epetra_SerialDenseMatrixWrapper> solver(
                      Schur_wrapper, B22_operator, B22_operator);

    RAILS::Epetra_MultiVectorWrapper V;
    RAILS::Epetra_SerialDenseMatrixWrapper T;

    if (V_ != Teuchos::null)
    {
        V = V_;
        if (!params->sublist("Lyapunov Solver").isParameter("Restart from solution"))
            params->sublist("Lyapunov Solver").set("Restart from solution", true);
    }
    else
        V = Teuchos::rcp(
            new Epetra_MultiVector(map2, 1000));

    T = Teuchos::rcp(
        new Epetra_SerialDenseMatrix(1000, 1000));

    solver.set_parameters(params->sublist("Lyapunov Solver"));

    solver.solve(V, T);

    if (!A->Comm().MyPID())
        RAILS_SAVE_PROFILES("");

    V_ = Teuchos::rcp(new Epetra_MultiVector(Copy, *V, 0, V.N()));

    Schur->SetSolution(*V, *T);
    Schur_wrapper.set_parameters(*params);

    RAILS::Epetra_MultiVectorWrapper eigenvectors;
    RAILS::Epetra_SerialDenseMatrixWrapper eigenvalues(0,0);

    Teuchos::ParameterList &eig_params = params->sublist("Eigenvalue Solver");
    eig_params.set("Verbosity", Anasazi::Errors +
                   // Anasazi::IterationDetails +
                   Anasazi::Warnings +
                   Anasazi::FinalSummary);

    Schur_wrapper.eigs(eigenvectors, eigenvalues,
                       eig_params.get("Number of Eigenvalues", 10));

    trace_ = Schur->Trace();

    int num_eigs = eigenvalues.M();
    eigenvalues_.clear();
    for (int i = 0; i < num_eigs; i++)
        eigenvalues_.push_back(eigenvalues(i));

    return 0;
}

template<typename Model>
std::vector<double> LyapunovModel<Model>::getEigenvalues()
{
    return eigenvalues_;
}

template<typename Model>
void LyapunovModel<Model>::postProcess()
{
    Model::postProcess();

    computeCovarianceMatrix();

    for (double eig: eigenvalues_)
        std::cout << std::setw(20) << eig
                  << std::setw(20) << eig / trace_
                  << std::endl;
}

template<typename Model>
std::string LyapunovModel<Model>::writeData(bool describe) const
{
    std::ostringstream datastring;
    datastring << Model::writeData(describe);

    if (describe)
    {
        datastring << std::setw(_FIELDWIDTH_)
                   << "max(ev)";

        return datastring.str();
    }
    else
    {
        datastring.precision(_PRECISION_);
        datastring << std::scientific << std::setw(_FIELDWIDTH_)
                   << eigenvalues_[0];

        return datastring.str();
    }
}
