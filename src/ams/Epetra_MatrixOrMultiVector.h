#ifndef EPETRA_MATRIXORVECTOR_H
#define EPETRA_MATRIXORVECTOR_H

#include "Teuchos_RCP.hpp"

#include "Epetra_LocalMap.h"
#include "Epetra_SerialComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"

class Epetra_MatrixOrMultiVector: public Epetra_MultiVector
{
    bool factored;
    Teuchos::RCP<Epetra_SerialDenseMatrix> fac;
    Teuchos::RCP<Epetra_SerialDenseSolver> sol;
public:
    Teuchos::RCP<Epetra_SerialDenseMatrix> mat;

    Epetra_MatrixOrMultiVector(int m, int n)
        :
        Epetra_MultiVector(Epetra_LocalMap(m, 0, Epetra_SerialComm()), n),
        factored(false),
        mat(Teuchos::rcp(new Epetra_SerialDenseMatrix(
                             View, Values(), Stride(), MyLength(), NumVectors())))
        {}

    Epetra_MatrixOrMultiVector(Epetra_SerialDenseMatrix const &in)
        :
        Epetra_MultiVector(Copy, Epetra_LocalMap(in.M(), 0, Epetra_SerialComm()),
                           in.A(), in.LDA(), in.N()),
        factored(false),
        mat(Teuchos::rcp(new Epetra_SerialDenseMatrix(
                             View, Values(), Stride(), MyLength(), NumVectors())))
        {}

    Epetra_MatrixOrMultiVector(Epetra_MultiVector const &in)
        :
        Epetra_MultiVector(in),
        factored(false),
        mat(Teuchos::rcp(new Epetra_SerialDenseMatrix(
                             View, Values(), Stride(), MyLength(), NumVectors())))
        {}

    Epetra_MatrixOrMultiVector(Epetra_MatrixOrMultiVector const &in)
        :
        Epetra_MultiVector(in),
        factored(false),
        mat(Teuchos::rcp(new Epetra_SerialDenseMatrix(
                             View, Values(), Stride(), MyLength(), NumVectors())))
        {}

    Epetra_MatrixOrMultiVector() = delete;

    ~Epetra_MatrixOrMultiVector() {}

    int Factor()
        {
            int ierr = 0;
            sol = Teuchos::rcp(new Epetra_SerialDenseSolver());
            fac = Teuchos::rcp(new Epetra_SerialDenseMatrix(*mat));
            ierr = sol->SetMatrix(*fac);
            EPETRA_CHK_ERR(ierr);
            ierr = sol->Invert();
            EPETRA_CHK_ERR(ierr);
            factored = true;
            return ierr;
        }

    int ApplyInverse(Epetra_MatrixOrMultiVector const &X, Epetra_MatrixOrMultiVector &Y)
        {
            int ierr = 0;
            if (!factored)
            {
                ierr = Factor();
                EPETRA_CHK_ERR(ierr);
            }

            ierr = sol->SetVectors(*Y.mat, *X.mat);
            EPETRA_CHK_ERR(ierr);

            ierr = sol->Solve();
            return ierr;
        }
};

#endif
