#ifndef EPETRA_MATRIXORVECTOR_H
#define EPETRA_MATRIXORVECTOR_H

#include "Teuchos_RCP.hpp"

#include "Epetra_LocalMap.h"
#include "Epetra_SerialComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"

class Epetra_MatrixOrMultiVector: public Epetra_MultiVector
{
public:
    Teuchos::RCP<Epetra_SerialDenseMatrix> mat;

    Epetra_MatrixOrMultiVector(int m, int n)
        :
        Epetra_MultiVector(Epetra_LocalMap(m, 0, Epetra_SerialComm()), n),
        mat(Teuchos::rcp(new Epetra_SerialDenseMatrix(
                             View, Values(), Stride(), MyLength(), NumVectors())))
        {}

    Epetra_MatrixOrMultiVector(Epetra_SerialDenseMatrix const &in)
        :
        Epetra_MultiVector(Copy, Epetra_LocalMap(in.M(), 0, Epetra_SerialComm()),
                           in.A(), in.LDA(), in.N()),
        mat(Teuchos::rcp(new Epetra_SerialDenseMatrix(
                             View, Values(), Stride(), MyLength(), NumVectors())))
        {}

    Epetra_MatrixOrMultiVector(Epetra_MultiVector const &in)
        :
        Epetra_MultiVector(in),
        mat(Teuchos::rcp(new Epetra_SerialDenseMatrix(
                             View, Values(), Stride(), MyLength(), NumVectors())))
        {}

    Epetra_MatrixOrMultiVector(Epetra_MatrixOrMultiVector const &in)
        :
        Epetra_MultiVector(in),
        mat(Teuchos::rcp(new Epetra_SerialDenseMatrix(
                             View, Values(), Stride(), MyLength(), NumVectors())))
        {}

    Epetra_MatrixOrMultiVector() = delete;

    ~Epetra_MatrixOrMultiVector() {}
};

#endif
