#ifndef COMBINED_MULTIVEC
#define COMBINED_MULTIVEC

#include <vector>

#include <Teuchos_RCP.hpp>

#include "Epetra_DataAccess.h"

#include "BelosMultiVec.hpp"

class Epetra_BlockMap;
class Epetra_MultiVector;

//!------------------------------------------------------------------
/*
//! This class combines rcp's to Epetra_MultiVectors and acts as a
//! large vector with the correct traits to use with Belos. The reason
//! for such a class is that we need a unified vector interface for
//! the collection of vectors that come from the submodels within the
//! coupling.

//! We require that the contained MultiVectors contain the same number
//! of ordinary (Epetra) vectors.
*/
//! ------------------------------------------------------------------

class Combined_MultiVec
{
    //! Number of combined Epetra_MultiVectors
    int size_;

    //! Number of vectors in each MultiVector
    int numVecs_;

    //! Pointers to multivectors
    std::vector<Teuchos::RCP<Epetra_MultiVector> > vectors_;

public:
    //! default constructor
    Combined_MultiVec();


    //! constructor using 2 maps
    Combined_MultiVec(const Epetra_BlockMap &map1,
                      const Epetra_BlockMap &map2,
                      int numVectors, bool zeroOut = true);

    //! constructor using 3 maps
    Combined_MultiVec(const Epetra_BlockMap &map1,
                      const Epetra_BlockMap &map2,
                      const Epetra_BlockMap &map3,
                      int numVectors, bool zeroOut = true);

    //! Copy constructor
    Combined_MultiVec(const Combined_MultiVec &source);

    //! constructor using 2 rcp's
    Combined_MultiVec(const Teuchos::RCP<Epetra_MultiVector> &mv1,
                      const Teuchos::RCP<Epetra_MultiVector> &mv2);

    //! constructor using 3 rcp's
    Combined_MultiVec(const Teuchos::RCP<Epetra_MultiVector> &mv1,
                      const Teuchos::RCP<Epetra_MultiVector> &mv2,
                      const Teuchos::RCP<Epetra_MultiVector> &mv3);

    //! constructor using 2 multivecs
    Combined_MultiVec(const Epetra_MultiVector &mv1,
                      const Epetra_MultiVector &mv2);

    //! constructor using 3 multivecs
    Combined_MultiVec(const Epetra_MultiVector &mv1,
                      const Epetra_MultiVector &mv2,
                      const Epetra_MultiVector &mv3);

    //! const
    Combined_MultiVec(Epetra_DataAccess CV, const Combined_MultiVec &source,
                      const std::vector<int> &index);

    //! nonconst
    Combined_MultiVec(Epetra_DataAccess CV, Combined_MultiVec &source,
                      const std::vector<int> &index);

    //! const
    Combined_MultiVec(Epetra_DataAccess CV, const Combined_MultiVec &source,
                      int startIndex, int numVectors);

    //! nonconst
    Combined_MultiVec(Epetra_DataAccess CV, Combined_MultiVec &source,
                      int startIndex, int numVectors);

    //! destructor
    ~Combined_MultiVec()
        {}

    //! Append a copy of an Epetra_MultiVector
    void AppendVector(const Epetra_MultiVector &mv);

    //! Append an rcp to an Epetra_MultiVector
    void AppendVector(const Teuchos::RCP<Epetra_MultiVector> &mv);

    //! Assignment operator, assigning the multivectors
    Combined_MultiVec &operator=(const Combined_MultiVec &source);

    //! Insertion operator
    friend std::ostream &operator<<(std::ostream &out, const Combined_MultiVec &mv);

    //! Index operator non-const
    Teuchos::RCP<Epetra_MultiVector> const &operator()(int index) const;

    //! Index operator non-const
    Teuchos::RCP<Epetra_MultiVector> &operator()(int index);

    //! Local const element access of the first vector in the multivectors
    double const &operator[](int index) const;

    //! Local non-const element access of the first vector in the multivectors
    double &operator[](int index);

    //! Get maps
    const Epetra_BlockMap &Map(int index) const;

    //! Get number of combined multivectors
    int Size() const;

    //! Get number of vectors in each multivector
    int NumVectors() const;

    //! Get the global length of the combined multivector
    int GlobalLength() const;

    //! Get the local length of the combined multivector
    int MyLength() const;

    //! Query the stride
    bool ConstantStride() const;

    //! this = alpha*A*B + scalarThis*this
    int Multiply(char transA, char transB, double scalarAB,
                 const Combined_MultiVec &A, const Epetra_MultiVector &B,
                 double scalarThis);

    //! this = scalarA*A + scalarThis*this
    int Update(double scalarA, const Combined_MultiVec &A, double scalarThis);

    //! this = scalarA*A + scalarB*B + scalarThis*this
    int Update(double scalarA, const Combined_MultiVec &A,
               double scalarB, const Combined_MultiVec &B, double scalarThis);

    //! b[j] := this[j]^T * A[j]
    int Dot(const Combined_MultiVec& A, std::vector<double> &b) const;

    // result[j] := this[j]^T * A[j]
    int Dot(const Combined_MultiVec& A, double *result) const;

    int Scale(double scalarValue);

    int Norm1(std::vector<double> &result) const;

    int Norm2(double *result) const;

    int Norm2(std::vector<double> &result) const;

    int NormInf(double *result) const;

    int NormInf(std::vector<double> &result) const;

    //! direct access to 2-norm
    double Norm() const;

    int SetSeed(unsigned int Seed_in);

    int Random();

    int PutScalar(double alpha);

    size_t hash() const;

    void Print(std::ostream &os) const;
};


namespace Belos
{
template<>
class MultiVecTraits <double, Combined_MultiVec>
{

public:

    static Teuchos::RCP<Combined_MultiVec>
    Clone (const Combined_MultiVec &mv, const int numVecs);

    static Teuchos::RCP<Combined_MultiVec>
    CloneCopy (const Combined_MultiVec &mv);

    static Teuchos::RCP<Combined_MultiVec>
    CloneCopy (const Combined_MultiVec &mv, const std::vector<int> &index);

    static Teuchos::RCP<Combined_MultiVec>
    CloneCopy (const Combined_MultiVec &mv, const Teuchos::Range1D &index);

    static Teuchos::RCP<Combined_MultiVec>
    CloneViewNonConst (Combined_MultiVec &mv, const std::vector<int> &index);

    static Teuchos::RCP<Combined_MultiVec>
    CloneViewNonConst (Combined_MultiVec& mv, const Teuchos::Range1D& index);

    static Teuchos::RCP<const Combined_MultiVec>
    CloneView(const Combined_MultiVec& mv, const std::vector<int>& index);

    static Teuchos::RCP<Combined_MultiVec>
    CloneView(const Combined_MultiVec &mv, const Teuchos::Range1D &index);

    static int  GetVecLength     (const Combined_MultiVec& mv );
    static int  GetGlobalLength  (const Combined_MultiVec& mv );
    static int  GetNumberVecs    (const Combined_MultiVec& mv );
    static bool HasConstantStride(const Combined_MultiVec& mv );

    //! Epetra style (we should compare this with just a bunch of updates)
    static void MvTimesMatAddMv (const double alpha,
                                 const Combined_MultiVec& A,
                                 const Teuchos::SerialDenseMatrix<int,double>& B,
                                 const double beta,
                                 Combined_MultiVec& mv);

    static void
    MvAddMv (const double alpha,
             const Combined_MultiVec& A,
             const double beta,
             const Combined_MultiVec& B,
             Combined_MultiVec& mv);

    static void
    MvScale (Combined_MultiVec& mv,
             const double alpha);

    //! For all columns j of  mv, set mv[j] = alpha[j] * mv[j].
    static void
    MvScale (Combined_MultiVec &mv,
             const std::vector<double> &alpha);
    //! B := alpha * A^T * mv.
    //! Epetra style
    static void MvTransMv(const double alpha, const Combined_MultiVec &A,
                          const Combined_MultiVec &mv, Teuchos::SerialDenseMatrix<int,double> &B);

    //! For all columns j of mv, set b[j] := mv[j]^T * A[j].
    static void
    MvDot (const Combined_MultiVec &mv,
           const Combined_MultiVec &A,
           std::vector<double> &b);

    //! For all columns j of mv, set normvec[j] = norm(mv[j]).
    static void
    MvNorm (const Combined_MultiVec &mv,
            std::vector<double> &normvec,
            NormType type = TwoNorm);

    static void
    SetBlock (const Combined_MultiVec &A,
              const std::vector<int> &index,
              Combined_MultiVec &mv);

    static void
    SetBlock (const Combined_MultiVec &A,
              const Teuchos::Range1D &index,
              Combined_MultiVec &mv);

    static void
    Assign (const Combined_MultiVec& A,
            Combined_MultiVec& mv);

    static void MvRandom (Combined_MultiVec& mv);

    static void MvInit (Combined_MultiVec& mv,
                        double alpha = Teuchos::ScalarTraits<double>::zero());

    static void MvPrint (const Combined_MultiVec& mv, std::ostream& os);

};

}

#endif
