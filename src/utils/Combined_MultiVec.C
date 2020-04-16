#include "Combined_MultiVec.H"

#include <math.h>
#include <vector>

#include <Teuchos_RCP.hpp>
#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"
#include "BelosTypes.hpp"
#include "BelosEpetraAdapter.hpp"

#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "Utils.H"

//! default constructor
Combined_MultiVec::Combined_MultiVec()
    :
    size_(0),
    numVecs_(0),
    vectors_(std::vector<Teuchos::RCP<Epetra_MultiVector> >(0))
{}


//! constructor using 2 maps
Combined_MultiVec::Combined_MultiVec(const Epetra_BlockMap &map1,
                                     const Epetra_BlockMap &map2,
                                     int numVectors, bool zeroOut)
    :
    size_(2),
    numVecs_(numVectors)
{
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(map1, numVectors, zeroOut)) );
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(map2, numVectors, zeroOut)) );
}

//! constructor using 3 maps
Combined_MultiVec::Combined_MultiVec(const Epetra_BlockMap &map1,
                                     const Epetra_BlockMap &map2,
                                     const Epetra_BlockMap &map3,
                                     int numVectors, bool zeroOut)
    :
    size_(3),
    numVecs_(numVectors)
{
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(map1, numVectors, zeroOut)) );
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(map2, numVectors, zeroOut)) );
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(map3, numVectors, zeroOut)) );
}

//! Copy constructor
Combined_MultiVec::Combined_MultiVec(const Combined_MultiVec &source)
    :
    size_(source.Size()),
    numVecs_(source.NumVectors())
{
    for (int i = 0; i != size_; ++i)
        vectors_.push_back(
            Teuchos::rcp(new Epetra_MultiVector(*source(i))) );
}

//! constructor using 2 rcp's
Combined_MultiVec::Combined_MultiVec(const Teuchos::RCP<Epetra_MultiVector> &mv1,
                                     const Teuchos::RCP<Epetra_MultiVector> &mv2)
    :
    size_(2),
    numVecs_(mv1->NumVectors())
{
    assert(mv1->NumVectors() == mv2->NumVectors());

    vectors_.push_back(mv1);
    vectors_.push_back(mv2);
}

//! constructor using 3 rcp's
Combined_MultiVec::Combined_MultiVec(const Teuchos::RCP<Epetra_MultiVector> &mv1,
                                     const Teuchos::RCP<Epetra_MultiVector> &mv2,
                                     const Teuchos::RCP<Epetra_MultiVector> &mv3)
    :
    size_(3),
    numVecs_(mv1->NumVectors())
{
    assert(mv1->NumVectors() == mv2->NumVectors());
    assert(mv2->NumVectors() == mv3->NumVectors());

    vectors_.push_back(mv1);

    vectors_.push_back(mv2);
    vectors_.push_back(mv3);
}

//! constructor using 2 multivecs
Combined_MultiVec::Combined_MultiVec(const Epetra_MultiVector &mv1,
                                     const Epetra_MultiVector &mv2)
    :
    size_(2),
    numVecs_(mv1.NumVectors())
{
    assert(mv1.NumVectors() == mv2.NumVectors());

    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(mv1)) );
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(mv2)) );
}

//! constructor using 3 multivecs
Combined_MultiVec::Combined_MultiVec(const Epetra_MultiVector &mv1,
                                     const Epetra_MultiVector &mv2,
                                     const Epetra_MultiVector &mv3)
    :
    size_(3),
    numVecs_(mv1.NumVectors())
{
    assert(mv1.NumVectors() == mv2.NumVectors());
    assert(mv2.NumVectors() == mv3.NumVectors());

    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(mv1)) );
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(mv2)) );
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(mv3)) );
}

//! const
Combined_MultiVec::Combined_MultiVec(Epetra_DataAccess CV, const Combined_MultiVec &source,
                                     const std::vector<int> &index)
    :
    size_(source.Size()),
    numVecs_(index.size())
{
    //! cast to nonconst for Epetra_MultiVector
    std::vector<int> &tmpInd = const_cast< std::vector<int>& >(index);

    for (int i = 0; i != size_; ++i)
        vectors_.push_back(
            Teuchos::rcp(new Epetra_MultiVector(CV, *source(i),
                                                &tmpInd[0], index.size())) );
}

//! nonconst
Combined_MultiVec::Combined_MultiVec(Epetra_DataAccess CV, Combined_MultiVec &source,
                                     const std::vector<int> &index)
    :
    size_(source.Size()),
    numVecs_(index.size())
{
    //! cast to nonconst for Epetra_MultiVector
    std::vector<int> &tmpInd = const_cast< std::vector<int>& >(index);

    for (int i = 0; i != size_; ++i)
        vectors_.push_back(
            Teuchos::rcp(new Epetra_MultiVector(CV, *source(i),
                                                &tmpInd[0], index.size())) );
}

//! const
Combined_MultiVec::Combined_MultiVec(Epetra_DataAccess CV, const Combined_MultiVec &source,
                                     int startIndex, int numVectors)
    :
    size_(source.Size()),
    numVecs_(numVectors)
{
    for (int i = 0; i != size_; ++i)
        vectors_.push_back(
            Teuchos::rcp(new Epetra_MultiVector(CV, *source(i),
                                                startIndex, numVectors)) );
}

//! nonconst
Combined_MultiVec::Combined_MultiVec(Epetra_DataAccess CV, Combined_MultiVec &source,
                                     int startIndex, int numVectors)
    :
    size_(source.Size()),
    numVecs_(numVectors)
{
    for (int i = 0; i != size_; ++i)
        vectors_.push_back(
            Teuchos::rcp(new Epetra_MultiVector(CV, *source(i),
                                                startIndex, numVectors)) );
}

//! Append a copy of an Epetra_MultiVector
void Combined_MultiVec::AppendVector(const Epetra_MultiVector &mv)
{
    if (size_ >= 1)
        assert(mv.NumVectors() == numVecs_);

    // adjust datamembers
    numVecs_ = mv.NumVectors();
    size_++;
    vectors_.push_back(
        Teuchos::rcp(new Epetra_MultiVector(mv)) );
}

//! Append an rcp to an Epetra_MultiVector
void Combined_MultiVec::AppendVector(const Teuchos::RCP<Epetra_MultiVector> &mv)
{
    if (size_ >= 1)
        assert(mv->NumVectors() == numVecs_);

    // adjust datamembers
    numVecs_ = mv->NumVectors();
    size_++;
    vectors_.push_back(mv);
}

//! Assignment operator, assigning the multivectors
Combined_MultiVec &Combined_MultiVec::operator=(const Combined_MultiVec &source)
{
    // Require that we have the same shape. This might be a
    // bit strong but Epetra_MultiVector asks the same.
    assert(size_    == source.Size());
    assert(numVecs_ == source.NumVectors());

    // Epetra_MultiVector checks whether the shapes of the
    // multivectors are equal.
    for (int i = 0; i != size_; ++i)
        *vectors_[i] = *source(i);

    return *this;
}

//! Insertion operator
std::ostream &operator<<(std::ostream &out, const Combined_MultiVec &mv)
{
    for (int i = 0; i != mv.Size(); ++i)
        out << *mv(i)  << std::endl;

    return out;
}

//! Index operator non-const
Teuchos::RCP<Epetra_MultiVector> const &Combined_MultiVec::operator()(int index) const
{
    assert( (index >= 0) && (index < size_) );
    return vectors_[index];
}

//! Index operator non-const
Teuchos::RCP<Epetra_MultiVector> &Combined_MultiVec::operator()(int index)
{
    assert( (index >= 0) && (index < size_) );
    return vectors_[index];
}

//! Local const element access of the first vector in the multivectors
double const &Combined_MultiVec::operator[](int index) const
{
    // bounds checking
    assert( (index >= 0) && index < MyLength() );

    int end   = 0; // range ending
    int start = 0; // range begin
    for (int i = 0; i != size_; ++i)
    {
        end += vectors_[i]->MyLength();

        if (index < end)
        {
            return (*(*vectors_[i])(0))[index - start];
        }

        start += vectors_[i]->MyLength();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(start == MyLength(), std::logic_error,
                               "undefined behaviour");

    return (*(*vectors_[0])(0))[0];
}

//! Local non-const element access of the first vector in the multivectors
double &Combined_MultiVec::operator[](int index)
{
    const Combined_MultiVec &constMe = *this;
    return const_cast<double &>(constMe.operator[](index));
}

//! Get maps
const Epetra_BlockMap &Combined_MultiVec::Map(int index) const
{
    assert(index >= 0); // input bounds check
    assert(size_ >= 1); // data contents check
    return vectors_[index]->Map();
}

//! Get number of combined multivectors
int Combined_MultiVec::Size() const {return size_;}

//! Get number of vectors in each multivector
int Combined_MultiVec::NumVectors() const {return numVecs_;}

//! Get the global length of the combined multivector
int Combined_MultiVec::GlobalLength() const
{
    int gLength = 0;
    for (int i = 0; i != size_; ++i)
        gLength += vectors_[i]->GlobalLength();

    return gLength;
}

//! Get the local length of the combined multivector
int Combined_MultiVec::MyLength() const
{
    int mLength = 0;
    for (int i = 0; i != size_; ++i)
        mLength += vectors_[i]->MyLength();

    return mLength;
}

//! Query the stride
bool Combined_MultiVec::ConstantStride() const
{
    bool cnstStride = true;
    for (int i = 0; i != size_; ++i)
        cnstStride = cnstStride && vectors_[i]->ConstantStride();

    return cnstStride;
}

//! this = alpha*A*B + scalarThis*this
int Combined_MultiVec::Multiply(char transA, char transB, double scalarAB,
                                const Combined_MultiVec &A, const Epetra_MultiVector &B,
                                double scalarThis )
{
    assert(size_ == A.Size());

    int info = 0;
    for (int i = 0; i != size_; ++i)
        info += vectors_[i]->Multiply(transA, transB,
                                      scalarAB, *A(i), B, scalarThis);
    return info;
}

//! this = scalarA*A + scalarThis*this
int Combined_MultiVec::Update(double scalarA, const Combined_MultiVec &A, double scalarThis)
{
    assert(size_ == A.Size());

    int info = 0;
    for (int i = 0; i != size_; ++i)
        info += vectors_[i]->Update(scalarA, *A(i), scalarThis);

    return info;
}

//! this = scalarA*A + scalarB*B + scalarThis*this
int Combined_MultiVec::Update(double scalarA, const Combined_MultiVec &A,
                              double scalarB, const Combined_MultiVec &B, double scalarThis)
{
    int info = 0;
    for (int i = 0; i != size_; ++i)
        info += vectors_[i]->Update(scalarA,  *A(i),  scalarB, *B(i),  scalarThis);

    return info;
}

//! b[j] := this[j]^T * A[j]
int Combined_MultiVec::Dot(const Combined_MultiVec& A, std::vector<double> &b) const
{
    assert(size_    == A.Size());
    assert(numVecs_ == A.NumVectors());
    assert(numVecs_ <= (int) b.size());

    // reset vector
    std::fill(b.begin(), b.end(), 0.0);

    std::vector<double> tmp = b;
    int info = 0;

    for (int i = 0; i != size_; ++i)
    {
        info += vectors_[i]->Dot(*A(i),  &tmp[0]);
        for (int j = 0; j != numVecs_; ++j)
            b[j] += tmp[j];
    }

    return info;
}

// result[j] := this[j]^T * A[j]
int Combined_MultiVec::Dot(const Combined_MultiVec& A, double *result) const
{
    std::vector<double> tmp(numVecs_, 0.0);
    int info = Dot(A, tmp);
    for (int i = 0; i != numVecs_; ++i)
        result[i] = tmp[i];
    return info;
}

int Combined_MultiVec::Scale(double scalarValue)
{
    int info = 0;
    for (int i = 0; i != size_; ++i)
        info += vectors_[i]->Scale(scalarValue);

    return info;
}

int Combined_MultiVec::Norm1(std::vector<double> &result) const
{
    assert(numVecs_ <= (int) result.size());

    // reset result vector
    std::fill(result.begin(), result.end(), 0.0);

    // copy result vector
    std::vector<double> result_tmp = result;

    int info = 0;
    for (int i = 0; i != size_; ++i)
    {
        info += vectors_[i]->Norm1(&result_tmp[0]);
        for (int j = 0; j != numVecs_; ++j)
            result[j] += result_tmp[j];
    }

    return info;
}

int Combined_MultiVec::Norm2(double *result) const
{
    std::vector<double> tmp(numVecs_, 0.0);
    int info = Norm2(tmp);
    for (int i = 0; i != numVecs_; ++i)
        result[i] = tmp[i];
    return info;
}

int Combined_MultiVec::Norm2(std::vector<double> &result) const
{
    assert(numVecs_ <= (int) result.size());

    // reset result vector
    std::fill(result.begin(), result.end(), 0.0);

    // copy result vector
    std::vector<double> result_tmp = result;

    int info = 0;
    for (int i = 0; i != size_; ++i)
    {
        info += vectors_[i]->Norm2(&result_tmp[0]);
        for (int j = 0; j != numVecs_; ++j)
            result[j] += pow(result_tmp[j], 2);
    }

    // take sqrt of summation per vec in multivec
    for (int j = 0; j != numVecs_; ++j)
        result[j] = sqrt(result[j]);

    return info;
}

int Combined_MultiVec::NormInf(double *result) const
{
    std::vector<double> tmp(numVecs_, 0.0);
    int info = NormInf(tmp);
    for (int i = 0; i != numVecs_; ++i)
        result[i] = tmp[i];
    return info;
}

int Combined_MultiVec::NormInf(std::vector<double> &result) const
{
    assert(numVecs_ <= (int) result.size());

    // reset result vector
    std::fill(result.begin(), result.end(), 0.0);

    // copy result vector
    std::vector<double> result_tmp = result;

    int info = 0;
    for (int i = 0; i != size_; ++i)
    {
        info += vectors_[i]->NormInf(&result_tmp[0]);
        for (int j = 0; j != numVecs_; ++j)
            result[j] = std::max(result[j], result_tmp[j]);
    }

    return info;
}

//! direct access to 2-norm
double Combined_MultiVec::Norm() const { return Utils::norm(this); }

int Combined_MultiVec::SetSeed(unsigned int Seed_in)
{
    int info = 0;
    for (int i = 0; i != size_; ++i)
        info += vectors_[i]->SetSeed(Seed_in);

    return info;
}

int Combined_MultiVec::Random()
{
    int info = 0;
    for (int i = 0; i != size_; ++i)
        info += vectors_[i]->Random();

    return info;
}

int Combined_MultiVec::PutScalar(double alpha)
{
    int info = 0;
    for (int i = 0; i != size_; ++i)
        info += vectors_[i]->PutScalar(alpha);

    return info;
}

size_t Combined_MultiVec::hash() const
{
    size_t seed = 0;
    for (int i = 0; i != size_; ++i)
        seed ^= Utils::hash(vectors_[i]);

    return seed;
}

void Combined_MultiVec::Print(std::ostream &os) const
{
    for (int i = 0; i != size_; ++i)
        vectors_[i]->Print(os);
}

//!------------------------------------------------------------------
//! Specialization of MultiVectorTraits for Belos,
//!  adapted from BelosEpetraAdapter.hpp, for better documentation go there.
//!------------------------------------------------------------------

namespace Belos
{
Teuchos::RCP<Combined_MultiVec>
MultiVecTraits <double, Combined_MultiVec>:: Clone(
    const Combined_MultiVec &mv, const int numVecs)
{
    TEUCHOS_TEST_FOR_EXCEPTION(
        numVecs <= 0, std::invalid_argument,
        "Belos::MultiVecTraits<double, Combined_MultiVec>::"
        "Clone(mv, numVecs = " << numVecs << "): "
        "outNumVecs must be positive.");

    // --> how to generalize this..?
    if (mv.Size() == 2)
    {
        return Teuchos::rcp
            (new Combined_MultiVec(mv(0)->Map(),
                                   mv(1)->Map(),
                                   numVecs, false));
    }
    else if (mv.Size() == 3)
    {
        return Teuchos::rcp
            (new Combined_MultiVec(mv(0)->Map(),
                                   mv(1)->Map(),
                                   mv(2)->Map(),
                                   numVecs, false));
    }
    else
    {
        TEUCHOS_TEST_FOR_EXCEPTION(
            (mv.Size() < 2) || (mv.Size() > 3),
            std::invalid_argument,
            " not yet implemented ");
        return Teuchos::null;
    }
}

Teuchos::RCP<Combined_MultiVec>
MultiVecTraits <double, Combined_MultiVec>::CloneCopy(
    const Combined_MultiVec &mv)
{
    return Teuchos::rcp(new Combined_MultiVec(mv));
}

Teuchos::RCP<Combined_MultiVec>
MultiVecTraits <double, Combined_MultiVec>::CloneCopy(
    const Combined_MultiVec &mv, const std::vector<int> &index)
{
    const int inNumVecs  = mv.NumVectors();
    const int outNumVecs = index.size();
    TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
                               "Belos::MultiVecTraits<double, Combined_MultiVec>::"
                               "CloneCopy(mv, index = {}): At least one vector must be"
                               " cloned from mv.");

    if (outNumVecs > inNumVecs)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double, Combined_Operator>::"
            "CloneCopy(mv, index = {";
        for (int k = 0; k < outNumVecs - 1; ++k)
            os << index[k] << ", ";
        os << index[outNumVecs-1] << "}): There are " << outNumVecs
           << " indices to copy, but only " << inNumVecs << " columns of mv.";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
    }

    return Teuchos::rcp(new Combined_MultiVec(Copy, mv, index));
}

Teuchos::RCP<Combined_MultiVec>
MultiVecTraits <double, Combined_MultiVec>::CloneCopy(
    const Combined_MultiVec &mv, const Teuchos::Range1D &index)
{
    const int inNumVecs   = mv.NumVectors();
    const int outNumVecs  = index.size();
    const bool validRange = outNumVecs > 0 && index.lbound() >= 0 &&
        index.ubound() < inNumVecs;

    if (! validRange)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double, Combined_MultiVec>::Clone(mv,"
            "index=[" << index.lbound() << ", " << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
                                   os.str() << "Column index range must be nonempty.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
                                   os.str() << "Column index range must be nonnegative.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= inNumVecs, std::invalid_argument,
                                   os.str() << "Column index range must not exceed "
                                   "number of vectors " << inNumVecs << " in the "
                                   "input multivector.");
    }

    return Teuchos::rcp
        (new Combined_MultiVec(Copy, mv, index.lbound(), index.size()));

}

Teuchos::RCP<Combined_MultiVec>
MultiVecTraits <double, Combined_MultiVec>::CloneViewNonConst(
    Combined_MultiVec &mv, const std::vector<int> &index)
{
    const int inNumVecs  = mv.NumVectors();
    const int outNumVecs = index.size();
    //! Simple, inexpensive tests of the index vector.

    TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
                               "Belos::MultiVecTraits<double,Combined_MultiVec>::"
                               "CloneViewNonConst(mv, index = {}): The output view "
                               "must have at least one column.");
    if (outNumVecs > inNumVecs)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double,Combined_MultiVec>::"
            "CloneViewNonConst(mv, index = {";
        for (int k = 0; k < outNumVecs - 1; ++k)
            os << index[k] << ", ";
        os << index[outNumVecs-1] << "}): There are " << outNumVecs
           << " indices to view, but only " << inNumVecs << " columns of mv.";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
    }
    return Teuchos::rcp
        (new Combined_MultiVec(View, mv, index));
}

Teuchos::RCP<Combined_MultiVec>
MultiVecTraits <double, Combined_MultiVec>::CloneViewNonConst(
    Combined_MultiVec& mv, const Teuchos::Range1D& index)
{
    const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 &&
        index.ubound() < mv.NumVectors();

    if (! validRange)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double,Combined_MultiVec>::CloneView"
            "NonConst(mv,index=[" << index.lbound() << ", " << index.ubound()
           << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
                                   os.str() << "Column index range must be nonempty.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
                                   os.str() << "Column index range must be nonnegative.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= mv.NumVectors(),
                                   std::invalid_argument,
                                   os.str() << "Column index range must not exceed "
                                   "number of vectors " << mv.NumVectors() << " in "
                                   "the input multivector.");
    }
    return Teuchos::rcp
        (new Combined_MultiVec(View, mv, index.lbound(), index.size()));
}

Teuchos::RCP<const Combined_MultiVec>
MultiVecTraits <double, Combined_MultiVec>::CloneView(
    const Combined_MultiVec& mv, const std::vector<int>& index)
{
    const int inNumVecs  = mv.NumVectors();
    const int outNumVecs = index.size();

    //! Simple, inexpensive tests of the index vector.
    TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
                               "Belos::MultiVecTraits<double,Combined_MultiVec>::"
                               "CloneView(mv, index = {}): The output view "
                               "must have at least one column.");
    if (outNumVecs > inNumVecs)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double,Combined_MultiVec>::"
            "CloneView(mv, index = {";
        for (int k = 0; k < outNumVecs - 1; ++k)
            os << index[k] << ", ";
        os << index[outNumVecs-1] << "}): There are " << outNumVecs
           << " indices to view, but only " << inNumVecs << " columns of mv.";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
    }

    return Teuchos::rcp (new Combined_MultiVec(View, mv, index));
}

Teuchos::RCP<Combined_MultiVec>
MultiVecTraits <double, Combined_MultiVec>::CloneView(
    const Combined_MultiVec &mv, const Teuchos::Range1D &index)
{
    const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 &&
        index.ubound() < mv.NumVectors();
    if (! validRange)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double,Combined_MultiVec>::CloneView"
            "(mv,index=[" << index.lbound() << ", " << index.ubound()
           << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
                                   os.str() << "Column index range must be nonempty.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
                                   os.str() << "Column index range must be nonnegative.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= mv.NumVectors(),
                                   std::invalid_argument,
                                   os.str() << "Column index range must not exceed "
                                   "number of vectors " << mv.NumVectors() << " in "
                                   "the input multivector.");
    }
    return Teuchos::rcp (new Combined_MultiVec(View, mv, index.lbound(), index.size()));
}

int  MultiVecTraits <double, Combined_MultiVec>::GetVecLength(
    const Combined_MultiVec& mv)
{
    return mv.GlobalLength();
}

int  MultiVecTraits <double, Combined_MultiVec>::GetGlobalLength(
    const Combined_MultiVec& mv)
{
    return mv.GlobalLength();
}

int  MultiVecTraits <double, Combined_MultiVec>::GetNumberVecs(
    const Combined_MultiVec& mv)
{
    return mv.NumVectors();
}

bool MultiVecTraits <double, Combined_MultiVec>::HasConstantStride(
    const Combined_MultiVec& mv)
{
    return mv.ConstantStride();
}

//! Epetra style (we should compare this with just a bunch of updates)
void MultiVecTraits <double, Combined_MultiVec>::MvTimesMatAddMv(
    const double alpha, const Combined_MultiVec& A,
    const Teuchos::SerialDenseMatrix<int,double>& B,
    const double beta, Combined_MultiVec& mv)
{
    //! Create Epetra_Multivector from SerialDenseMatrix
    Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map(0).Comm());
    Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());

    const int info = mv.Multiply ('N', 'N', alpha, A, B_Pvec, beta);

    TEUCHOS_TEST_FOR_EXCEPTION(
        info != 0, EpetraMultiVecFailure,
        "Belos::MultiVecTraits<double,Combined_MultiVec>::MvTimesMatAddMv: "
        "Combined_MultiVec::multiply() returned a nonzero value info=" << info
        << ".");
}

void
MultiVecTraits <double, Combined_MultiVec>::MvAddMv(
    const double alpha,const Combined_MultiVec& A,
    const double beta,const Combined_MultiVec& B,
    Combined_MultiVec& mv)
{
    const int info = mv.Update (alpha, A, beta, B, 0.0);

    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
                               "Belos::MultiVecTraits<double, Combined_MultiVec>::MvAddMv: Call to "
                               "update() returned a nonzero value " << info << ".");
}

void
MultiVecTraits <double, Combined_MultiVec>::MvScale(
    Combined_MultiVec& mv, const double alpha)
{
    const int info = mv.Scale(alpha);

    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
                               "Belos::MultiVecTraits<double,Combined_MultiVec>::MvScale: "
                               "Combined_MultiVec::scale() returned a nonzero value info="
                               << info << ".");
}

//! For all columns j of  mv, set mv[j] = alpha[j] * mv[j].
void
MultiVecTraits <double, Combined_MultiVec>::MvScale(
    Combined_MultiVec &mv, const std::vector<double> &alpha)
{
    //! Check to make sure the vector has the same number of entries
    //! as the multivector has columns.
    const int numvecs = mv.NumVectors();

    TEUCHOS_TEST_FOR_EXCEPTION(
        (int) alpha.size () != numvecs, EpetraMultiVecFailure,
        "Belos::MultiVecTraits<double,Combined_MultiVec>::MvScale: "
        "Array alpha of scaling coefficients has " << alpha.size ()
        << " entries, which is not the same as the number of columns "
        << numvecs << " in the input multivector mv.");

    int info = 0;
    int startIndex = 0;
    for (int i = 0; i < numvecs; ++i)
    {
        Combined_MultiVec temp_vec(::View, mv, startIndex, 1);
        info = temp_vec.Scale(alpha[i]);

        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
                                   "Belos::MultiVecTraits<double,Combined_MultiVec>::MvScale: "
                                   "On column " << (i+1) << " of " << numvecs << ", Epetra_Multi"
                                   "Vector::Scale() returned a nonzero value info=" << info << ".");
        startIndex++;
    }
}

//! B := alpha * A^T * mv.
//! Epetra style
void MultiVecTraits <double, Combined_MultiVec>::MvTransMv(
    const double alpha, const Combined_MultiVec &A,
    const Combined_MultiVec &mv, Teuchos::SerialDenseMatrix<int,double> &B)
{
    //! Create Epetra_MultiVector from SerialDenseMatrix
    Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map(0).Comm());
    Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());


    int info = B_Pvec.Multiply('T', 'N', alpha, *A(0), *mv(0), 0.0);
    for (int i = 1; i != mv.Size(); ++i)
        info += B_Pvec.Multiply('T', 'N', alpha, *A(i), *mv(i), 1.0);

    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
                               "Belos::MultiVecTraits<double,Combined_MultiVec>::MvTransMv: "
                               "Combined_MultiVec::multiply() returned a nonzero value info="
                               << info << ".");
}

//! For all columns j of mv, set b[j] := mv[j]^T * A[j].
void
MultiVecTraits <double, Combined_MultiVec>::MvDot(
    const Combined_MultiVec &mv,const Combined_MultiVec &A, std::vector<double> &b)
{
    const int info = mv.Dot(A, b);

    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
                               "Belos::MultiVecTraits<double,Combined_MultiVec>::MvDot: "
                               "Combined_MultiVec::dot() returned a nonzero value info="
                               << info << ".");
}

//! For all columns j of mv, set normvec[j] = norm(mv[j]).
void
MultiVecTraits <double, Combined_MultiVec>::MvNorm(
    const Combined_MultiVec &mv, std::vector<double> &normvec, NormType type)
{
    if ((int) normvec.size() >= mv.NumVectors())
    {
        int info = 0;
        switch( type )
        {
        case ( OneNorm ) :
            info = mv.Norm1(normvec);
            break;
        case ( TwoNorm ) :
            info = mv.Norm2(normvec);
            break;
        case ( InfNorm ) :
            info = mv.NormInf(normvec);
            break;
        default:
            break;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
                                   "Belos::MultiVecTraits<double,Combined_MultiVec>::MvNorm: "
                                   "Combined_MultiVec::Norm() returned a nonzero value info="
                                   << info << ".");
    }
}

void
MultiVecTraits <double, Combined_MultiVec>::SetBlock(
    const Combined_MultiVec &A, const std::vector<int> &index, Combined_MultiVec &mv)
{
    const int inNumVecs  = GetNumberVecs(A);
    const int outNumVecs = index.size();

    if (inNumVecs < outNumVecs)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double,Combined_MultiVec>::"
            "SetBlock(A, mv, index = {";
        if (outNumVecs > 0)
        {
            for (int k = 0; k < outNumVecs - 1; ++k)
                os << index[k] << ", ";
            os << index[outNumVecs-1];
        }
        os << "}): A has only " << inNumVecs << " columns, but there are "
           << outNumVecs << " indices in the index vector.";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
    }

    //! Make a view of the columns of mv indicated by the index std::vector.
    Teuchos::RCP<Combined_MultiVec> mv_view = CloneViewNonConst(mv, index);

    //! View of columns [0, outNumVecs-1] of the source multivector A.
    //! If A has fewer columns than mv_view, then create a view of
    //! the first outNumVecs columns of A.
    Teuchos::RCP<const Combined_MultiVec> A_view;
    if (outNumVecs == inNumVecs)
        A_view = Teuchos::rcpFromRef(A); //! Const, non-owning RCP
    else
        A_view = CloneView(A, Teuchos::Range1D(0, outNumVecs - 1));

    *mv_view = *A_view;
}

void
MultiVecTraits <double, Combined_MultiVec>::SetBlock(
    const Combined_MultiVec &A, const Teuchos::Range1D &index, Combined_MultiVec &mv)
{
    const int numColsA  = A.NumVectors();
    const int numColsMv = mv.NumVectors();

    //! 'index' indexes into mv; it's the index set of the target.
    const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;

    //! We can't take more columns out of A than A has.
    const bool validSource = index.size() <= numColsA;

    if (! validIndex || ! validSource)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double, Combined_MultiVec>::SetBlock"
            "(A, index=[" << index.lbound() << ", " << index.ubound() << "], "
            "mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
                                   os.str() << "Range lower bound must be nonnegative.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= numColsMv, std::invalid_argument,
                                   os.str() << "Range upper bound must be less than "
                                   "the number of columns " << numColsA << " in the "
                                   "'mv' output argument.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.size() > numColsA, std::invalid_argument,
                                   os.str() << "Range must have no more elements than"
                                   " the number of columns " << numColsA << " in the "
                                   "'A' input argument.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
    }

    //! View of columns [index.lbound(), index.ubound()] of the
    //! target multivector mv.  We avoid view creation overhead by
    //! only creating a view if the index range is different than [0,
    //! (# columns in mv) - 1].
    Teuchos::RCP<Combined_MultiVec> mv_view;
    if (index.lbound() == 0 && index.ubound()+1 == numColsMv)
        mv_view = Teuchos::rcpFromRef(mv); //! Non-const, non-owning RCP
    else
        mv_view = CloneViewNonConst(mv, index);

    //! View of columns [0, index.size()-1] of the source multivector
    //! A.  If A has fewer columns than mv_view, then create a view
    //! of the first index.size() columns of A.
    Teuchos::RCP<const Combined_MultiVec> A_view;
    if (index.size() == numColsA)
        A_view = Teuchos::rcpFromRef(A); //! Const, non-owning RCP
    else
        A_view = CloneView(A, Teuchos::Range1D(0, index.size()-1));

    *mv_view = *A_view;
}

void
MultiVecTraits <double, Combined_MultiVec>::Assign(
    const Combined_MultiVec& A, Combined_MultiVec& mv)
{
    const int numColsA  = GetNumberVecs(A);
    const int numColsMv = GetNumberVecs(mv);
    if (numColsA > numColsMv)
    {
        std::ostringstream os;
        os << "Belos::MultiVecTraits<double, Combined_MultiVec>::Assign"
            "(A, mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(numColsA > numColsMv, std::invalid_argument,
                                   os.str() << "Input multivector 'A' has "
                                   << numColsA << " columns, but output multivector "
                                   "'mv' has only " << numColsMv << " columns.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
    }

    //! View of the first [0, numColsA-1] columns of mv.
    Teuchos::RCP<Combined_MultiVec> mv_view;
    if (numColsMv == numColsA)
        mv_view = Teuchos::rcpFromRef(mv); //! Non-const, non-owning RCP
    else //! numColsMv > numColsA
        mv_view = CloneView(mv, Teuchos::Range1D(0, numColsA - 1));

    *mv_view = A;
}

void MultiVecTraits <double, Combined_MultiVec>::MvRandom(Combined_MultiVec& mv)
{
    TEUCHOS_TEST_FOR_EXCEPTION( mv.Random()!= 0, EpetraMultiVecFailure,
                                "Belos::MultiVecTraits<double,Combined_MultiVec>::"
                                "MvRandom() call to random() returned a nonzero value.");
}

void MultiVecTraits <double, Combined_MultiVec>::MvInit(
    Combined_MultiVec& mv, double alpha)
{
    TEUCHOS_TEST_FOR_EXCEPTION( mv.PutScalar(alpha) != 0, EpetraMultiVecFailure,
                                "Belos::MultiVecTraits<double,Combined_MultiVec>::"
                                "MvInit() call to putScalar() returned a nonzero value.");
}

void MultiVecTraits <double, Combined_MultiVec>::MvPrint(
    const Combined_MultiVec& mv, std::ostream& os)
{
    os << mv << std::endl;
}

} //! end of Belos namespace
