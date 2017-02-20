#ifndef COMBINED_MULTIVEC
#define COMBINED_MULTIVEC

#include <math.h>
#include <vector>

#include <Teuchos_RCP.hpp>
#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"
#include "BelosTypes.hpp"

//------------------------------------------------------------------
// This class combines rcp's to two multivectors
//  --> this should generalize to n multivectors
//------------------------------------------------------------------

class Combined_MultiVec
{
	//! Number of vectors in each multivector
	int numVecs_;

	//! Pointers to multivectors
	Teuchos::RCP<Epetra_MultiVector> first_;
	Teuchos::RCP<Epetra_MultiVector> second_;

public:
	// default constructor
	Combined_MultiVec()
		:
		numVecs_(0),
		first_(NULL),
		second_(NULL)
		{}
	
	Combined_MultiVec(const Epetra_BlockMap &map1, const Epetra_BlockMap &map2,
		              int numVectors, bool zeroOut = true)
		:
		numVecs_(numVectors)		
		{
			first_   = Teuchos::rcp(new Epetra_MultiVector(map1, numVectors, zeroOut));
			second_  = Teuchos::rcp(new Epetra_MultiVector(map2, numVectors, zeroOut));
		}
	
	// Copy constructor
	Combined_MultiVec(const Combined_MultiVec &source)
		:
		numVecs_(source.NumVecs())
		{
			first_   = Teuchos::rcp(new Epetra_MultiVector(*source.First()));
			second_  = Teuchos::rcp(new Epetra_MultiVector(*source.Second()));
		}

	Combined_MultiVec(const Teuchos::RCP<Epetra_MultiVector> &mv1,
					  const Teuchos::RCP<Epetra_MultiVector> &mv2)
		:
		numVecs_(std::min(mv1->NumVectors(), mv2->NumVectors())),
		first_(mv1),
		second_(mv2)
		{}
	
	Combined_MultiVec(const Epetra_MultiVector &mv1, const Epetra_MultiVector &mv2)
		:
		numVecs_(std::min(mv1.NumVectors(), mv2.NumVectors()))
		{
			first_  = Teuchos::rcp(new Epetra_MultiVector(mv1));
			second_ = Teuchos::rcp(new Epetra_MultiVector(mv2));
		}
	
	// const
	Combined_MultiVec(Epetra_DataAccess CV, const Combined_MultiVec &source,
					  const std::vector<int> &index)
		:
		numVecs_(index.size())
		{
			// cast to nonconst for Epetra_MultiVector
			std::vector<int> &tmpInd = const_cast< std::vector<int>& >(index);
			first_  = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.First(),  &tmpInd[0], index.size()));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.Second(), &tmpInd[0], index.size()));
		}

	// nonconst
	Combined_MultiVec(Epetra_DataAccess CV, Combined_MultiVec &source,
					  const std::vector<int> &index)
		:
		numVecs_(index.size())
		{
			// cast to nonconst for Epetra_MultiVector
			std::vector<int> &tmpInd = const_cast< std::vector<int>& >(index);
			first_  = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.First(),  &tmpInd[0], index.size()));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.Second(), &tmpInd[0], index.size()));
		}

	// const 
	Combined_MultiVec(Epetra_DataAccess CV, const Combined_MultiVec &source,
					  int startIndex, int numVectors)
		:
		numVecs_(numVectors)
		{
			first_  = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.First(),  startIndex, numVectors));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.Second(), startIndex, numVectors));	
		}

	// nonconst
	Combined_MultiVec(Epetra_DataAccess CV, Combined_MultiVec &source,
					  int startIndex, int numVectors)
		:
		numVecs_(numVectors)
		{
			first_  = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.First(),  startIndex, numVectors));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.Second(), startIndex, numVectors));	
		}

	// Assignment operator
	Combined_MultiVec &operator=(const Combined_MultiVec &Source)
		{
			*first_  = *Source.First();
			*second_ = *Source.Second();
			return *this;
		}

	// Insertion operator
	friend std::ostream &operator<<(std::ostream &out, const Combined_MultiVec &mv)
		{
			out << *mv.First()  << std::endl;
			out << *mv.Second() << std::endl;
			return out;
		}

	// Get maps
	const Epetra_BlockMap &Map1() const {return first_->Map();  }
	const Epetra_BlockMap &Map2() const {return second_->Map(); }
	
	// Get the rcpointers
	Teuchos::RCP<Epetra_MultiVector> First()  { return first_ ;}
	Teuchos::RCP<Epetra_MultiVector> Second() { return second_;}

	Teuchos::RCP<Epetra_MultiVector> First()  const { return first_ ;}
	Teuchos::RCP<Epetra_MultiVector> Second() const { return second_;}

	// Get number of vectors in each multivector
	int NumVecs() const {return numVecs_;}

	// Get the length of a vector in the combined multivector
	int GlobalLength() const {return first_->GlobalLength() + second_->GlobalLength();}

	// Query the stride
	bool ConstantStride() const {return first_->ConstantStride() && second_->ConstantStride();}

	// this = alpha*A*B + scalarThis*this
	int Multiply(char transA, char transB, double scalarAB,
				 const Combined_MultiVec &A, const Epetra_MultiVector &B,
				 double scalarThis )
		{
			int info = 0;
			info =  first_->Multiply(transA, transB, scalarAB, *A.First(), B, scalarThis);
			info += second_->Multiply(transA, transB, scalarAB, *A.Second(), B, scalarThis);
			return info;
		}

	// this = scalarThis*this + scalarA*A + scalarB*B
	int Update(double scalarA, const Combined_MultiVec &A,
			   double scalarB, const Combined_MultiVec &B, double scalarThis)
		{
			int info = 0;
			info =  first_->Update(scalarA,  *A.First(),  scalarB, *B.First(),  scalarThis);
			info += second_->Update(scalarA, *A.Second(), scalarB, *B.Second(), scalarThis);
			return info;
		}
	
	// b[j] := this[j]^T * A[j]
	// Not 100% sure
	int Dot(const Combined_MultiVec& A, std::vector<double> &b) const
		{
			// copy b
			std::vector<double> b2 = b;

			int info = 0;
			info += first_->Dot(*A.First(), &b[0]);
			info += second_->Dot(*A.Second(), &b2[0]);

			for (int i = 0; i != numVecs_; ++i)
				b[i] += b2[i];

			return info;
		}

	int Scale(double scalarValue)
		{
			int info = 0;
			info =  first_->Scale(scalarValue);
			info += second_->Scale(scalarValue);
			return info;
		}

	int Norm1(std::vector<double> &result) const
		{
			// copy result vector
			std::vector<double> result_tmp = result;

			int info = 0;
			info =  first_->Norm1(&result[0]);
			info += second_->Norm1(&result_tmp[0]);

			// combine results
			for (int i = 0; i != numVecs_; ++i)
				result[i] += result_tmp[i];
			return info;
		}

	int Norm2(std::vector<double> &result) const
		{
			// copy result vector
			std::vector<double> result_tmp = result;

			int info = 0;
			info =  first_->Norm2(&result[0]);
			info += second_->Norm2(&result_tmp[0]);

			// combine results
			for (int i = 0; i != numVecs_; ++i)
				result[i] = sqrt(pow(result[i],2) + pow(result_tmp[i],2));
			
			return info;
		}

	int NormInf(std::vector<double> &result) const
		{
			std::vector<double> result_tmp = result;

			int info = 0;
			info =  first_->NormInf(&result[0]);
			info += second_->NormInf(&result_tmp[0]);

			// combine results
			for (int i = 0; i != numVecs_; ++i)
				result[i] = std::max(result[i], result_tmp[i]);
			
			return info;
		}

	int Random()
		{
			int info = 0;
			info += first_->Random();
			info += second_->Random();
			return info;
		}
	
	int PutScalar(double alpha)
		{
			int info = 0;
			info += first_->PutScalar(alpha);
			info += second_->PutScalar(alpha);
			return info;
		}
};


//------------------------------------------------------------------
// Specialization of MultiVectorTraits for Belos,
//  adapted from BelosEpetraAdapter.hpp, for better documentation go there.
//------------------------------------------------------------------
namespace Belos
{
	template<>
	class MultiVecTraits <double, Combined_MultiVec>
	{
		
	public:

		static Teuchos::RCP<Combined_MultiVec>
		Clone (const Combined_MultiVec &mv, const int numVecs)
			{
				TEUCHOS_TEST_FOR_EXCEPTION(
					numVecs <= 0, std::invalid_argument,
					"Belos::MultiVecTraits<double, Combined_MultiVec>::"
					"Clone(mv, numVecs = " << numVecs << "): "
					"outNumVecs must be positive.");

				return Teuchos::rcp
					(new Combined_MultiVec(mv.First()->Map(),
										   mv.First()->Map(), numVecs, false));
			}

		static Teuchos::RCP<Combined_MultiVec>
		CloneCopy (const Combined_MultiVec &mv)
			{
				return Teuchos::rcp(new Combined_MultiVec(mv));
			}

		static Teuchos::RCP<Combined_MultiVec>
		CloneCopy (const Combined_MultiVec &mv, const std::vector<int> &index)
			{
				const int inNumVecs  = mv.NumVecs();
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

		static Teuchos::RCP<Combined_MultiVec>
		CloneCopy (const Combined_MultiVec &mv, const Teuchos::Range1D &index)
			{
				const int inNumVecs   = mv.NumVecs();
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

		static Teuchos::RCP<Combined_MultiVec>
		CloneViewNonConst (Combined_MultiVec &mv, const std::vector<int> &index)
			{
				const int inNumVecs  = mv.NumVecs();
				const int outNumVecs = index.size();
				// Simple, inexpensive tests of the index vector.

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

		static Teuchos::RCP<Combined_MultiVec>
		CloneViewNonConst (Combined_MultiVec& mv, const Teuchos::Range1D& index)
			{
				const bool validRange = index.size() > 0 &&
					index.lbound() >= 0 &&
					index.ubound() < mv.NumVecs();
				
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
					TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= mv.NumVecs(),
											   std::invalid_argument,
											   os.str() << "Column index range must not exceed "
											   "number of vectors " << mv.NumVecs() << " in "
											   "the input multivector.");
				}
				return Teuchos::rcp
					(new Combined_MultiVec(View, mv, index.lbound(), index.size()));
			}

		static Teuchos::RCP<const Combined_MultiVec>
		CloneView(const Combined_MultiVec& mv, const std::vector<int>& index)
			{
				const int inNumVecs  = mv.NumVecs();
				const int outNumVecs = index.size();

				// Simple, inexpensive tests of the index vector.
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

		static Teuchos::RCP<Combined_MultiVec>
		CloneView(const Combined_MultiVec &mv, const Teuchos::Range1D &index)
			{
				const bool validRange = index.size() > 0 &&
					index.lbound() >= 0 &&
					index.ubound() < mv.NumVecs();
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
					TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= mv.NumVecs(),
											   std::invalid_argument,
											   os.str() << "Column index range must not exceed "
											   "number of vectors " << mv.NumVecs() << " in "
											   "the input multivector.");
				}
				return Teuchos::rcp (new Combined_MultiVec(View, mv, index.lbound(), index.size()));
			}
		static int  GetVecLength     (const Combined_MultiVec& mv ) { return mv.GlobalLength(); }
		static int  GetNumberVecs    (const Combined_MultiVec& mv ) { return mv.NumVecs(); }
		static bool HasConstantStride(const Combined_MultiVec& mv ) { return mv.ConstantStride(); }

		// Epetra style (we should compare this with just a bunch of updates)
		static void MvTimesMatAddMv (const double alpha,
									 const Combined_MultiVec& A,
									 const Teuchos::SerialDenseMatrix<int,double>& B,
									 const double beta,
									 Combined_MultiVec& mv)
			{
				// Create Epetra_Multivector from SerialDenseMatrix
				Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map1().Comm());
				Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
				
				const int info = mv.Multiply ('N', 'N', alpha, A, B_Pvec, beta);

				TEUCHOS_TEST_FOR_EXCEPTION(
					info != 0, EpetraMultiVecFailure,
					"Belos::MultiVecTraits<double,Combined_MultiVec>::MvTimesMatAddMv: "
					"Combined_MultiVec::multiply() returned a nonzero value info=" << info
					<< ".");
			}
		
		static void
		MvAddMv (const double alpha,
				 const Combined_MultiVec& A,
				 const double beta,
				 const Combined_MultiVec& B,
				 Combined_MultiVec& mv)
			{
				const int info = mv.Update (alpha, A, beta, B, 0.0);								

				TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
										   "Belos::MultiVecTraits<double, Combined_MultiVec>::MvAddMv: Call to "
										   "update() returned a nonzero value " << info << ".");
			}

		static void
		MvScale (Combined_MultiVec& mv,
				 const double alpha)
			{
				const int info = mv.Scale(alpha);

				TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
										   "Belos::MultiVecTraits<double,Combined_MultiVec>::MvScale: "
										   "Combined_MultiVec::scale() returned a nonzero value info="
										   << info << ".");
			}

		//! For all columns j of  mv, set mv[j] = alpha[j] * mv[j].
		static void
		MvScale (Combined_MultiVec &mv,
				 const std::vector<double> &alpha)
			{
				// Check to make sure the vector has the same number of entries
				// as the multivector has columns.
				const int numvecs = mv.NumVecs();

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
		static void MvTransMv(const double alpha, const Combined_MultiVec &A,
							  const Combined_MultiVec &mv, Teuchos::SerialDenseMatrix<int,double> &B)
			{
				// Create Epetra_MultiVector from SerialDenseMatrix
				Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map1().Comm());
				Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
				int info = B_Pvec.Multiply('T', 'N', alpha, *A.First(), *mv.First(), 0.0);
				info    += B_Pvec.Multiply('T', 'N', alpha, *A.Second(), *mv.Second(), 1.0);
				
				TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
										   "Belos::MultiVecTraits<double,Combined_MultiVec>::MvTransMv: "
										   "Combined_MultiVec::multiply() returned a nonzero value info="
										   << info << ".");				
			}
		
		//! For all columns j of mv, set b[j] := mv[j]^T * A[j].
		static void
		MvDot (const Combined_MultiVec &mv,
			   const Combined_MultiVec &A,
			   std::vector<double> &b)
			{
				const int info = mv.Dot(A, b);
				
				TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraMultiVecFailure,
										   "Belos::MultiVecTraits<double,Combined_MultiVec>::MvDot: "
										   "Combined_MultiVec::dot() returned a nonzero value info="
										   << info << ".");
			}

		//! For all columns j of mv, set normvec[j] = norm(mv[j]).
		static void
		MvNorm (const Combined_MultiVec &mv,
				std::vector<double> &normvec,
				NormType type = TwoNorm)
			{
				if ((int) normvec.size() >= mv.NumVecs())
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

		static void
		SetBlock (const Combined_MultiVec &A,
				  const std::vector<int> &index,
				  Combined_MultiVec &mv)
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
				
				// Make a view of the columns of mv indicated by the index std::vector.
				Teuchos::RCP<Combined_MultiVec> mv_view = CloneViewNonConst(mv, index);

				// View of columns [0, outNumVecs-1] of the source multivector A.
				// If A has fewer columns than mv_view, then create a view of
				// the first outNumVecs columns of A.
				Teuchos::RCP<const Combined_MultiVec> A_view;
				if (outNumVecs == inNumVecs)
					A_view = Teuchos::rcpFromRef(A); // Const, non-owning RCP
				else
					A_view = CloneView(A, Teuchos::Range1D(0, outNumVecs - 1));

				*mv_view = *A_view;
			}
		
		static void
		SetBlock (const Combined_MultiVec &A,
				  const Teuchos::Range1D &index,
				  Combined_MultiVec &mv)
			{
				const int numColsA  = A.NumVecs();
				const int numColsMv = mv.NumVecs();

				// 'index' indexes into mv; it's the index set of the target.
				const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;
				
				// We can't take more columns out of A than A has.
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

				// View of columns [index.lbound(), index.ubound()] of the
				// target multivector mv.  We avoid view creation overhead by
				// only creating a view if the index range is different than [0,
				// (# columns in mv) - 1].
				Teuchos::RCP<Combined_MultiVec> mv_view;
				if (index.lbound() == 0 && index.ubound()+1 == numColsMv)
					mv_view = Teuchos::rcpFromRef(mv); // Non-const, non-owning RCP
				else
					mv_view = CloneViewNonConst(mv, index);

				// View of columns [0, index.size()-1] of the source multivector
				// A.  If A has fewer columns than mv_view, then create a view
				// of the first index.size() columns of A.
				Teuchos::RCP<const Combined_MultiVec> A_view;
				if (index.size() == numColsA)
					A_view = Teuchos::rcpFromRef(A); // Const, non-owning RCP
				else
					A_view = CloneView(A, Teuchos::Range1D(0, index.size()-1));

				*mv_view = *A_view;
			}

		static void
		Assign (const Combined_MultiVec& A,
				Combined_MultiVec& mv)
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
				
				// View of the first [0, numColsA-1] columns of mv.
				Teuchos::RCP<Combined_MultiVec> mv_view;
				if (numColsMv == numColsA)
					mv_view = Teuchos::rcpFromRef(mv); // Non-const, non-owning RCP
				else // numColsMv > numColsA
					mv_view = CloneView(mv, Teuchos::Range1D(0, numColsA - 1));

				*mv_view = A;
			}
		
		static void MvRandom (Combined_MultiVec& mv)
			{
				TEUCHOS_TEST_FOR_EXCEPTION( mv.Random()!= 0, EpetraMultiVecFailure,
											"Belos::MultiVecTraits<double,Combined_MultiVec>::"
											"MvRandom() call to random() returned a nonzero value.");
			}

		static void MvInit (Combined_MultiVec& mv,
							double alpha = Teuchos::ScalarTraits<double>::zero())
			{
				TEUCHOS_TEST_FOR_EXCEPTION( mv.PutScalar(alpha) != 0, EpetraMultiVecFailure,
											"Belos::MultiVecTraits<double,Combined_MultiVec>::"
											"MvInit() call to putScalar() returned a nonzero value.");
			}

		static void MvPrint (const Combined_MultiVec& mv, std::ostream& os)
			{
				os << mv << std::endl;
			}
		
	}; // end of specialization

	
} // end of Belos namespace
#endif
