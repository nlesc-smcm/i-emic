#ifndef COMBINED_MULTIVEC
#define COMBINED_MULTIVEC

//------------------------------------------------------------------
// This class combines rcp's to two multivectors
//------------------------------------------------------------------

class Combined_MultiVec
{
	Teuchos::RCP<Epetra_MultiVector> first_;
	Teuchos::RCP<Epetra_MultiVector> second_;

	//! Number of vectors in each multivector
	int numVecs_;

public:
	
	Combined_MultiVec(const Epetra_BlockMap &map1, const Epetra_BlockMap &map2,
		              int numVectors, bool zeroOut = true)
		:
		numVecs_(numVectors)		
		{
			first_   = Teuchos::rcp(new Epetra_MultiVector(Epetra_BlockMap1, numVectors, zeroOut));
			second_  = Teuchos::rcp(new Epetra_MultiVector(Epetra_BlockMap2, numVectors, zeroOut));
		}
	
	// Copy constructor
	Combined_MultiVec(const Combined_MultiVec &source)
		:
		numVecs_(source.numVecs())
		{
			first_   = Teuchos::rcp(new Epetra_MultiVector(*source.first()));
			second_  = Teuchos::rcp(new Epetra_MultiVector(*source.second()));
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
				(new Epetra_MultiVector(CV, *source.first(),  &tmpInd[0], index.size()));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.second(), &tmpInd[0], index.size()));
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
				(new Epetra_MultiVector(CV, *source.first(),  &tmpInd[0], index.size()));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.second(), &tmpInd[0], index.size()));
		}

	// const 
	Combined_MultiVec(Epetra_DataAccess CV, const Combined_MultiVec &source,
					  int startIndex, int numVectors)
		:
		numVecs_(numVectors)
		{
			first_  = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.first(),  startIndex, numVectors));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.second(), startIndex, numVectors));	
		}

	// nonconst
	Combined_MultiVec(Epetra_DataAccess CV, Combined_MultiVec &source,
					  int startIndex, int numVectors)
		:
		numVecs_(numVectors)
		{
			first_  = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.first(),  startIndex, numVectors));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.second(), startIndex, numVectors));	
		}
	
	// Get the rcpointers
	Teuchos::RCP<Epetra_MultiVector> first()  { return first_ ;}
	Teuchos::RCP<Epetra_MultiVector> second() { return second_;}

	// Get number of vectors in each multivector
	int numVecs() {return numVecs_;}

	// Get the length of a vector in the combined multivector
	int globalLength() {return first_->GlobalLength() + second_->GlobalLength();}

	// Query the stride
	bool constantStride() {return first_->ConstantStride() && second_->ConstantStride();}

	
};



//------------------------------------------------------------------
// Specialization of MultiVectorTraits for Belos,
//  adapted from BelosEpetraAdapter.hpp, for better documentation go there.
//------------------------------------------------------------------
namespace Belos
{
	template<>
	class MultiVecTraits<double, Combined_MultiVec>
	{
		static Teuchos::RCP<Combined_MultiVec>
		Clone (const Combined_MultiVec &mv, const int numVecs)
			{
				TEUCHOS_TEST_FOR_EXCEPTION(
					numVecs <= 0, std::invalid_argument,
					"Belos::MultiVecTraits<double, Combined_MultVec>::"
					"Clone(mv, numVecs = " << numVecs << "): "
					"outNumVecs must be positive.");

				return Teuchos::rcp
					(new Combined_MultiVec(mv.first()->Map(),
										   mv.first()->Map(), numVecs, false));
			}

		static Teuchos::RCP<Combined_MultiVec>
		CloneCopy (const Combined_MultiVec &mv)
			{
				return Teuchos::rcp(new Combined_MultiVec(mv));
			}

		static Teuchos::RCP<Combined_MultiVec>
		CloneCopy (const Combined_MultiVec &mv, const std::vector<int> &index)
			{
				const int inNumVecs  = mv.numVecs();
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
				const int inNumVecs   = mv.numVecs();
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
		CloneViewNonConst (Combined_MultiVec &mv, const std::vector<int>& index)
			{
				const int inNumVecs  = mv.numVecs();
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
					index.ubound() < mv.numVecs();
				
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
					TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= mv.numVecs(),
											   std::invalid_argument,
											   os.str() << "Column index range must not exceed "
											   "number of vectors " << mv.numVecs() << " in "
											   "the input multivector.");
				}
				return Teuchos::rcp
					(new Combined_MultiVec(View, mv, index.lbound(), index.size()));
			}

		static Teuchos::RCP<const Combined_MultiVec>
		CloneView (const Combined_MultiVec& mv, const std::vector<int>& index)
			{
				const int inNumVecs  = mv.numVecs();
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
		CloneView (const Combined_& mv, const Teuchos::Range1D& index)
			{
				const bool validRange = index.size() > 0 &&
					index.lbound() >= 0 &&
					index.ubound() < mv.numVecs();
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
					TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= mv.numVecs(),
											   std::invalid_argument,
											   os.str() << "Column index range must not exceed "
											   "number of vectors " << mv.numVecs() << " in "
											   "the input multivector.");
				}
				return Teuchos::rcp (new Combined_MultiVec(View, mv, index.lbound(), index.size()));
			}

		static int  GetVecLength     (const Combined_MultiVec& mv ) { return mv.globalLength(); }
		static int  GetNumberVecs    (const Combined_MultiVec& mv ) { return mv.numVecs(); }
		static bool HasConstantStride(const Combined_MultiVec& mv ) { return mv.constantStride(); }



		// Experimental!! --> HIER VERDER
		static void MvTimesMatAddMv (const double alpha,
									 const Combined_MultiVec& A,
									 const Teuchos::SerialDenseMatrix<int,double>& B,
									 const double beta,
									 Combined_MultiVec& mv)
			{
				int rowsB = B.numRows();
				int colsB = B.numCols();

				// We can do B(i,j) 
			}

		
	};

	
} // end of Belos namespace





#endif
