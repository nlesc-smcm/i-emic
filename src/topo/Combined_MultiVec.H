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

	Combined_MultiVec(Epetra_DataAccess CV, const Combined_MultiVec &source,
					  const std::vector<int> &index)
		:
		numVecs_(source.numVecs())
		{
			// cast to nonconst for Epetra_MultiVector
			std::vector<int> &tmpInd = const_cast< std::vector<int>& >(index);
			first_  = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.first(),  &tmpInd[0], index.size()));
			second_ = Teuchos::rcp
				(new Epetra_MultiVector(CV, *source.second(), &tmpInd[0], index.size()));
		}	
	
	// Get the rcpointers
	Teuchos::RCP<Epetra_MultiVector> first()  { return first_ ;}
	Teuchos::RCP<Epetra_MultiVector> second() { return second_;}

	// Get number of vectors in each multivector
	int numVecs() {return numVecs_;}
};

//------------------------------------------------------------------
// Specialization of MultiVectorTraits for Belos
//  adapted from BelosEpetraAdapter
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
				
				return Teuchos::rcp(new Combined_MultVec(Copy, mv, index));
			}


				
				
				
	};

	

	
} // end of Belos namespace





#endif
