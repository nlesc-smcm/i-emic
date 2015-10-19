#include "SuperVector.H"
#include <functional> // for std::hash
#include <cstdlib>    // for rand();

//------------------------------------------------------------------
// Default constructor:
SuperVector::SuperVector()
	:
	oceanVector_(Teuchos::null),
	atmosVector_(std::shared_ptr<std::vector<double> >()),
	haveOceanVector_(false),
	haveAtmosVector_(false),
	isInitialized_(false)
{
	init();
}

//------------------------------------------------------------------
// Constructor 1:
SuperVector::SuperVector(Teuchos::RCP<Epetra_Vector> vector)
	:
	oceanVector_(vector),
	atmosVector_(std::shared_ptr<std::vector<double> >()),
	haveOceanVector_(true),
	haveAtmosVector_(false),
	isInitialized_(false)
{
	init();
}

//------------------------------------------------------------------
// Constructor 2:
SuperVector::SuperVector(std::shared_ptr<std::vector<double> > vector)
	:
	oceanVector_(Teuchos::null),
	atmosVector_(vector),
	haveOceanVector_(false),
	haveAtmosVector_(true),
	isInitialized_(false)
{
	init();
}

//------------------------------------------------------------------
// Constructor 3:
SuperVector::SuperVector(Teuchos::RCP<Epetra_Vector> vector1,
						 std::shared_ptr<std::vector<double> > vector2)
	:
	oceanVector_(vector1),
	atmosVector_(vector2),
	haveOceanVector_(true),
	haveAtmosVector_(true),
	isInitialized_(false)
{
	init();
}

//------------------------------------------------------------------
// Copy constructor 
SuperVector::SuperVector(SuperVector const &other)
{
	if (other.haveOceanVector())
	{
		oceanVector_ = Teuchos::rcp
			(new Epetra_Vector(*(other.getOceanVector())));
		haveOceanVector_ = true;
	}
	else
	{
		oceanVector_ = Teuchos::null;
		haveOceanVector_ = false;
	}
	if (other.haveAtmosVector())
	{
		atmosVector_ = std::make_shared<std::vector<double> >
			(*(other.getAtmosVector()));
		haveAtmosVector_ = true;
	}
	else
	{
		atmosVector_ = std::shared_ptr<std::vector<double> >();
		haveAtmosVector_ = false;
	}
	init();
}

//------------------------------------------------------------------
// Destructor
SuperVector::~SuperVector()
{}

//------------------------------------------------------------------
int SuperVector::length() const
{
	if (isInitialized_)
		return length_;
	else
		return 1;
}

//------------------------------------------------------------------
void SuperVector::update(double scalarA, SuperVector const &A, double scalarThis)
{
	if (length_ != A.length())
	{
		ERROR("Wrong dimensions!", __FILE__, __LINE__);
		return;
	}
	
	if (haveOceanVector_)
		oceanVector_->Update(scalarA, *(A.getOceanVector()), scalarThis);
	
	if (haveAtmosVector_)
	{
		for (size_t idx = 0; idx != A.getAtmosVector()->size(); ++idx)
		{
			(*atmosVector_)[idx] =
				scalarA * (*A.getAtmosVector())[idx]
				+ scalarThis * (*atmosVector_)[idx];
		}
	}
	
	if (!haveOceanVector_ && !haveAtmosVector_)
	{
		ERROR("Undefined behaviour ahead!!", __FILE__, __LINE__);
		return;
	}
}

//----------------------------------------------------------------
double SuperVector::dot(SuperVector const &A) const
{
	if (length_ != A.length())
	{
		std::cout << "Wrong dimensions!" << std::endl;
		return 1;
	}
	
	double dot1 = 0;
	if (haveOceanVector_)
		oceanVector_->Dot(*(A.getOceanVector()), &dot1);
	
	double dot2 = 0;
	if (haveAtmosVector_)
		for (size_t idx = 0; idx != A.getAtmosVector()->size(); ++idx)
			dot2 += (*A.getAtmosVector())[idx] * (*atmosVector_)[idx];
	
	return dot1 + dot2;
}

//------------------------------------------------------------------
double SuperVector::norm(char mode) const
{
	if ((mode == 'V') && haveAtmosVector_ && haveOceanVector_)
	{
		haveAtmosVector_ = false;
		INFO(" ||ocean vector|| : " << sqrt(dot(*this)));
		haveAtmosVector_ = true;
		haveOceanVector_ = false;
		INFO(" ||atmos vector|| : " << sqrt(dot(*this)));
		haveOceanVector_ = true;
	}
			
	double nrm2 = 0.0;
	nrm2 = dot(*this);
	return sqrt(nrm2);
}

//------------------------------------------------------------------
void SuperVector::scale(double scale) const
{
	if (haveOceanVector_)
		oceanVector_->Scale(scale);
	if (haveAtmosVector_)
		for (auto &i : *atmosVector_)
			i *= scale;
}

//------------------------------------------------------------------
Teuchos::RCP<Epetra_Vector> SuperVector::getOceanVector() const
{
	if (!haveOceanVector_)
	{
		ERROR("This wrapper does not contain an EpetraVector",
			  __FILE__, __LINE__);
		return Teuchos::null;
	}
	return oceanVector_;
}

//------------------------------------------------------------------
std::shared_ptr<std::vector<double> > SuperVector::getAtmosVector() const
{
	if (!haveAtmosVector_)
	{
		ERROR("This wrapper does not contain a std::vector",
			  __FILE__, __LINE__);
		return std::shared_ptr<std::vector<double> >();
	}
	return atmosVector_;
}

// ------------------------------------------------------------------
void SuperVector::random()
{
	if (haveAtmosVector_)
	{
		for (size_t i = 0; i < atmosVector_->size(); ++i)
			(*atmosVector_)[i] = (std::rand() / (double) RAND_MAX);
	}
	if (haveOceanVector_)
		oceanVector_->Random();
}

//------------------------------------------------------------------
void SuperVector::zero()
{
	zeroAtmos();
	zeroOcean();
}

//------------------------------------------------------------------
void SuperVector::zeroAtmos()
{
	if (haveAtmosVector_)
		atmosVector_ = std::make_shared<std::vector<double> >
			(atmosVector_->size(), 0.0);		
}

//------------------------------------------------------------------
void SuperVector::zeroOcean()
{
	if (haveOceanVector_)
		oceanVector_->PutScalar(0.0);
}

//------------------------------------------------------------------
void SuperVector::print() const
{
	if (haveOceanVector_)
	{
		oceanVector_->Print(*outFile);  // see GlobalDefinitions.H
		oceanVector_->Print(std::cout); // see GlobalDefinitions.H
	}
			
	if (haveAtmosVector_)
	{
		for (auto &it : *atmosVector_)
			std::cout << it << " ";
		std::cout << std::endl;
	}						
}

//------------------------------------------------------------------
void SuperVector::linearTransformation(std::vector<double> const &diagonal,
									   std::vector<int> const &indices,
									   char domain, char range)
{
	if (domain == 'O' && range == 'A')
	{
		TIMER_START("SuperVector: linearTransformation O->A...");
		int dstLength = diagonal.size();
		int srcLength = indices.size();
		
		// Re-initialize atmosVector				
		atmosVector_ = std::make_shared<std::vector<double> >
			(dstLength, 0.0);

		// Get the part of the oceanVector restricted to the supplied indices
		Teuchos::RCP<Epetra_Vector> restricted =
			Utils::RestrictVector(*oceanVector_, indices);		

		// Gather the restricted ocean
		Teuchos::RCP<Epetra_MultiVector> gathered =
			Utils::AllGather(*restricted);
		
		// FullSol should be allocated
		double *fullSol = new double[srcLength];

		// Get the values in the epetraVector
		gathered->ExtractCopy(fullSol, srcLength);

		// Calculate the scaled values in the destination stdVector
		for (size_t i = 0; i != atmosVector_->size(); ++i)
			(*atmosVector_)[i] = diagonal[i] * fullSol[i];
		// Cleanup
		delete fullSol;
		
		TIMER_STOP("SuperVector: linearTransformation O->A...");
	}
	else if (domain == 'A' && range == 'O')
	{
		TIMER_START("SuperVector: linearTransformation A->O...");
		// calculate diagonal scaling
		std::vector<double> values;
		for (size_t i = 0; i != atmosVector_->size(); ++i)
			values.push_back(diagonal[i] * (*atmosVector_)[i]);
				
		// re-initialize oceanVector
		oceanVector_->PutScalar(0.0);
				
		// fill the oceanVector
		oceanVector_->ReplaceGlobalValues(
			atmosVector_->size(),
			&values[0], &indices[0]);
		TIMER_STOP("SuperVector: linearTransformation A->O...");
	}
}
//------------------------------------------------------------------
void SuperVector::linearTransformation(Teuchos::RCP<Epetra_CrsMatrix> mat)
{
	TIMER_START("SuperVector: linearTransformation O->O...");

	Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(oceanVector_->Map()));
	if (haveOceanVector_)
	{
		CHECK_ZERO(mat->Apply(*oceanVector_, *tmp));
		oceanVector_ = tmp;
	}
	else
		WARNING("No oceanVector...", __FILE__, __LINE__);
	
	TIMER_STOP("SuperVector: linearTransformation O->O...");
}

//------------------------------------------------------------------
void SuperVector::linearTransformation(std::shared_ptr<std::map<std::string,
									   std::vector<double> > > mat)
{
	TIMER_START("SuperVector: linearTransformation A->A...");

	int first;
	int last;
	if (haveAtmosVector_)
	{
		std::vector<double> result(atmosVector_->size(), 0.0);
		// Perform matrix vector product
		// 1->0 based... horrible... 
		for (size_t row = 1; row <= atmosVector_->size(); ++row)
		{
			first   = (*mat)["beg"][row-1];
			last    = (*mat)["beg"][row] - 1;
			for (int col = first; col <= last; ++col)
			{
				result[row-1] += (*mat)["ico"][col-1] *
					(*atmosVector_)[(*mat)["jco"][col-1]-1];
			}			
		}
		atmosVector_ = std::make_shared<std::vector<double> >(result);
	}
	else
		WARNING("No atmosVector...", __FILE__, __LINE__);
	
	TIMER_STOP("SuperVector: linearTransformation A->A...");
}

//------------------------------------------------------------------
// Using an XOR and rotate hash on std::hash
std::size_t SuperVector::hash() const
{
	std::size_t seed = 0;
	std::hash<double> double_hash;
	if (haveOceanVector_)
	{
		int numMyElements = oceanVector_->Map().NumMyElements();
		for (int i = 0; i < numMyElements; ++i)
			seed ^= double_hash((*oceanVector_)[i]) + (seed << 6);
	}
	if (haveAtmosVector_)
	{
		for (size_t i = 0; i < atmosVector_->size(); ++i)
			seed ^= double_hash((*atmosVector_)[i]) + (seed << 6);
	}
	return seed;		
}

//------------------------------------------------------------------
void SuperVector::init()
{
	length_ = 0;
	length_ += (haveOceanVector_) ? oceanVector_->GlobalLength() : 0;
	length_ += (haveAtmosVector_) ? atmosVector_->size() : 0;
			
	isInitialized_ = true;
}

//------------------------------------------------------------------
void SuperVector::info()
{
	std::cout << "SuperVector info:" << std::endl;
	std::cout << "  haveOceanVector " << haveOceanVector_ << std::endl;
	std::cout << "  haveAtmosVector " << haveAtmosVector_ << std::endl;
	std::cout << "  length          " << length_ << std::endl;
}
