#include "SuperVector.H"

#include <functional> // for std::hash

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
		std::cout << "Wrong dimensions!" << std::endl;
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

//------------------------------------------------------------------
void SuperVector::print() const
{
	if (haveOceanVector_)
	{
		std::cout << "\nPrinting epetraVector to outFile stream" << std::endl;
		oceanVector_->Print(*outFile);  // see GlobalDefinitions.H
		oceanVector_->Print(std::cout); // see GlobalDefinitions.H
	}
			
	if (haveAtmosVector_)
	{
		std::cout << "\nPrinting stdVector std::cout" << std::endl;
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
	// Do work
	if (domain == 'O' && range == 'A')
	{
		int dstLength = diagonal.size();
		int srcLength = oceanVector_->GlobalLength();

		// re-initialize stdVector				
		atmosVector_ = std::make_shared<std::vector<double> >
			(dstLength, 0.0);

		// gather the epetraVector
		Teuchos::RCP<Epetra_MultiVector> gathered =
			Utils::AllGather(*oceanVector_);

		// fullSol should be allocated
		double *fullSol = new double[srcLength];

		// get the values in the epetraVector
		gathered->ExtractCopy(fullSol, srcLength);

		// calculate the values in the destination stdVector
		for (size_t i = 0; i != atmosVector_->size(); ++i)
			(*atmosVector_)[i] = diagonal[i] * fullSol[indices[i]];

		// cleanup
		delete fullSol;				
	}
	else if (domain == 'A' && range == 'O')
	{
		std::vector<double> values;
		for (size_t i = 0; i != atmosVector_->size(); ++i)
			values.push_back(diagonal[i] * (*atmosVector_)[i]);
				
		// re-initialize epetraVector
		oceanVector_->PutScalar(0.0);
				
		// Fill the epetraVector
		oceanVector_->ReplaceGlobalValues(
			atmosVector_->size(),
			&values[0], &indices[0]);
	}
}

//------------------------------------------------------------------
// Using an XOR and rotate hash on std::hash
size_t SuperVector::hash() const
{
	size_t seed = 0;
	std::hash<double> double_hash;
	if (haveOceanVector_)
	{
		int numMyElements = oceanVector_->Map().NumMyElements();
		for (int i = 0; i < numMyElements; ++i)
			seed ^= double_hash((*oceanVector_)[i]) + (seed << 6);
	}
	if (haveAtmosVector_)
	{
		for (int i = 0; i < atmosVector_->size(); ++i)
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

