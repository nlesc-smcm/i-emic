#include "SuperVector.H"

//------------------------------------------------------------------
// Constructor 1:
SuperVector::SuperVector(Teuchos::RCP<Epetra_Vector> vector)
	:
	epetraVector_(vector),
	stdVector_(nullptr),
	haveEpetraVector_(true),
	haveStdVector_(false),
	isInitialized_(false)
{
	init();
}

//------------------------------------------------------------------
// Constructor 2:
SuperVector::SuperVector(std::shared_ptr<std::vector<double> > vector)
	:
	epetraVector_(Teuchos::null),
	stdVector_(vector),
	haveEpetraVector_(false),
	haveStdVector_(true),
	isInitialized_(false)
{
	init();
}

//------------------------------------------------------------------
// Constructor 3:
SuperVector::SuperVector(Teuchos::RCP<Epetra_Vector> vector1,
						 std::shared_ptr<std::vector<double> > vector2)
	:
	epetraVector_(vector1),
	stdVector_(vector2),
	haveEpetraVector_(true),
	haveStdVector_(true),
	isInitialized_(false)
{
	init();
}

//------------------------------------------------------------------
// Destructor
SuperVector::~SuperVector()
{}

//------------------------------------------------------------------
int SuperVector::length()
{
	if (isInitialized_)
		return length_;
	else
		return 1;
}

//------------------------------------------------------------------
void SuperVector::update(double scalarA,	SuperVector &A, double scalarThis)
{
	if (length_ != A.length())
	{
		std::cout << "Wrong dimensions!" << std::endl;
		return;
	}
	if (haveEpetraVector_)
		epetraVector_->Update(scalarA, *(A.getEpetraVector()), scalarThis);
	
	if (haveStdVector_)
	{
		for (size_t idx = 0; idx != A.getStdVector()->size(); ++idx)
		{
			(*stdVector_)[idx] =
				scalarA * (*A.getStdVector())[idx]
				+ scalarThis * (*stdVector_)[idx];
		}
	}
}

//----------------------------------------------------------------
double SuperVector::dot(SuperVector &A)
{
	if (length_ != A.length())
	{
		std::cout << "Wrong dimensions!" << std::endl;
		return 1;
	}
			
	double dot1 = 0;
	if (haveEpetraVector_)
		epetraVector_->Dot(*(A.getEpetraVector()), &dot1);
			
	double dot2 = 0;
	if (haveStdVector_)
		for (size_t idx = 0; idx != A.getStdVector()->size(); ++idx)
			dot2 += (*A.getStdVector())[idx] * (*stdVector_)[idx];
			
	return dot1 + dot2;
}

//------------------------------------------------------------------
double SuperVector::norm(char mode)
{
	if ((mode == 'V') && haveStdVector_ && haveEpetraVector_)
	{
		haveStdVector_ = false;
		INFO("Norm vector 1 (Epetra): " << sqrt(dot(*this)));
		haveStdVector_ = true;
		haveEpetraVector_ = false;
		INFO("Norm vector 2    (STL): " << sqrt(dot(*this)));
		haveEpetraVector_ = true;
	}
			
	double nrm2 = 0.0;
	nrm2 = dot(*this);
	return sqrt(nrm2);
}

//------------------------------------------------------------------
void SuperVector::random(double scale)
{
	if (haveEpetraVector_)
	{
		epetraVector_->Random();
		epetraVector_->Scale(scale);
	}
	if (haveStdVector_)
		std::cout << "Not implemented for std::vector" << std::endl;
}

//------------------------------------------------------------------
void SuperVector::scale(double scale)
{
	if (haveEpetraVector_)
		epetraVector_->Scale(scale);
	if (haveStdVector_)
		for (auto &i : *stdVector_)
			i *= scale;
}

//------------------------------------------------------------------
Teuchos::RCP<Epetra_Vector> SuperVector::getEpetraVector()
{
	if (!haveEpetraVector_)
	{
		ERROR("This wrapper does not contain an EpetraVector",
			  __FILE__, __LINE__);
		return Teuchos::null;
	}
	return epetraVector_;
}

//------------------------------------------------------------------
std::shared_ptr<std::vector<double> > SuperVector::getStdVector()
{
	if (!haveStdVector_)
	{
		ERROR("This wrapper does not contain a std::vector",
			  __FILE__, __LINE__);
		return nullptr;
	}
	return stdVector_;
}

//------------------------------------------------------------------
void SuperVector::print()
{
	if (haveEpetraVector_)
	{
		std::cout << "\nPrinting epetraVector to outFile stream" << std::endl;
		epetraVector_->Print(*outFile);  // see GlobalDefinitions.H
		epetraVector_->Print(std::cout); // see GlobalDefinitions.H
	}
			
	if (haveStdVector_)
	{
		std::cout << "\nPrinting stdVector std::cout" << std::endl;
		for (auto &it : *stdVector_)
			std::cout << it << " ";
		std::cout << std::endl;
	}						
}

//------------------------------------------------------------------
void SuperVector::linearTransformation(std::vector<double> &diagonal,
									   std::vector<int> &indices,
									   char domain, char range)
{
	if (domain == 'O' && range == 'A')
	{
		int dstLength = diagonal.size();
		int srcLength = epetraVector_->GlobalLength();

		// re-initialize stdVector				
		stdVector_ = std::make_shared<std::vector<double> >
			(dstLength, 0.0);

		// gather the epetraVector
		Teuchos::RCP<Epetra_MultiVector> gathered =
			Utils::AllGather(*epetraVector_);

		// fullSol should be allocated
		double *fullSol = new double[srcLength];

		// get the values in the epetraVector
		gathered->ExtractCopy(fullSol, srcLength);

		// calculate the values in the destination stdVector
		for (size_t i = 0; i != stdVector_->size(); ++i)
			(*stdVector_)[i] = diagonal[i] * fullSol[indices[i]];

		// cleanup
		delete fullSol;				
	}
	else if (domain == 'A' && range == 'O')
	{
		std::vector<double> values;
		for (size_t i = 0; i != stdVector_->size(); ++i)
			values.push_back(diagonal[i] * (*stdVector_)[i]);
				
		// re-initialize epetraVector
		epetraVector_->PutScalar(0.0);
				
		// Fill the epetraVector
		epetraVector_->ReplaceGlobalValues(
			stdVector_->size(),
			&values[0], &indices[0]);
	}
}

//------------------------------------------------------------------
void SuperVector::init()
{
	length_ = 0;
	length_ += (haveEpetraVector_) ? epetraVector_->GlobalLength() : 0;
	length_ += (haveStdVector_) ? stdVector_->size() : 0;
	isInitialized_ = true;
}
