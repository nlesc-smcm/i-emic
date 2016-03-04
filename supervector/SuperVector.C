#include "SuperVector.H"
#include <functional> // for std::hash
#include <cstdlib>    // for rand();
#include <assert.h>
#include "THCMdefs.H"

//------------------------------------------------------------------
// We need some BLAS 
extern "C" void dscal_(int* N, double *DA, double *X, int *INCX);


extern "C" void daxpy_(int* N,    double *DA, double *X,
					   int *INCX, double *Y,  int *INCY);

extern "C" double ddot_(int *N, double *X, int *INCX, double *Y, int *INCY);


//------------------------------------------------------------------
// Default constructor:
SuperVector::SuperVector()
	:
	oceanVector_(Teuchos::null),
	atmosVector_(std::shared_ptr<std::vector<double> >()),
	haveOceanVector_(false),
	haveAtmosVector_(false),
	isInitialized_(false),
	isView_(false)
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
	isInitialized_(false),
	isView_(false)
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
	isInitialized_(false),
	isView_(false)
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
	isInitialized_(false),
	isView_(false)
{
	init();
}

//------------------------------------------------------------------
// Copy constructor, does not maintain view unless specified!
SuperVector::SuperVector(SuperVector const &other)
{
	assign(other);
	init();
}

//------------------------------------------------------------------
// Destructor
SuperVector::~SuperVector()
{}

//------------------------------------------------------------------
// Assignment operator -> copy construction, does not maintain view!
void SuperVector::operator=(SuperVector const &other)
{
	assign(other);
	init();
}

//------------------------------------------------------------------
double &SuperVector::operator[](int index)
{
	return operatorIndex(index);
}

//------------------------------------------------------------------
double const &SuperVector::operator[](int index) const
{
	return operatorIndex(index);
}

//------------------------------------------------------------------
// we assume the ordering {OCEAN, ATMOSPHERE}
double &SuperVector::operatorIndex(int index) const
{
	if ((index < 0) || (index > length_))
		ERROR("INVALID index", __FILE__, __LINE__);
	
	int oceanLength = (haveOceanVector_) ? oceanVector_->GlobalLength() : 0;
	
	if (haveOceanVector_ && index < oceanLength) 
	{
		// Find index in parallel vector
		int lid = oceanVector_->Map().LID(index);
		if (lid >= 0) // means our proc has it
			return (*getOceanVector())[lid];
	}
	
	if (haveAtmosVector_ && index >= oceanLength)  // Serial update
		return (*atmosVector_)[index];
	
	ERROR("UNDEFINED BEHAVIOUR", __FILE__, __LINE__);
	return (*atmosVector_)[index];	
}

//------------------------------------------------------------------
// Assign a copy
void SuperVector::assign(Teuchos::RCP<Epetra_Vector> vector)
{
	oceanVector_ = Teuchos::rcp(new Epetra_Vector(*vector));
	haveOceanVector_ = true;
	init();
}

//------------------------------------------------------------------
// Assign a copy
void SuperVector::assign(std::shared_ptr<std::vector<double> > vector)
{
	atmosVector_ = std::make_shared<std::vector<double> >(*vector);
	haveAtmosVector_ = true;
	init();
}

//------------------------------------------------------------------
//------------------------------------------------------------------
void SuperVector::assign(SuperVector const &other)
{
	TIMER_START("SuperVector: assign");
	if (!other.isView())
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
	}
	else
	{
		if (other.haveOceanVector())
		{
			oceanVector_ = other.getOceanVector();
			haveOceanVector_ = true;
		}
		else
		{
			oceanVector_ = Teuchos::null;
			haveOceanVector_ = false;
		}
		if (other.haveAtmosVector())
		{
			atmosVector_ = other.getAtmosVector();
			haveAtmosVector_ = true;
		}
		else
		{
			atmosVector_ = std::shared_ptr<std::vector<double> >();
			haveAtmosVector_ = false;
		}
	}
	
	isInitialized_ = false;
	TIMER_STOP("SuperVector: assign");
}

//------------------------------------------------------------------
SuperVector SuperVector::createView()
{
	if (haveOceanVector_ && haveAtmosVector_)
	{
		SuperVector view(oceanVector_, atmosVector_);
		view.setView(true);
		return view;
	}
	else if (haveOceanVector_)
	{
		SuperVector view(oceanVector_);
		view.setView(true);
		return view;
		
	}
	else if (haveAtmosVector_)
	{
		SuperVector view(atmosVector_);
		view.setView(true);
		return view;
	}
	else
	{
		SuperVector view;
		view.setView(true);
		return view;
	}
}

//------------------------------------------------------------------
int SuperVector::length() const
{
	if (isInitialized_)
		return length_;
	else
		return -1;
}

//------------------------------------------------------------------
// this = scalarA * A + scalarThis * this

void SuperVector::update(double scalarA, SuperVector const &A, double scalarThis)
{
	assert(length_ == A.length());

	// Update the ocean component
	TIMER_START("SuperVector: update (ocean)");
	if (haveOceanVector_)
		oceanVector_->Update(scalarA, *(A.getOceanVector()), scalarThis);
	TIMER_STOP("SuperVector: update (ocean)");
	
	// Update the atmosphere component
	TIMER_START("SuperVector: update (atmos)");
	if (haveAtmosVector_)
	{
		int N = A.getAtmosVector()->size();
		assert(N == (int) atmosVector_->size());
		
		int incX = 1; int incY = 1;		
		// scale our vector with scalarThis
		dscal_(&N, &scalarThis, &(*atmosVector_)[0], &incY);

		// update 
		daxpy_(&N, &scalarA, &(*A.getAtmosVector())[0],
			   &incX, &(*atmosVector_)[0], &incY);
	}
	TIMER_STOP("SuperVector: update (atmos)");
	
	if (!haveOceanVector_ && !haveAtmosVector_)
	{
		ERROR("Undefined behaviour ahead!!", __FILE__, __LINE__);
		return;
	}
}

//----------------------------------------------------------------
// --> This routine needs some checking I think
void SuperVector::updateElement(int index, double scalar, double scalarThis)
{

	if ((index < 0) || (index > length_))
		ERROR("INVALID index", __FILE__, __LINE__);
	
	int oceanLength = (haveOceanVector_) ? oceanVector_->GlobalLength() : 0;
	
	if (haveOceanVector_ && index < oceanLength) // Parallel update
	{
		int lid = getOceanVector()->Map().LID(index);
		if (lid >= 0)
		{
			(*getOceanVector())[lid] =
				scalarThis * (*getOceanVector())[lid] + scalar;
		}
	}

	if (haveAtmosVector_ && index >= oceanLength)  // Serial update
	{
		int atmosIndex = index - oceanLength;
 		INFO(" 1 update element " << atmosIndex << " "
			 << (*atmosVector_)[atmosIndex]);
		(*atmosVector_)[atmosIndex] =
			scalarThis * (*atmosVector_)[atmosIndex] + scalar;
		INFO(" 2 update element " << atmosIndex << " "
			 << (*atmosVector_)[atmosIndex]);				
	}
	
}

//----------------------------------------------------------------
double SuperVector::dot(SuperVector const &A) const
{
	if (length_ != A.length())
	{
		ERROR("Wrong dimensions!", __FILE__, __LINE__);
		return 1;
	}
	
	TIMER_START("Supervector: dot (ocean)");
	double dot1 = 0;
	if (haveOceanVector_)
		oceanVector_->Dot(*(A.getOceanVector()), &dot1);
	TIMER_STOP("Supervector: dot (ocean)");

	TIMER_START("Supervector: dot (atmos)");
	double dot2 = 0;
	if (haveAtmosVector_)
	{
		int N = A.getAtmosVector()->size();
		assert(N == (int) atmosVector_->size());

		int incX = 1;
		int incY = 1;

		dot2 = ddot_(&N, &(*atmosVector_)[0],
					 &incX, &(*A.getAtmosVector())[0], &incY);
		
		// for (size_t idx = 0; idx != A.getAtmosVector()->size(); ++idx)
		// 	dot2 += (*A.getAtmosVector())[idx] * (*atmosVector_)[idx];
	}
	TIMER_STOP("Supervector: dot (atmos)");
	
	return dot1 + dot2;
}

//------------------------------------------------------------------
double SuperVector::norm(char mode, std::string const msg) const
{
	// Simple verbose component outputting
	if (mode == 'V')
	{
		INFO('\n' << "  " << msg);
		if (haveOceanVector_ && haveAtmosVector_)
		{
			// temporarily disable the other component
			haveAtmosVector_ = false;
			INFO(" ||ocean vector|| : " << sqrt(dot(*this)));
			haveAtmosVector_ = true;
			haveOceanVector_ = false;
			INFO(" ||atmos vector|| : " << sqrt(dot(*this)));
			haveOceanVector_ = true;
		}
		else if (haveOceanVector_)
		{
			INFO(" ||ocean vector|| : " << sqrt(dot(*this)));
		}
		else if (haveAtmosVector_)
		{
			INFO(" ||atmos vector|| : " << sqrt(dot(*this)));
		}
	}
	else if (mode == 'E') // More elaborate component outputting
	{
 		INFO('\n' << "  " << msg);
		componentNorms();
	}
			
	double nrm2 = 0.0;
	nrm2 = dot(*this);
	return sqrt(nrm2);
}

//------------------------------------------------------------------
void SuperVector::componentNorms() const
{
	TIMER_START("SuperVector: component norms");
	double total = 0.0;
	double nrm   = 0.0;
	if (haveOceanVector_)
	{
		// assuming the ordering u,v,w,p,T,S in oceanVector
		int const nun = _NUN_; // 6
		int length    = oceanVector_->GlobalLength();
		std::array<std::vector<int>, nun> indices;
		std::array<char, nun> names = {'u','v','w','p','T','S'};
		for (int i = 0; i != nun; ++i)
		{
			for (int j = 0; j != length / nun; ++j)
				indices[i].push_back(i+j*nun);
			Utils::RestrictVector(*oceanVector_, indices[i])->Norm2(&nrm);
			INFO("  |  " << names[i] << " :  " << nrm);
			total += nrm*nrm;
		}
	}
	if (haveAtmosVector_)
	{
		haveOceanVector_ = false;
		nrm = sqrt(dot(*this));
		total += nrm*nrm;
		INFO("  |  " << "Ta" << ":  " << nrm );
		haveOceanVector_ = true;
	}
	INFO("  |  " << "total = " << sqrt(total) << '\n');
	TIMER_STOP("SuperVector: component norms");
}

//------------------------------------------------------------------
void SuperVector::scale(double scale) const
{
	TIMER_START("Supervector: scale (ocean)");
	if (haveOceanVector_)
		oceanVector_->Scale(scale);
	TIMER_STOP("Supervector: scale (ocean)");

	TIMER_START("Supervector: scale (atmos)");
	if (haveAtmosVector_)
		for (auto &i : *atmosVector_)
			i *= scale;
	TIMER_STOP("Supervector: scale (atmos)");
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
	TIMER_START("SuperVector: zeroAtmos");
	if (haveAtmosVector_)
		atmosVector_->assign(atmosVector_->size(),0.0);
	TIMER_STOP("SuperVector: zeroAtmos");
}

//------------------------------------------------------------------
void SuperVector::zeroOcean()
{
	TIMER_START("SuperVector: zeroOcean");
	if (haveOceanVector_)
		oceanVector_->PutScalar(0.0);
	TIMER_STOP("SuperVector: zeroOcean");
}

//------------------------------------------------------------------
void SuperVector::removeAtmos()
{
	if (haveAtmosVector_)
	{
		atmosVector_ = std::shared_ptr<std::vector<double> >(); // nullptr
		haveAtmosVector_ = false;
	}
}

//------------------------------------------------------------------
void SuperVector::removeOcean()
{
	if (haveOceanVector_)
	{
		oceanVector_ = Teuchos::null; // nullptr
		haveOceanVector_ = false;
	}
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
		{
			INFO(it);
			std::cout << it << " ";
		}
		std::cout << std::endl;
	}						
}

//------------------------------------------------------------------
void SuperVector::print(std::string const filename) const
{
	std::stringstream ocean_fname, atmos_fname;
	ocean_fname << filename << ".ocean";
	atmos_fname << filename << ".atmos";
	
	std::ofstream ocean_ofstream, atmos_ofstream;
	ocean_ofstream.open(ocean_fname.str());
	atmos_ofstream.open(atmos_fname.str());

	if (haveOceanVector_)
		oceanVector_->Print(ocean_ofstream);
	ocean_ofstream.close();
			
	if (haveAtmosVector_)
		for (auto &it : *atmosVector_)
			atmos_ofstream << std::setprecision(12) << it << '\n';

	atmos_ofstream.close();
}

//------------------------------------------------------------------
void SuperVector::linearTransformation(std::vector<double> const &diagonal,
									   std::vector<int> const &indices,
									   char domain, char range)
{
	if (domain == 'O' && range == 'A')
	{
		WARNING("DEPRECATED", __FILE__, __LINE__);
	}
	else if (domain == 'A' && range == 'O')
	{
		TIMER_START("SuperVector: A->O");
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
		TIMER_STOP("SuperVector: A->O");
	}
}
//------------------------------------------------------------------
void SuperVector::linearTransformation(Teuchos::RCP<Epetra_CrsMatrix> mat)
{
	TIMER_START("SuperVector: LINTRANS O->O");

	Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(oceanVector_->Map()));
	if (haveOceanVector_)
	{
		CHECK_ZERO(mat->Apply(*oceanVector_, *tmp));
		oceanVector_ = tmp;
	}
	else
		WARNING("No oceanVector", __FILE__, __LINE__);
	
	TIMER_STOP("SuperVector: LINTRANS O->O");
}

//------------------------------------------------------------------
void SuperVector::linearTransformation(std::shared_ptr<std::map<std::string,
									   std::vector<double> > > mat)
{
	TIMER_START("SuperVector: LINTRANS A->A");

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
		WARNING("No atmosVector", __FILE__, __LINE__);
	
	TIMER_STOP("SuperVector: LINTRANS A->A");
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
	// Right now this is just the computation of our length
	length_ = 0;
	length_ += (haveOceanVector_) ? oceanVector_->GlobalLength() : 0;
	length_ += (haveAtmosVector_) ? atmosVector_->size() : 0;
			
	isInitialized_ = true;
}

//------------------------------------------------------------------
void SuperVector::info() const
{
	std::cout << "*****************************" << std::endl;
	std::cout << " SuperVector info:" << std::endl;
	std::cout << "  haveOceanVector " << haveOceanVector_ << std::endl;
	std::cout << "  haveAtmosVector " << haveAtmosVector_ << std::endl;
	std::cout << "  length          " << length_ << std::endl;
	std::cout << "  norm            " << norm() << std::endl;
}
