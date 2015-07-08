#ifndef SUPERVECTOR_H
#define SUPERVECTOR_H

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "Teuchos_RCP.hpp"
#include <memory> // shared ptr
#include <vector>
#include <math.h>
#include "GlobalDefinitions.H"

// This class acts as a normal vector with a few standard methods:
//   length, update, dot, norm, scale, randomize
// It is in fact a combined 'super' vector and should be able to contain
// multiple distributed Epetra vectors and serial STL vector containers.

// In the future we probably need datamembers of the form:
//   std::vector<std::shared_ptr<std::vector<double> > >
//   std::vector<Teuchos::RCP<Epetra_Vector> > 

class SuperVector
{
	Teuchos::RCP<Epetra_Vector> epetraVector_;
	std::shared_ptr<std::vector<double> > stdVector_;

	bool haveEpetraVector_;
	bool haveStdVector_;

	int length_;
	
	
public:
	SuperVector(Teuchos::RCP<Epetra_Vector> vector)
		:
		epetraVector_(vector),
		stdVector_(nullptr),
		haveEpetraVector_(true),
		haveStdVector_(false)
		{
			init();
		}

	//------------------------------------------------------------------
	SuperVector(std::shared_ptr<std::vector<double> > vector)
		:
		epetraVector_(Teuchos::null),
		stdVector_(vector),
		haveEpetraVector_(false),
		haveStdVector_(true)
		{
			init();
		}
	
	//------------------------------------------------------------------
	SuperVector(Teuchos::RCP<Epetra_Vector> vector1,
		   std::shared_ptr<std::vector<double> > vector2)
		:
		epetraVector_(vector1),
		stdVector_(vector2),
		haveEpetraVector_(true),
		haveStdVector_(true)
		{
			init();
		}
	
	//------------------------------------------------------------------
	~SuperVector()
		{}
	
	//------------------------------------------------------------------
	void init()
		{
			length_ = 0;
			length_ += (haveEpetraVector_) ? epetraVector_->GlobalLength() : 0;
			length_ += (haveStdVector_) ? stdVector_->size() : 0;
		}

	//------------------------------------------------------------------
	int length() { return length_; }

	// this = scalarA * A + scalarThis * this
	void update(double scalarA,	SuperVector &A, double scalarThis)
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

	//------------------------------------------------------------------
	double dot(SuperVector &A)
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
	// mode: 'S': silent
	//       'V': verbose
	double norm(char mode = 'S')
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
	void random(double scale = 1.0)
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
	void scale(double scale)
		{
			if (haveEpetraVector_)
				epetraVector_->Scale(scale);
			if (haveStdVector_)
				for (auto &i : *stdVector_)
					i *= scale;
		}

	//------------------------------------------------------------------
	Teuchos::RCP<Epetra_Vector> getEpetraVector()
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
	std::shared_ptr<std::vector<double> > getStdVector()
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
	void print()
		{
			if (haveEpetraVector_)
			{
				std::cout << "Printing epetraVector to outFile stream" << std::endl;
				epetraVector_->Print(*outFile); // see GlobalDefinitions.H
			}
			
			if (haveStdVector_)
			{
				for (auto &it : *stdVector_)
					std::cout << it << " ";
				std::cout << std::endl;
			}						
		}

	//------------------------------------------------------------------
	void linearTransformation(std::shared_ptr<std::vector<double> > > diagonal,
							  char domain, char range)
		{
			if (domain == 'O' && range == 'A')
			{
				// lookup Trilinos stuff

				haveEpetraVector_ = false;
			}
			else if (domain == 'A' && range == 'O')
			{

				// ...

				haveStdVector_ = false;
			}
		}
};
#endif
