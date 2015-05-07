#ifndef VECTOR_H
#define VECTOR_H

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "GlobalDefinitions.H"

class Vector
{
	Teuchos::RCP<Epetra_Vector> vector_;
public:
	Vector(Teuchos::RCP<Epetra_Vector> vector)
		:
		vector_(vector)
	 	{
	 		INFO("Vector constructor called");
	 	}
	
	~Vector() { INFO("Vector destructor called" ); }
	
	int Length()  { return vector_->GlobalLength(); }

	// this = scalarThis * this + scalarA * A
	void Update(double scalarA,	Vector &A,
				double scalarThis)
		{
			vector_->Update(scalarA, *(A.GetRCP()), scalarThis);
		}
	
	double Norm()
		{
			double nrm;
			vector_->Norm2(&nrm);
			return nrm;
		}

	void Random(double scale = 1.0)
		{
			vector_->Random();
			vector_->Scale(scale);
		}

	void Scale(double scale) { vector_->Scale(scale); }

	Teuchos::RCP<Epetra_Vector> GetRCP() { return vector_; }	
};

#endif
