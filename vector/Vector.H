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
	 	{ DEBUG("Vector constructor called"); }
	
	~Vector() { DEBUG("Vector destructor called" ); }
	
	int length()  { return vector_->GlobalLength(); }

	// this = scalarA * A + scalarThis * this
	void update(double scalarA,	Vector &A, double scalarThis)
		{ vector_->Update(scalarA, *(A.getRCP()), scalarThis); }
	
	double dot(Vector &A)
		{
			double dot;
			vector_->Dot(*(A.getRCP()), &dot);
			return dot;
		}
	
	double norm()
		{
			double nrm;
			vector_->Norm2(&nrm);
			return nrm;
		}

	void random(double scale = 1.0)
		{
			vector_->Random();
			vector_->Scale(scale);
		}

	void scale(double scale) { vector_->Scale(scale); }

	Teuchos::RCP<Epetra_Vector> getRCP() { return vector_; }	
};

#endif
