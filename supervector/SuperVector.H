#ifndef SUPERVECTOR_H
#define SUPERVECTOR_H

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>
#include <memory> // shared ptr
#include <vector>
#include <math.h>
#include "GlobalDefinitions.H"
#include "Utils.H"

// This class acts as a normal vector with a few standard methods:
//   length, update, dot, norm, scale, randomize
//
// It is in fact a combined 'super' vector and should be able to contain
// multiple distributed Epetra vectors and serial STL vector containers.
//
// In the future we probably need datamembers of the form:
//   std::vector<std::shared_ptr<std::vector<double> > >
//   std::vector<Teuchos::RCP<Epetra_Vector> > 

class SuperVector
{
	Teuchos::RCP<Epetra_Vector>           oceanVector_;
	std::shared_ptr<std::vector<double> > atmosVector_;

	bool haveOceanVector_;
	bool haveAtmosVector_;
	bool isInitialized_;
	
	int length_;	
	
public:
	// Constructors
	SuperVector(Teuchos::RCP<Epetra_Vector> vector);
	SuperVector(std::shared_ptr<std::vector<double> > vector);
	SuperVector(Teuchos::RCP<Epetra_Vector> vector1,
				std::shared_ptr<std::vector<double> > vector2);

	// Destructor
	~SuperVector();	

	// Get the total length of the combined vector
	int length();

	// Apply vector update:
	// this = scalarA * A + scalarThis * this
	void update(double scalarA,	SuperVector &A, double scalarThis);

	// Apply dot product
	double dot(SuperVector &A);

	// Obtain norm of vector
	// mode: 'S': silent
	//       'V': verbose
	//            prints the separate norms of each vector,
	//            returns complete norm
	double norm(char mode = 'S');

	// Randomize and scale the vector\
	// not implemented for STL vector
	void random(double scale = 1.0);

	// Apply scaling to vectors
	void scale(double scale);

	// Get a pointer to the Trilinos/Ocean vector
	Teuchos::RCP<Epetra_Vector> getOceanVector();

	// Get a pointer to the STL/Atmos vector
	std::shared_ptr<std::vector<double> > getAtmosVector();

	// Print the vectors to various output streams
	void print();

	// Apply a linear mapping from to ocean to the atmosphere and vice versa
	// domain/range: 'O': ocean
	//               'A': atmosphere
	// diagonal elements define the scaling
	// indices define the grid points
	//  (surface temperature points in the ocean)
	void linearTransformation(std::vector<double> &diagonal,
							  std::vector<int> &indices,
							  char domain, char range);
	
private:	
	// Perform a few initializations
	// right now it only initializes length_
	void init();
	
};
#endif
