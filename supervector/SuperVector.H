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

	// mutable bools: methods may not alter the the vectors but these
	// flags might be allowed to change temporarily
	mutable bool haveOceanVector_;
	mutable bool haveAtmosVector_;
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
	int length() const;

	// Apply vector update:
	// this = scalarA * A + scalarThis * this
	void update(double scalarA,	SuperVector const &A, double scalarThis);
	
	// Apply dot product
	double dot(SuperVector const &A) const;

	// Obtain norm of vector
	// mode: 'S': silent
	//       'V': verbose
	//            prints the separate norms of each vector,
	//            returns complete norm
	double norm(char mode = 'S') const;

	// Apply scaling to vectors
	void scale(double scale) const;

	// Get a pointer to the Trilinos/Ocean vector
	Teuchos::RCP<Epetra_Vector> getOceanVector() const;

	// Get a pointer to the STL/Atmos vector
	std::shared_ptr<std::vector<double> > getAtmosVector() const;

	// Print the vectors to various output streams
	void print() const;

	// Create a hash of the SuperVector
	std::size_t hash() const;

	// Apply a linear mapping from to ocean to the atmosphere and vice versa
	// domain/range: 'O': ocean
	//               'A': atmosphere
	// diagonal elements define the scaling
	// indices define the grid points
	//  (surface temperature points in the ocean)
	void linearTransformation(std::vector<double> const &diagonal,
							  std::vector<int> const &indices,
							  char domain, char range);
	
private:	
	// Perform a few initializations
	// right now it only initializes length_
	void init();
	
};
#endif
