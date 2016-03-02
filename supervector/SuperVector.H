#ifndef SUPERVECTOR_H
#define SUPERVECTOR_H

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>
#include <memory> // shared ptr
#include <vector>
#include <string>
#include <map>
#include <math.h>
#include "GlobalDefinitions.H"
#include "Utils.H"

// This class acts as a normal vector with a few standard methods:
//   length, update, dot, norm, scale, random, zero
//
// It is in fact a combined 'super' vector and should be able to contain
// multiple distributed Epetra vectors and serial STL vector containers.
//
// In the future we might need datamembers of the form:
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
	// Default constructor
	SuperVector();
	
	// Constructors
	SuperVector(Teuchos::RCP<Epetra_Vector> vector);
	SuperVector(std::shared_ptr<std::vector<double> > vector);
	SuperVector(Teuchos::RCP<Epetra_Vector> vector1,
				std::shared_ptr<std::vector<double> > vector2);

	// Copy constructor
	SuperVector(SuperVector const &other);

	// Overloaded assignment operator
	void operator=(SuperVector const &other);
	
	// Destructor
	~SuperVector();
	
	// Get info on contents SuperVector
	bool haveOceanVector() const { return haveOceanVector_; }
	bool haveAtmosVector() const { return haveAtmosVector_; }

	// Assign Ocean and Atmos vectors
	void assign(Teuchos::RCP<Epetra_Vector> vector);
	void assign(std::shared_ptr<std::vector<double> > vector);
	
	// Get the total length of the combined vector
	int length() const;

	// Apply vector update:
	// this = scalarA * A + scalarThis * this
	void update(double scalarA,	SuperVector const &A, double scalarThis);

	// Update element in the vector
	// X[i] = scalarThis*X[i] + scalar
	void updateElement(int index, double scalar, double scalarThis);
	
	// Apply dot product
	double dot(SuperVector const &A) const;

	// Obtain norm of vector
	// mode: 'S': silent
	//       'V': verbose
	//            prints the separate norms of each vector,
	//            returns complete norm
	//            supply message <msg> to improve output
	//       'E': elaborate
	//            prints components of the vectors
	//            supply message <msg> to improve output	
	double norm(char mode = 'S', std::string const msg = "") const;

	// Apply scaling to vectors
	void scale(double scale) const;

	// Get a pointer to the Trilinos/Ocean vector
	Teuchos::RCP<Epetra_Vector> getOceanVector() const;

	// Get a pointer to the STL/Atmos vector
	std::shared_ptr<std::vector<double> > getAtmosVector() const;

	// Randomize vectors
	void random();

	// Zero vectors
	void zero();
	void zeroAtmos();
	void zeroOcean();

	// Remove vectors
	void removeAtmos();
	void removeOcean();
		
	// Print the vectors to various output streams
	void print() const;
	void print(std::string const filename) const;

	// Create a hash of the SuperVector
	std::size_t hash() const;

	// Apply a linear mapping given by a diagonal matrix.
	//   -maps the ocean to the atmosphere and vice versa
	//   -options domain/range: 'O': ocean
	//                          'A': atmosphere
	// The result will be available in *this object.
	void linearTransformation(std::vector<double> const &values,
							  std::vector<int> const &indices,
							  char domain, char range);

	// Apply a linear mapping given by an Epetra_CrsMatrix,
	//  only capable of mapping the ocean to the ocean.
	void linearTransformation(Teuchos::RCP<Epetra_CrsMatrix> mat);

	// Apply a linear mapping given by a CRS matrix in a
	// std::shared_ptr<map<std::string, vector<double> > >
	//  only capable of mapping atmos to atmos.
	void linearTransformation(std::shared_ptr<std::map<std::string,
							  std::vector<double> > > mat);

	void info() const;
	
private:	
	// Perform a few initializations
	// right now it only initializes length_
	void init();
	void assign(SuperVector const &other);

	// Print detailed norm information on all components
	void componentNorms() const;
	
};
#endif
