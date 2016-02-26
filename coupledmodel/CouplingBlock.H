#ifndef COUPLING_BLOCK_H
#define COUPLING_BLOCK_H

#include <vector>
#include <string>
#include <memory> // shared_ptr

#include "SuperVector.H"

#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Import.h>
#include <Teuchos_RCP.hpp>

class CouplingBlock
{
// string defining the coupling type
	std::string coupling_;
	
// Storage for sparse coupling block
	std::vector<double> values_;
	std::vector<int>    row_ind_;
	std::vector<int>    col_ind_;

// Dimensions
	int n_, m_;

	Teuchos::RCP<Epetra_BlockMap> indexMap_;    // Submap based on restriction indices
	Teuchos::RCP<Epetra_Vector> restrVec_;      // Distributed restricted vector
	Teuchos::RCP<Epetra_Import> restrImp_;      // Importer for restricted vector
	Teuchos::RCP<Epetra_MultiVector> gathered_; // Gathered restricted vector

	bool initialized_;

public:
// Default constructor
	CouplingBlock(){}
	
// Constructor
	CouplingBlock(std::string const &coupling,
				  std::vector<double> const &values,
				  std::vector<int> const &row_ind,
				  std::vector<int> const &col_ind);

// Destructor
	~CouplingBlock(){};

// Fill our storage
	void setValues(std::vector<double> const &values);
	void setColInd(std::vector<int> const &col_ind);
	void setRowInd(std::vector<int> const &row_ind);
	
// Apply the coupling block
	void applyMatrix(SuperVector const &v, SuperVector &out);

	void info();

private:
	void initialize();
	void initializeOA(Epetra_BlockMap const &map);
	void initializeAO();
	void applyOA(SuperVector const &v, SuperVector &out);
	void applyAO(SuperVector const &v, SuperVector &out);
};

#endif
