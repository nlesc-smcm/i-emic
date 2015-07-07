/**********************************************************************
 * Copyright by T.E. Mulder, Utrecht University 2014/15.              *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: t.e.mulder@uu.nl                                          *
 **********************************************************************/
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_RCP.hpp"
#include "Teko_BlockedEpetraOperator.hpp"

#include "GlobalDefinitions.H"

#include <vector>
#include <array>
#include <unordered_map>

using Teuchos::RCP;
using Teuchos::rcp;

//! This class rearranges the matrix from THCM into a better conditioned form.
class Rearranger
{
	// The matrix subjected to the rearrangement
	RCP<Epetra_CrsMatrix> matrix_;
	
    // Nested vector describing the ordering of the unknowns (u,v,w,p,T,S)
	// It is a vector of vectors of ints describing how the matrix should be blocked.
	std::vector< std::vector<int> > ordering_;

	// A BlockedEpetraoperator should be the main datamember of this class. It takes
	// care of creating the block reordering.
	RCP<Teko::Epetra::BlockedEpetraOperator> blockOperator_;

	// Number of nonzero blocks;
	int numNonzBlocks_ = 10;
	
	// An unordered_map with RCP's pointing to the 10 matrix blocks.
	// The keys for the matrix blocks are given below.
	std::unordered_map< std::string, RCP<Epetra_CrsMatrix> > blocks_;

	// vector of strings used for identifying contents of the blocks_ vector.
    // Notation similar to [A. de Niet et al. 2007], see constructor.
	//** This initialization may not work on other compilers (std=c++11) **//
	std::vector< std::string > keys_  = {"A_uv", "E_uv", "G_uv", "G_w", "B_ST",
										 "D_uv", "D_w", "B_uv", "B_w", "A_ST"};

	// Indices giving the locations of the nonzero blocks
	int blockLocations_[10][2]  = { {0,0}, {0,1}, {0,2},
									{1,2}, {1,3},
									{2,0}, {2,1},
									{3,0}, {3,1}, {3,3,} };
	
	// The rebuild and reordered system
	RCP<Epetra_CrsMatrix> finalMatrix_;

	// some flags controlling the status
	bool matrixFilled_;
	bool orderingFilled_;
	bool blockOperatorFilled_;

public:
	Rearranger()
		:
		matrixFilled_(false),
		orderingFilled_(false),
		blockOperatorFilled_(false)
		{
			INFO("Entering Rearranger constructor");

			INFO("Leaving Rearranger constructor");
		}
	
	~Rearranger(){ INFO("Rearranger destructor called"); }
	
	void setMatrix(RCP<Epetra_CrsMatrix> matrix);
	void buildOrdering();
	void setBlockOperator();
	void fillBlocks();
	void test();

	//! find dummy rows in a subset of rows of the matrix A.
    
	/*! Given a matrix A and a (sub-)map M, this function
	    sets the flag array is_dummy to true if the corresponding
	    local node is a dummy row, i.e. if it has only one entry
	    and that entry is 1.0 on the diagonal of the matrix.
	    Memory for the is_dummy array is allocated in this function.
          
		returns the number of dummies found.
	*/
	// Copyright Jonas Thies, Univ. of Groningen 2006/7/8.          
   	int detectDummies(const Epetra_CrsMatrix& A,
					  const Epetra_Map& M, bool *is_dummy) const;
	
};
