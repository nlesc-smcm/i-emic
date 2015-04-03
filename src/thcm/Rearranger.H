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

using Teuchos::RCP;
using Teuchos::rcp;

//! This class rearranges the matrix from THCM into a better conditioned form.
class Rearranger
{
	// The matrix subjected to the rearrangement
	RCP<Epetra_CrsMatrix> matrix_;
	
    // Nested vector describing the ordering of the unknowns (u,v,w,p,T,S)
	// It is a vector of vectors describing how the matrix should be blocked.
	std::vector< std::vector<int> > ordering_;

	// A BlockedEpetraoperator should be the main datamember of this class. It takes
	// care of creating the block reordering.
	RCP<Teko::Epetra::BlockedEpetraOperator> blockOperator_;

	// Number of blocks in rows and columns respectively;
	int numRowBlocks_;
	int numColBlocks_;

	// 10 nonzero matrix blocks using a notation similar to [A. de Niet et al. 2007]
	RCP<Epetra_CrsMatrix> A_uv_, E_uv_, G_uv_, D_uv_, B_uv_,
		G_w_, D_w_, B_w_, B_TS_, A_TS_;

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
			INFO("Rearranger constructor called");			
		}
	
	~Rearranger(){ INFO("Rearranger destructor called"); }
	
	void setMatrix(RCP<Epetra_CrsMatrix> matrix);
	void buildOrdering();
	void setBlockOperator();
	void fillBlocks();
	
};
