#include "Rearranger.H"
#include "GlobalDefinitions.H"
#include <EpetraExt_RowMatrixOut.h>


void Rearranger::setMatrix(RCP<Epetra_CrsMatrix> matrix)
{
	INFO("Entering setMatrix()");
	matrix_ = rcp(&(*matrix), false);
	matrixFilled_ = true;
	EpetraExt::RowMatrixToMatrixMarketFile("matrix_.txt", *matrix_);
	INFO("Leaving setMatrix() ");
}

void Rearranger::buildOrdering()
{
	INFO("Entering buildOrdering()");
	if (!matrixFilled_)
	{
		INFO(" Matrix not set, returning");
		return;
	}
	int dof = 6; // Number of unknowns per grid cell (degrees of freedom)
	int N   = matrix_->NumMyRows() / dof; // Local number of rows
	INFO(" N = " << N);
	
	// We want to have an ordering with ((u,v),w,p,(T,S)), so a 4x4 block
	// structure. For this we need vectors of size 2*N, N, N and 2*N respectively:
	ordering_.push_back(std::vector<int>(2*N)); // Will give the (u,v) indices.
	ordering_.push_back(std::vector<int>(N));   // Will give the (w) indices.
	ordering_.push_back(std::vector<int>(N));   // Will give the (p) indices.
	ordering_.push_back(std::vector<int>(2*N)); // Will give the (T,S) indices.

	// Obtain the RowMap from matrix_
    Epetra_Map const &map = matrix_->RowMatrixRowMap();

	// Fill ordering_ with the GID's of ((u,v),w,p,(T,S))
	for (int i = 0; i != N; ++i)
    {
        ordering_[0][i*2]    =  map.GID(i*6);   // u
		ordering_[0][i*2+1]  =  map.GID(i*6+1); // v
		ordering_[1][i]      =  map.GID(i*6+2); // w
		ordering_[2][i]      =  map.GID(i*6+3); // p
		ordering_[3][i*2]    =  map.GID(i*6+4); // T
		ordering_[3][i*2+1]  =  map.GID(i*6+5); // S
    }
	orderingFilled_ = true;
	INFO("Leaving buildOrdering()");
}

void Rearranger::setBlockOperator()
{
	INFO("Entering setBlockOperator()");
	if (!orderingFilled_)
	{
		INFO(" Ordering not build, returning");
		return;
	}
	// Set the BlockedEpetraOperator
	blockOperator_ = rcp(new Teko::Epetra::BlockedEpetraOperator(ordering_,matrix_));
	blockOperatorFilled_ = true;
	INFO("Leaving setBlockOperator()");
}

void Rearranger::fillBlocks()
{
	INFO("Entering fillBlocks()");
	if (!blockOperatorFilled_)
	{
		INFO(" blockOperator not build, returning");
		return;
	}
	int i,j; 	// block indices
	std::string fname;
	for (int idx = 0; idx != numNonzBlocks_; ++idx)
	{
		i = blockLocations_[idx][0];
		j = blockLocations_[idx][1];
		INFO("Creating " << keys_[idx]);
		blocks_[keys_[idx]] = rcp(new Epetra_CrsMatrix(
									  *(Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>
										(blockOperator_->GetBlock(i,j))
										  )));
		fname = keys_[idx];
		fname += ".txt";
		EpetraExt::RowMatrixToMatrixMarketFile(fname.c_str(), *(blocks_[keys_[idx]]));
	}
	INFO("Leaving fillBlocks()");
}
