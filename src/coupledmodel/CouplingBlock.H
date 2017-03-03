#ifndef COUPLING_BLOCK_H
#define COUPLING_BLOCK_H

#include <vector>
#include <string>
#include <memory> // shared_ptr

#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Import.h>
#include <Teuchos_RCP.hpp>

template<typename ModelRow, typename ModelCol>
class CouplingBlock
{	
	ModelRow modelRow_;
	ModelCol modelCol_;

	Teuchos::RCP<Epetra_CrsMatrix> block_;

	Teuchos::RCP<TRIOS::Domain> modelRowDomain_;
	Teuchos::RCP<TRIOS::Domain> modelColDomain_;
	
public:
	// Default constructor
	CouplingBlock() {}
	
	// Constructor  
	CouplingBlock(ModelRow modelRow, ModelCol modelCol)
		:
		modelRow_(modelRow),
		modelCol_(modelCol)
		{
			modelRowDomain_ = modelRow_->getDomain();
			modelColDomain_ = modelCol_->getDomain();
			
			block_ =
				Teuchos::rcp(
					new Epetra_CrsMatrix(Copy,
										 *modelRowDomain_->GetSolveMap(),
										 *modelColDomain_->GetColMap(), 0));

			
			std::shared_ptr<Utils::CRSMat> blockCRS = modelCol_->getBlock(modelRow_);
			
		}

// Destructor
	~CouplingBlock()
		{
			INFO("CouplingBlock destructor");
		};

private:

	
};

#endif
