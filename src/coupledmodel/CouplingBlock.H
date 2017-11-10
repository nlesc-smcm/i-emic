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

    //------------------------------------------------------------------
    // Default constructor
    CouplingBlock() {}

    //------------------------------------------------------------------
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
                                         *modelColDomain_->GetColMap(), 0) );

            // initial compute
            computeBlock();
        }

    //------------------------------------------------------------------
    void computeBlock()
        {
            TIMER_START("CouplingBlock: compute block");

            // obtain 0-based CRS matrix from modelRow
            std::shared_ptr<Utils::CRSMat> blockCRS =
                modelRow_->getBlock(modelCol_);

            // construct CRS matrix

            // max nonzeros per row using dof in models
            // and max neighbours=27
            const int maxnnz =
                (modelRow_->dof() + modelCol_->dof()) * 27 + 1;
            
            // indices array
            int indices[maxnnz];

            // values array
            double values[maxnnz];

            // local size
            int numMyElements = block_->RowMap().NumMyElements();
            int numGlElements = block_->RowMap().NumGlobalElements();

            // global case
            if (numGlElements == (int) blockCRS->beg.size() - 1)
            {
                int gRow, index, numentries;

                // fill Epetra CRS from global CRSMat
                for (int i = 0; i < numMyElements; ++i)
                {
                    gRow       = block_->RowMap().GID(i);
                    index      = blockCRS->beg[gRow];
                    numentries = blockCRS->beg[gRow+1] - index;
                    for (int j = 0; j < numentries; ++j)
                    {
                        indices[j] = blockCRS->jco[index+j];
                        values[j]  = blockCRS->co[index+j];
                    }

                    int ierr;
                    if (block_->Filled())
                    {
                        ierr =
                            block_->ReplaceGlobalValues(gRow, numentries,
                                                        values, indices);
                    }
                    else
                    {
                        ierr =
                            block_->InsertGlobalValues(gRow, numentries,
                                                       values, indices);
                    }

                    if (ierr != 0)
                        ERROR("Error in InsertGlobalValues", __FILE__, __LINE__);
                }
            }
            else
                ERROR("Local CRS case not implemented!", __FILE__, __LINE__);

            // Finalize
            CHECK_ZERO(block_->FillComplete(
                           *modelColDomain_->GetSolveMap(),
                           *modelRowDomain_->GetSolveMap()));

            TIMER_STOP("CouplingBlock: compute block");
        }

    //------------------------------------------------------------------
    // Destructor
    ~CouplingBlock()
        {
            INFO("CouplingBlock destructor");
        }

    //------------------------------------------------------------------
    void applyMatrix(Epetra_MultiVector const &in, Epetra_MultiVector &out)
        {
            TIMER_START("CouplingBlock: apply matrix...");
            assert(block_->DomainMap().SameAs(in.Map()));
            assert(block_->RangeMap().SameAs(out.Map()));
            CHECK_ZERO(block_->Apply(in, out));
            TIMER_STOP("CouplingBlock: apply matrix...");
        }

    //------------------------------------------------------------------
    // Get RCP to block
    Teuchos::RCP<Epetra_CrsMatrix> getBlock() { return block_; } 
        
};

#endif
