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
#include <Epetra_CrsMatrix.h>
#include <Teuchos_RCP.hpp>

#include "TRIOS_Domain.H"

template<typename ModelRow, typename ModelCol>
class CouplingBlock
{
    ModelRow modelRow_;
    ModelCol modelCol_;
    
    Teuchos::RCP<Epetra_CrsMatrix> block_;
    Teuchos::RCP<TRIOS::Domain> modelRowDomain_;
    Teuchos::RCP<TRIOS::Domain> modelColDomain_;

    std::string name_;

    bool initialized_, computed_;

public:

    //------------------------------------------------------------------
    // Default constructor
    CouplingBlock()
        :
        name_("None"),
        initialized_(false),
        computed_(false)
        {}
    
    //------------------------------------------------------------------
    // Constructor
    CouplingBlock(ModelRow modelRow, ModelCol modelCol)
        :
        modelRow_(modelRow),
        modelCol_(modelCol)
        {
            // create name
            std::stringstream ss;
            ss << modelRow->name() << "-" << modelCol->name();
            name_ = ss.str();

            // obtain domains
            modelRowDomain_ = modelRow_->getDomain();
            modelColDomain_ = modelCol_->getDomain();

            // initialize block
            block_ =
                Teuchos::rcp(
                    new Epetra_CrsMatrix(Copy,
                                         *modelRowDomain_->GetSolveMap(),
                                         *modelColDomain_->GetColMap(), 0) );

            computed_    = false;
            initialized_ = true;

            // compute the block
            computeBlock();
        }

    //------------------------------------------------------------------
    void computeBlock()
        {
            INFO("CouplingBlock: computing " << name_);
            assert(initialized_);

            // local size
            int numMyElements = block_->RowMap().NumMyElements();
            int numGlElements = block_->RowMap().NumGlobalElements();

            // obtain 0-based global CRS matrix from modelRow
            std::shared_ptr<Utils::CRSMat> blockCRS =
                modelRow_->getBlock(modelCol_);
            
            // inspect CRS struct
            if (blockCRS->beg.empty())
            {
                WARNING(name_ << ": Empty CRS struct, not computing coupling block"
                        << "   flags: " << computed_ << " "
                        << initialized_,  __FILE__, __LINE__);
                return;
            }

            // construct CRS matrix

            // max nonzeros per row using dof in models
            // and max neighbours=27
            const int maxnnz =
                (modelRow_->dof() + modelCol_->dof()) * 27 + 1;
            
            // indices array
            std::vector<int> indices(maxnnz, 0);

            // values array
            std::vector<double> values(maxnnz, 0.0);
            

            // global case
            if (numGlElements == (int) blockCRS->beg.size() - 1)
            {
                TIMER_START("CouplingBlock: compute block");
                
                int gRow, index, numentries;
                
                // fill Epetra CRS from global CRSMat
                for (int i = 0; i < numMyElements; ++i)
                {                    
                    gRow       = block_->RowMap().GID(i);
                    index      = blockCRS->beg[gRow];
                    numentries = blockCRS->beg[gRow+1] - index;

                    // if we encounter a dense row (probably an integral equation)
                    if (numentries > maxnnz)
                    {
                        indices = std::vector<int>(numentries, 0);
                        values  = std::vector<double>(numentries, 0);
                    }
                    
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
                                                        &values[0], &indices[0]);
                    }
                    else
                    {
                        ierr =
                            block_->InsertGlobalValues(gRow, numentries,
                                                       &values[0], &indices[0]);
                    }
                    
                    if (ierr != 0)
                    {
                        INFO (name_ << ": Error in InsertGlobalValues: " << ierr);
                        INFO (" Filled ? " << block_->Filled());
                        std::cout << name_ << ": Error in Insert/ReplaceGlobalValues: "
                                  << ierr << std::endl;
                        std::cout << "Filled = " << block_->Filled() << std::endl;
                        std::cout << "  GRID = " << gRow << std::endl;
                        std::cout << "  LRID = " << block_->LRID(gRow) << std::endl;
                        std::cout << " graph inds in LRID:   " 
                              << block_->Graph().NumMyIndices(block_->LRID(gRow)) << std::endl;

                        std::cout << "indices : ";
                        for (int ii = 0; ii != numentries; ++ii)
                        {
                            std::cout << indices[ii] << " ";
                        }
                        std::cout << std::endl;

                        ERROR("Error in InsertGlobalValues", __FILE__, __LINE__);
                    }
                }
                
                TIMER_STOP("CouplingBlock: compute block");
            }
            else
            {
                WARNING(name_ << ": local CRS case not implemented! Continue with empty coupling block.",
                        __FILE__, __LINE__);
                return;
            }

            // Finalize
            CHECK_ZERO(block_->FillComplete(
                           *modelColDomain_->GetSolveMap(),
                           *modelRowDomain_->GetSolveMap()));

            computed_ = true;

        }

    //------------------------------------------------------------------
    // Destructor
    ~CouplingBlock() {}

    //------------------------------------------------------------------
    void applyMatrix(Epetra_MultiVector const &in, Epetra_MultiVector &out)
        {
            if (!computed_ || !initialized_)
            {
                WARNING(name_ << ": trying to apply empty coupling block", __FILE__, __LINE__);
                out.PutScalar(0.0);
                return;
            }
            
            TIMER_START("CouplingBlock: apply matrix...");
            assert(block_->DomainMap().SameAs(in.Map()));
            assert(block_->RangeMap().SameAs(out.Map()));
            CHECK_ZERO(block_->Apply(in, out));
            TIMER_STOP("CouplingBlock: apply matrix...");
        }

    //------------------------------------------------------------------
    // Get RCP to block
    Teuchos::RCP<Epetra_CrsMatrix> getBlock()
        {
            if (!computed_ || !initialized_)
            {
                WARNING(name_ << ": returning empty coupling block", __FILE__, __LINE__);
            }

            return block_;
        }

    std::string const name() { return name_; }
        
};

#endif
