#ifndef TESTDEFINITIONS_H
#define TESTDEFINITIONS_H

//=============================================================================
// General definitions for tests
//=============================================================================

#include <Epetra_config.h>

#  ifdef HAVE_MPI
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif // HAVE_MPI

#include <EpetraExt_HDF5.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <ios>
#include <iomanip>
#include <memory>
#include <vector>
#include <array>
#include <stack>
#include <string>
#include <math.h> // sqrt
#include <time.h> // for seeding
#include <sstream>

#include "../main/RunDefinitions.H"

#include "THCMdefs.H"
#include "ComplexVector.H"
#include "JDQZInterface.H"
#include "MyParameterList.H"
#include "NumericalJacobian.H"

#include "jdqz.H"

#include "gtest/gtest.h" // google test

// We need this to compute norms
extern "C" double ddot_(int *N, double *X, int *INCX, double *Y, int *INCY);

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;
using std::shared_ptr;

//------------------------------------------------------------------
double norm(Teuchos::RCP<Epetra_Vector> vec)
{
    double nrm;
    vec->Norm2(&nrm);
    return nrm;
}

//------------------------------------------------------------------
double norm(std::shared_ptr<std::vector<double> > vec)
{
    int dim = (int) vec->size();
    int incX = 1;
    int incY = 1;
    double dot;
    dot = ddot_(&dim, &(*vec)[0], &incX, &(*vec)[0], &incY);
    return sqrt(dot);
}

//------------------------------------------------------------------
std::shared_ptr<std::vector<double> >
getGatheredVector(Teuchos::RCP<Epetra_Vector> vec)
{
    // get size
    int dim = (int) vec->GlobalLength();

    // Create gather scheme: map and importer
    Teuchos::RCP<Epetra_BlockMap> gmap =
        Utils::AllGather(vec->Map());
    Teuchos::RCP<Epetra_Import> imp =
        Teuchos::rcp(new Epetra_Import(*gmap, vec->Map()));

    // Create vector to hold the gathered state
    Teuchos::RCP<Epetra_Vector> gathered =
        Teuchos::rcp(new Epetra_Vector(*gmap));

    // Gather the state into <gathered>
    gathered->Import(*vec, *imp, Insert);

    // Create full state for serial atmosphere
    std::shared_ptr< std::vector<double> > gvec =
        std::make_shared<std::vector<double> >(dim, 0.0);

    // Extract gvec from gathered
    gathered->ExtractCopy(&(*gvec)[0], dim);

    return gvec;
}

//------------------------------------------------------------------
template<typename ModelPtr, typename VectorPtr, typename CCS>
void testEntries(ModelPtr model, CCS &ccs, VectorPtr x)
{
    x->PutScalar(0.0);
    int myLength = x->MyLength();
    EXPECT_EQ(myLength, ccs.beg.size() - 1);
            
    typename VectorPtr::element_type e_i(*x);
    typename VectorPtr::element_type e_j(*x);
    typename VectorPtr::element_type tmp(*x);
            
    // assuming we run on a single core
    int    ico  = 0;
    double Jij  = 0.0;

    double error, diff;
    for (int jco = 0; jco != myLength; ++jco)
    {
        e_j[jco] = 1;
        tmp.PutScalar(0.0);

        // tmp = jth column 
        model->applyMatrix(e_j, tmp);
        for (int i = ccs.beg[jco]; i != ccs.beg[jco+1]; ++i)
        {
            ico = ccs.ico[i];
            e_i[ico] = 1;

            // get Jij component
            e_i.Dot(tmp, &Jij);
            e_i[ico] = 0;

            // difference between analytical and numerical value
            diff  = std::abs(Jij - ccs.co[i]);                    
                    
            error = std::max(1e-2 * std::abs(Jij), 1e-10);
            EXPECT_NEAR(diff, 0.0, error);
                                                  
            if (diff > error) // some additional outputting
            {
                std::cout << "(" << ico << "," << jco << ") ";
                std::cout << "nuJac: " << ccs.co[i] << " " << "anJac: " << Jij << std::endl;
                int mdl,i,j,k,xx;
                        
                model->gid2coord(ico,mdl,i,j,k,xx);
                std::cout << "\nrow: mdl=" << mdl << " i=" << i << " j=" << j
                          << " k=" << k << " xx=" << xx << std::endl;
                        
                model->gid2coord(jco,mdl,i,j,k,xx);
                std::cout << "col: mdl=" << mdl << " i=" << i << " j=" << j
                          << " k=" << k << " xx=" << xx << std::endl << std::endl;
            }
        }
        e_j[jco] = 0; // reset value
    }            
}

#endif
