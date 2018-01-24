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

#endif
