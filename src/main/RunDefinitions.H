#ifndef RUNDEFINITIONS_H
#define RUNDEFINITIONS_H

//=============================================================================
// General definitions for runs
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

#include "CoupledModel.H"
#include "Atmosphere.H"
#include "Ocean.H"
#include "SeaIce.H"
#include "Continuation.H"
#include "GlobalDefinitions.H"
#include "Topo.H"
#include "Utils.H"

#include "ComplexVector.H"
#include "JDQZInterface.H"
#include "jdqz.hpp"

//-----------------------------------------------------------------------------

using Teuchos::RCP;
using Teuchos::rcp;

#endif
