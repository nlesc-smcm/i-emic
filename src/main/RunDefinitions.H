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

#include "ComplexVector.H"
#include "JDQZInterface.H"
#include "jdqz.hpp"

//-----------------------------------------------------------------------------

using Teuchos::RCP;
using Teuchos::rcp;

//-----------------------------------------------------------------------------
Teuchos::RCP<Teuchos::ParameterList> obtainParams(std::string const &str,
                                                  std::string const &name)
{
    RCP<Teuchos::ParameterList> pars = rcp(new Teuchos::ParameterList);
    
    std::ifstream file(str);
    if (!file)
    {
        WARNING(str << ", " << name <<
                " not found, continuing with defaults at your own risk!",
                __FILE__, __LINE__);
    }
    else if (str == "dummy")
    {
        INFO("Continuing with dummy " << name);
    }
    else
    {
        updateParametersFromXmlFile(str.c_str(), pars.ptr());
    }    
    pars->setName(name.c_str());
    return pars;
}

#endif
