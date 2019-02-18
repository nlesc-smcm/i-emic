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
#include "jdqz.H"

//-----------------------------------------------------------------------------

using Teuchos::RCP;
using Teuchos::rcp;

RCP<std::ostream> outFile;      // output file
RCP<std::ostream> cdataFile;    // cdata file
RCP<std::ostream> tdataFile;    // tdata file

//-----------------------------------------------------------------------------
void  outputFiles(RCP<Epetra_Comm> Comm)
{
    // Setup output files "fname_#.txt" for P==0 && P==1, other processes
    // will get a blackholestream.
    
  if (Comm->MyPID() < 10)
    {
        std::ostringstream infofile;   // setting up a filename

        infofile    << "info_"    << Comm->MyPID()   << ".txt";

        std::cout << "info for CPU" << Comm->MyPID() << " is written to "
                  << infofile.str().c_str() << std::endl;

        outFile =
            Teuchos::rcp(new std::ofstream(infofile.str().c_str()));
    }
    else
        outFile = Teuchos::rcp(new Teuchos::oblackholestream());

    // Write continuation data on proc 0
    if (Comm->MyPID() == 0)
    {
        std::ostringstream cdatafile;  // setting up a filename

        cdatafile << "cdata.txt";

        std::cout << "continuation data is written to " <<  cdatafile.str().c_str()
                  << " by CPU " << Comm->MyPID() << std::endl;            
        
        cdataFile =
            Teuchos::rcp(new std::ofstream(cdatafile.str().c_str()));
    }
    else
        cdataFile = Teuchos::rcp(new Teuchos::oblackholestream());
        
    // Write timestepping data on proc 0
    if (Comm->MyPID() == 0)
    {
        std::ostringstream tdatafile;  // setting up a filename

        tdatafile << "tdata.txt";

        std::cout << "continuation data is written to " <<  tdatafile.str().c_str()
                  << " by CPU " << Comm->MyPID() << std::endl;            
        
        tdataFile =
            Teuchos::rcp(new std::ofstream(tdatafile.str().c_str()));
    }
    else
        tdataFile = Teuchos::rcp(new Teuchos::oblackholestream());
        
    
}

//-----------------------------------------------------------------------------
RCP<Epetra_Comm>  initializeEnvironment(int argc, char **argv)
{
    // Setup MPI communicator

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    RCP<Epetra_MpiComm> Comm =
        rcp(new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
    RCP<Epetra_SerialComm> Comm =
        rcp(new Epetra_SerialComm() );
#endif
    
    // Specify output files
    outputFiles(Comm);
    return Comm;
}

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
