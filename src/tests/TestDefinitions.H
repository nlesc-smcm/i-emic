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

#include "SuperVector.H"
#include "CoupledModel.H"
#include "Continuation.H"
#include "Atmosphere.H"
#include "GlobalDefinitions.H"
#include "THCMdefs.H"
#include "JDQZInterface.H"
#include "MyParameterList.H"

#include "jdqz.H"

#include "gtest/gtest.h" // google test

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
// A few globals (see GlobalDefinitions.H)
//------------------------------------------------------------------
RCP<std::ostream> outFile;      // output file
ProfileType       profile;      // profile
std::stack<Timer> timerStack;   // timing stack
RCP<Epetra_Comm>  comm;         // communicator object

//------------------------------------------------------------------
RCP<std::ostream> outputFiles()
{
	Teuchos::RCP<std::ostream> outFile;
	if (comm->MyPID() < 1)
	{
		std::ostringstream infofile;     // setting up a filename

		infofile    << "info_"    << comm->MyPID()   << ".txt";

		std::cout << "info for CPU" << comm->MyPID() << " is written to "
				  << infofile.str().c_str() << std::endl;

		outFile = Teuchos::rcp(new std::ofstream(infofile.str().c_str()));
	}
	else
	{
		outFile = Teuchos::rcp(new Teuchos::oblackholestream());
	}
	return outFile;
}

//------------------------------------------------------------------
void initializeEnvironment(int argc, char **argv)
{
#ifdef HAVE_MPI           // Initialize communicator
	MPI_Init(&argc, &argv);
	comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	comm = rcp(new Epetra_SerialComm());
#endif
	outFile = outputFiles(); 	// Initialize output files
}
