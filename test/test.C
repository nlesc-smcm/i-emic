//==============================================================================
// Test routine
//==============================================================================
// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>
#include <iomanip>

//==============================================================================
#  ifdef HAVE_MPI
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif // HAVE_MPI

//==============================================================================
#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "Vector.H"
#include "OceanCont.H"
#include "ThetaStepper.H"
#include "Continuation.H"
#include "GlobalDefinitions.H"

//==============================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//==============================================================================
// --> ugly!
// A few global things, see ocean/GlobalDefinitions.H 
std::map<std::string, double> profile;  // profiler
RCP<std::ostream> outFile;              // output file

//==============================================================================
RCP<std::ostream> outputFiles(RCP<Epetra_Comm> Comm);
RCP<Epetra_Comm>  initializeEnvironment(int argc, char **argv);
void printProfile(std::map<std::string, double> profile);

//==============================================================================
int main(int argc, char **argv)
{
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);
	
 	// Create ocean model OceanCont
	//  (based on THCM):
	RCP<OceanCont> ocean = rcp(new OceanCont(Comm));

	// Create continuation
	Continuation<RCP<OceanCont>, RCP<Vector> >
		continuation(ocean);

	continuation.Run();
	
	ocean->DumpState();

	printProfile(profile);
    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	Comm->Barrier();
	MPI_Finalize();	
}

//==============================================================================
// Auxiliary stuff
//==============================================================================
RCP<Epetra_Comm> initializeEnvironment(int argc, char **argv)
{
	
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	RCP<Epetra_MpiComm> Comm =
		rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	RCP<Epetra_SerialComm> Comm =
		rcp(new Epetra_SerialComm());
#endif
	// Specify output files
	outFile = outputFiles(Comm);
	return Comm;
}

//==============================================================================
Teuchos::RCP<std::ostream> outputFiles(Teuchos::RCP<Epetra_Comm> Comm)
{
	// Setup output files "fname_#.txt" for P==0 && P==1, other processes
	// will get a blackholestream.
	Teuchos::RCP<std::ostream> outFile;
	if (Comm->MyPID() == 0)
	{
		std::ostringstream infofile;  // setting up a filename
		infofile << "info_" << Comm->MyPID() << ".txt";
		std::cout << "info for P=" << Comm->MyPID() << " is written to "
				  << infofile.str().c_str() << std::endl;
		outFile = 
			Teuchos::rcp(new std::ofstream(infofile.str().c_str()));
	}
	else
		outFile = 
			Teuchos::rcp(new Teuchos::oblackholestream());
	return outFile;
}

//===========================================================================
void printProfile(std::map<std::string, double> profile)
{
	INFO("==================================================================")
	INFO(" Profile:");
	for (std::map<std::string, double>::iterator it = profile.begin();
		 it != profile.end(); ++it)
	{
		INFO( std::setw(40) << it->first << " -> " <<
			  std::setw(12) << std::setprecision(8) << it->second);
	}
	INFO("==================================================================")
}
