//==============================================================================
// Test routine
//==============================================================================
// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>

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

#include <Epetra_Vector.h>

#include "GlobalDefinitions.H"
#include "OceanTheta.H"
#include "ThetaStepper.H"

//==============================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//==============================================================================
// A few declarations
RCP<std::ostream> outFile;
RCP<std::ostream> outputFiles(RCP<Epetra_Comm> Comm);
RCP<Epetra_Comm>  initializeEnvironment(int argc, char **argv);

//==============================================================================
int main(int argc, char **argv)
{
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);
	
 	// Create timestepping (theta method) ocean model OceanTheta
	//  (based on THCM):
	RCP<OceanTheta> ocean = rcp(new OceanTheta(Comm));

	// Create timestepping object ThetaStepper using an OceanTheta model
	//  and an Epetra_Vector
	ThetaStepper<RCP<OceanTheta>, RCP<Epetra_Vector> >
		thetaStepper(ocean, ocean->GetState());
	
	thetaStepper.Run();	
	ocean->DumpState();
	
    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
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
