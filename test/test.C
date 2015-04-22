// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>

#  ifdef HAVE_MPI
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif // HAVE_MPI

#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Epetra_Vector.h>

#include "GlobalDefinitions.H"
#include "Ocean.H"
#include "Newton.H"
#include "Singleton.H"

using Teuchos::RCP;
using Teuchos::rcp;

RCP<std::ostream> outFile;
RCP<std::ostream> outputFiles(RCP<Epetra_Comm> Comm);

int main(int argc, char **argv)
{
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	RCP<Epetra_MpiComm> Comm =
		rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	RCP<Epetra_SerialComm> Comm =
		rcp(new Epetra_SerialComm());
#endif
	// Specify output streams
	outFile = outputFiles(Comm);
	// Initialize ocean model (THCM)
	RCP<Ocean> ocean = rcp(new Ocean(Comm));
	Newton<RCP<Ocean>, RCP<Epetra_Vector> >
		newtonSolver(ocean);
	
	newtonSolver.run();	
	ocean->dumpState();
	
    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	MPI_Finalize();	
}

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
