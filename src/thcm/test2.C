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
#include <Teuchos_Version.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#include <sstream>
#include "GlobalDefinitions.H"

// for gethostname in pardebug
#include <unistd.h>

using Teuchos::RCP;
using Teuchos::rcp;


Teuchos::RCP<std::ostream> outputFiles(Teuchos::RCP<Epetra_Comm> Comm);
void parDebug();
RCP<std::ostream> outFile;


int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	RCP<Epetra_MpiComm> Comm =
		rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	RCP<Epetra_SerialComm> Comm =
		rcp(new Epetra_SerialComm());
#endif
	// -------------------------------------------------------------------
	// Setup stream for INFO, DEBUG, ERROR and WARNING Macros
	// -------------------------------------------------------------------
	outFile = outputFiles(Comm);

	long long int numGlobalElemenst = 64;
	long long int indexBase = 0;
	Epetra_Map map(numGlobalElemenst, indexBase, *Comm);
	RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(Copy, map, 8));
	RCP<Epetra_CrsMatrix> B = rcp(new Epetra_CrsMatrix(Copy, map, 8));
	int numMyElmntsA = A->map().NumMyElements();
	DEBVAR(A->RowMap());
	DEBVAR(A->ColMap());
	DEBVAR(A->DomainMap());
	DEBVAR(A->RangeMap());
		
	//------------------------------------------------------------------
	// Finalize MPI
	//------------------------------------------------------------------
	MPI_Finalize();
	return 0;
}

Teuchos::RCP<std::ostream> outputFiles(Teuchos::RCP<Epetra_Comm> Comm)
{
	// Setup output files "fname_#.txt" for P==0 && P==1, other processes
	// will get a blackholestream.
	Teuchos::RCP<std::ostream> outFile;
	if (Comm->MyPID() < 2)
	{
		std::ostringstream outfile;  // setting up a filename
		outfile << "out_" << Comm->MyPID() << ".txt";
		std::cout << "Output for Process " << Comm->MyPID() << " is written to "
				  << outfile.str().c_str() << std::endl;
		outFile = Teuchos::rcp(new std::ofstream(outfile.str().c_str()));
	}
	else
	{
		std::cout << "Output for Process " << Comm->MyPID()
				  << " is neglected" << std::endl;
		outFile = Teuchos::rcp(new Teuchos::oblackholestream());
	}
	return outFile;
}
void parDebug()
{
	int i = 0;
	char hostname[256];
	gethostname(hostname, sizeof(hostname));
	printf("PID %d on %s ready for attach\n", getpid(), hostname);
	fflush(stdout);
	while (0 == i)
		sleep(5);
}
