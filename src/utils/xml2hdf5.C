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

#include "Ocean.H"
#include "SuperVector.H"
#include "GlobalDefinitions.H"
#include "Continuation.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
// A few declarations (see GlobalDefinitions.H)
// --> put them in a namespace
RCP<std::ostream> outFile;      // output file
ProfileType       profile;      // profile
std::stack<Timer> timerStack;   // timing stack

//------------------------------------------------------------------
RCP<std::ostream> outputFiles(RCP<Epetra_Comm> Comm);
RCP<Epetra_Comm>  initializeEnvironment(int argc, char **argv);
void testOcean(RCP<Epetra_Comm> Comm);
//------------------------------------------------------------------
int main(int argc, char **argv)
{
	if (argc != 2)
		throw std::runtime_error("ERROR: Supply xml file");		
	
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);

	// Check if outFile is specified
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");


	// Create parameter object for Ocean
	RCP<Teuchos::ParameterList> oceanParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());

	// Create parallelized Ocean object
	RCP<Ocean> ocean = Teuchos::rcp(new Ocean(Comm, oceanParams));

	INFO("*****************************************************");
	INFO("    Replacing the ocean hdf5 parameters with ");
	INFO("    the parameters supplied in " << argv[1]   );
	INFO("*****************************************************");	
	
	ocean->loadState();

	INFO("*****************************************************");
	INFO("*****************************************************");	

	// Create additional parameter object for Ocean
	RCP<Teuchos::ParameterList> extraParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile(argv[1], extraParams.ptr());

	ocean->setParameters(extraParams);
	
	ocean->saveState();

	INFO("*****************************************************");
	INFO("*****************************************************");	

    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	MPI_Finalize();	
}

//------------------------------------------------------------------
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

//------------------------------------------------------------------
Teuchos::RCP<std::ostream> outputFiles(Teuchos::RCP<Epetra_Comm> Comm)
{
	// Setup output files "fname_#.txt" for P==0 && P==1, other processes
	// will get a blackholestream.
	Teuchos::RCP<std::ostream> outFile;
	if (Comm->MyPID() < 2)
	{
		std::ostringstream infofile;     // setting up a filename

		infofile    << "info_"    << Comm->MyPID()   << ".txt";

		std::cout << "info for CPU" << Comm->MyPID() << " is written to "
				  << infofile.str().c_str() << std::endl;

		outFile =
			Teuchos::rcp(new std::ofstream(infofile.str().c_str()));
	}
	else
	{
		outFile =
			Teuchos::rcp(new Teuchos::oblackholestream());
	}
	return outFile;
}
