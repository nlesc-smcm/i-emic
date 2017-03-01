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

#include "SuperVector.H"
#include "CoupledModel.H"
#include "Continuation.H"
#include "Atmosphere.H"
#include "AtmospherePar.H"
#include "GlobalDefinitions.H"
#include "THCMdefs.H"
#include "JDQZInterface.H"
#include "MyParameterList.H"

#include "jdqz.H"

#include "gtest/gtest.h" // google test

// We need this to compute norms
extern "C" double ddot_(int *N, double *X, int *INCX, double *Y, int *INCY);

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;
using std::shared_ptr;

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
	if (comm->MyPID() < 4)
	{
		std::ostringstream infofile;     // setting up a filename

		infofile  << "info_"    << comm->MyPID()   << ".txt";
		
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

	// seed random number generator
	srand(time(NULL));
}

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
