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
#include "Ocean.H"
#include "Continuation.H"
#include "GlobalDefinitions.H"
#include "THCMdefs.H"
#include "Topo.H"

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

		infofile  << "info_" << comm->MyPID() << ".txt";
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

//------------------------------------------------------------------
void printProfile(ProfileType profile)
{
	if (timerStack.empty() == false)
		WARNING("Unequal amount of TIMER_START and TIMER_STOP uses",
				__FILE__, __LINE__);
	
	std::ostringstream profilefile("profile_output");   // setting up a filename
	std::ofstream file(profilefile.str().c_str());      // setup output file

	// Set format flags
	file << std::left;

	// Define line format
#ifndef LINE
# define LINE(s1, s2, s3, s4, s5, s6, s7, s8, s9)						\
	{																	\
		int sp = 3;  int it = 5;  int id = 5;							\
		int db = 12; int st = 45;										\
		file << std::setw(id) << s1	<< std::setw(sp) << s2				\
			 << std::setw(st) << s3 << std::setw(sp) << s4				\
			 << std::setw(db) << s5	<< std::setw(sp) << s6				\
			 << std::setw(it) << s7	<< std::setw(sp) << s8				\
			 << std::setw(db) << s9	<< std::endl;						\
	}
#endif

	// Header
	LINE("", "", "", "", "cumul.", "", "calls", "", "average");
	
	// Display timings of the separate models, summing
	int counter = 0;
	for (auto const &map : profile)
		if (map.first.compare(0,5,"(itr)") != 0)
		{
			counter++;
			std::stringstream s;
			s << " (" << counter << ")";
			LINE(s.str(), "", map.first, ":", map.second[0], "",
				 map.second[1], "", map.second[2]);
		}

	// Newline
	file << std::endl;
	
	// Display iteration information
	for (auto const &map : profile)
		if (map.first.compare(0,5,"(itr)") == 0 )
		{
			counter++;
			std::stringstream s;
			s << " (" << counter << ")";
			LINE(s.str(), "", map.first.substr(5), ":", map.second[0], "",
				 map.second[1], "", map.second[2]);
		}	
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
	initializeEnvironment(argc, argv);
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");

	// Create parallel Ocean 
	RCP<Teuchos::ParameterList> oceanParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());

	Teuchos::RCP<Ocean> ocean = Teuchos::rcp(new Ocean(comm, oceanParams));	

	// Create topography class
	RCP<Teuchos::ParameterList> topoParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("topo_params.xml", topoParams.ptr());

	RCP<Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> > > topo =
		rcp(new Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> >
			(ocean, topoParams));
	
	// Create parameter object for continuation
	RCP<Teuchos::ParameterList> continuationParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("continuation_params.xml",
								continuationParams.ptr());
	
	// Create continuation
	Continuation<RCP<Topo<RCP<Ocean>, RCP<Teuchos::ParameterList> > >,
				 RCP<Teuchos::ParameterList> >
		continuation(topo, continuationParams);
	
	// Run
	int nMasks    = topo->nMasks();
	int startMask = topo->startMaskIdx();
	int status = 0;
	for (int maskIdx = startMask; maskIdx != nMasks-1; maskIdx++)
	{
		topo->setMaskIndex(maskIdx);		
		topo->setPar(0.0);

		TIMER_START("  TOPO:  Predictor I");
		topo->predictor();
		TIMER_STOP ("  TOPO:  Predictor I");

		TIMER_START("  TOPO:  Predictor II");
		continuation.run();
		TIMER_STOP ("  TOPO:  Predictor II");

		topo->preProcess();
		
		TIMER_START("  TOPO:  Corrector");
		status = topo->corrector();
		TIMER_STOP ("  TOPO:  Corrector");
		
		if (status)
		{			
			INFO(" Corrector failed! Abort");
			break;
		}
		
		topo->setPar(1.0);

		if (argc == 1)
			topo->postProcess();
		else
			std::cout << " arguments given, no usage yet" << std::endl;		
		
	}
	
	ocean = Teuchos::null;
	topo  = Teuchos::null;
	
	// print the profile
	if (comm->MyPID() == 0)
		printProfile(profile);
	
	comm->Barrier();
	MPI_Finalize();
}
