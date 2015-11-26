//=======================================================================
// Main continuation of the coupled model 
//=======================================================================
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
#include "Atmosphere.H"
#include "ThetaStepper.H"
#include "Continuation.H"
#include "GlobalDefinitions.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
// A few declarations (see GlobalDefinitions.H)
// --> put them in a namespace
RCP<std::ostream> outFile;      // output file
ProfileType       profile;      // profile
std::stack<Timer> timerStack;   // timing stack


void runCoupledModel(RCP<Epetra_Comm> Comm);

// These are duplicated from test.C -> factorize and put them in a namespace
//------------------------------------------------------------------
RCP<std::ostream> outputFiles(RCP<Epetra_Comm> Comm);
RCP<Epetra_Comm>  initializeEnvironment(int argc, char **argv);
void              printProfile(ProfileType profile);

//------------------------------------------------------------------
int main(int argc, char **argv)
{
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);

	runCoupledModel(Comm);

	// print the profile
	if (Comm->MyPID() == 0)
		printProfile(profile);
	
    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	MPI_Finalize();	
}

//------------------------------------------------------------------
void runCoupledModel(RCP<Epetra_Comm> Comm)
{
	TIMER_START("Total time...");

	//------------------------------------------------------------------
	// Check if outFile is specified
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");	
	
	//------------------------------------------------------------------
	// Create parameter object for Ocean
	RCP<Teuchos::ParameterList> oceanParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());

	// Create parallelized Ocean object
	RCP<Ocean> ocean = Teuchos::rcp(new Ocean(Comm, oceanParams));

	//------------------------------------------------------------------
	// Create parameter object for Atmosphere
	RCP<Teuchos::ParameterList> atmosphereParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());

    // Create Atmosphere object
	// Let the ocean model dictate the horizontal resolution for the atmosphere
	std::shared_ptr<Atmosphere> atmos =
		std::make_shared<Atmosphere>(atmosphereParams);

	//------------------------------------------------------------------
    // Create parameter object for coupledmodel
	RCP<Teuchos::ParameterList> coupledmodelParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("coupledmodel_params.xml", coupledmodelParams.ptr());

	// Create CoupledModel
	std::shared_ptr<CoupledModel> coupledModel =
		std::make_shared<CoupledModel>(ocean, atmos, coupledmodelParams);

	//------------------------------------------------------------------
 	// Create parameter object for continuation
	RCP<Teuchos::ParameterList> continuationParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("continuation_params.xml", continuationParams.ptr());
	
	// Create Continuation
	Continuation<std::shared_ptr<CoupledModel>, RCP<Teuchos::ParameterList> >
		continuation(coupledModel, continuationParams);

	// Run continuation
	continuation.run();

	//------------------------------------------------------------------
	TIMER_STOP("Total time...");		
}

//------------------------------------------------------------------
// Auxiliary stuff
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
		if (map.first.compare(0,5,"_NOTIME_") != 0)
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
		if (map.first.compare(0,5,"_NOTIME_") == 0 )
		{
			counter++;
			std::stringstream s;
			s << " (" << counter << ")";
			LINE(s.str(), "", map.first.substr(5), ":", map.second[0], "",
				 map.second[1], "", map.second[2]);
		}
	
}
