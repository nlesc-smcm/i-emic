//=======================================================================
// Continuation test routine
//=======================================================================
#include <Epetra_config.h>
#include <iomanip>

#  ifdef HAVE_MPI
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif // HAVE_MPI


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "Vector.H"
#include "Atmosphere.H"
#include "OceanCont.H"
#include "ThetaStepper.H"
#include "Continuation.H"
#include "GlobalDefinitions.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
// A few declarations
RCP<std::ostream> outFile;               // output file
std::map<std::string, double> profile;   // profile map
RCP<std::ostream> outputFiles(RCP<Epetra_Comm> Comm);
RCP<Epetra_Comm>  initializeEnvironment(int argc, char **argv);
void printProfile(std::map<std::string, double> profile,
				  RCP<Epetra_Comm> Comm);

//------------------------------------------------------------------
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

	// Create Atmosphere model
	RCP<Atmosphere> atmos = rcp(new Atmosphere());
	atmos->test();
	atmos->computeJacobian();	
	getchar();

	// Create parameter object for continuation
	RCP<Teuchos::ParameterList> continuationParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("continuation_params.xml",
								continuationParams.ptr());
	// Create continuation
	Continuation<RCP<OceanCont>,
				 RCP<Vector>,
				 RCP<Teuchos::ParameterList> >
		continuation(ocean, continuationParams);

	continuation.Run();
	
	ocean->DumpState();

	printProfile(profile, Comm);
    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	Comm->Barrier();
	MPI_Finalize();
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
	if (Comm->MyPID() == 0)
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
void printProfile(std::map<std::string, double> profile,
				  RCP<Epetra_Comm> Comm)
{
	
	std::ostringstream profilefile;  // setting up a filename
	RCP<std::ostream> file;          // output file
	
	if (Comm->MyPID() == 0)
	{
		profilefile << "profile_"    << Comm->NumProc() <<  ".txt";
		INFO("profile for #procs = " << Comm->NumProc()
			 << " is written to "    << profilefile.str().c_str());
		file =
			Teuchos::rcp(new std::ofstream(profilefile.str().c_str()));	
	}
	else
	{
		file = Teuchos::rcp(new Teuchos::oblackholestream());
	}
	
	double sum = 0;
	(*file) << "=================================================================="
			<< std::endl;
	(*file) << "  Profile #CPU   = " << Comm->NumProc() << std::endl;
	for (std::map<std::string, double>::iterator it = profile.begin();
		 it != profile.end(); ++it)
	{
		sum += it->second;
		(*file) << " " << std::setw(35) << std::left << it->first
				<< std::setw(6)  << " -> "
				<< std::setw(12) << std::left << std::setprecision(8) << it->second
				<< std::endl;
	}
	(*file) << " " << std::setw(35) << std::left << " "  
			<< std::setw(6)  << " tot " 
			<< std::setw(12) << std::left << std::setprecision(8) << sum
			<< std::endl;
	(*file) << "=================================================================="
			<< std::endl;
}
