//=======================================================================
// Continuation test routines
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

#include <memory>
#include <vector>
#include <array>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_FancyOStream.hpp>

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
RCP<std::ostream> outFile;               // output file
ProfileType profile;                     // profile 

//------------------------------------------------------------------
void testVecWrap();
void testCoupling(int argc, char **argv);
void testOcean(int argc, char **argv);

//------------------------------------------------------------------
RCP<std::ostream> outputFiles(RCP<Epetra_Comm> Comm);
RCP<Epetra_Comm>  initializeEnvironment(int argc, char **argv);
void printProfile(ProfileType profile, RCP<Epetra_Comm> Comm);

//------------------------------------------------------------------
int main(int argc, char **argv)
{
//	testOcean(argc, argv);
	testCoupling(argc, argv);
}

//------------------------------------------------------------------
void testCoupling(int argc, char **argv)
{
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);

 	// Create parameter object for coupledmodel
	RCP<Teuchos::ParameterList> coupledmodelParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("coupledmodel_params.xml",
								coupledmodelParams.ptr());
	
	std::shared_ptr<CoupledModel> coupledModel =
		std::make_shared<CoupledModel>(Comm, coupledmodelParams);

 	// Create parameter object for continuation
	RCP<Teuchos::ParameterList> continuationParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("continuation_params.xml",
								continuationParams.ptr());
	
	// Create continuation
	Continuation<std::shared_ptr<CoupledModel>,
				 std::shared_ptr<SuperVector>,
				 RCP<Teuchos::ParameterList> >
		continuation(coupledModel, continuationParams);
	
	continuation.run();	
	
	printProfile(profile, Comm);
    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	Comm->Barrier();
	MPI_Finalize();	
}

void testOcean(int argc, char **argv)
{
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);

 	// Create ocean model Ocean
	//  (based on THCM):
	RCP<Ocean> ocean = rcp(new Ocean(Comm));

	// Create parameter object for continuation
	RCP<Teuchos::ParameterList> continuationParams =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("continuation_params.xml",
								continuationParams.ptr());
	// Create continuation
	Continuation<RCP<Ocean>,
				 RCP<SuperVector>,
				 RCP<Teuchos::ParameterList> >
		continuation(ocean, continuationParams);

	continuation.run();
	
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
	if (Comm->MyPID() < 1)
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
void printProfile(ProfileType profile, RCP<Epetra_Comm> Comm)
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
	(*file) << "==========================================================================="
			<< std::endl;
	(*file) << "  Profile #CPU   = " << Comm->NumProc() << std::endl;
	(*file) << " " << std::setw(35) << std::left << " "
			<< std::setw(6)  << " "
			<< std::setw(10) << std::left << "time"
			<< std::setw(2)  << " "
			<< std::setw(6)  << std::left << "calls"
			<< std::setw(2)  << " "
			<< std::setw(10) << std::left << "avg"		
			<< std::endl;

	for (ProfileType::iterator it = profile.begin();
		 it != profile.end(); ++it)
	{
		sum += it->second[0];
		(*file) << " " << std::setw(35) << std::left << it->first
				<< std::setw(6)  << " -> "
				<< std::setw(10) << std::left << std::setprecision(5) << it->second[0]
				<< std::setw(2)  << " "
				<< std::setw(6)  << std::left << std::setprecision(3) << it->second[1]
				<< std::setw(2)  << " "
				<< std::setw(10) << std::left << std::setprecision(5)
				<< it->second[0] / it->second[1]
				<< std::endl;
	}
	(*file) << " " << std::setw(35) << std::left << " "  
			<< std::setw(6)  << " tot " 
			<< std::setw(10) << std::left << std::setprecision(5) << sum
			<< std::endl;
	(*file) << "==========================================================================="
			<< std::endl;
}

//------------------------------------------------------------------
void testVecWrap()
{
	std::cout << "Testing the SuperVector Wrapper..." << std::endl;
	
	std::vector<double> vec1 = {1,2,3,4,5,6,7,8,9,10};
	std::vector<double> vec2 = {2,2,3,2,5,2,2,2,9,2};
	std::shared_ptr<std::vector<double> > spvec1 =
		std::make_shared<std::vector<double> >(vec1);
	std::shared_ptr<std::vector<double> > spvec2 =
		std::make_shared<std::vector<double> >(vec2);
	
 	SuperVector wrvec1(spvec1);
	SuperVector wrvec2(spvec2);

	std::cout << "v1 length: " <<  wrvec1.length() << std::endl;
	std::cout << "v2 length: " <<  wrvec2.length() << std::endl;

	std::cout << "v1 norm: " << wrvec1.norm() << std::endl;
	std::cout << "v2 norm: " << wrvec2.norm() << std::endl;
		
	std::cout << "(v1,v2): " <<	wrvec1.dot(wrvec2) << std::endl;
	
	wrvec1.update(4, wrvec2, 3);
	wrvec1.print();

	wrvec2.scale(3);
	wrvec2.print();
	
	std::cout << "Testing the SuperVector Wrapper...done" << std::endl;
}
