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

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>

#include "THCM.H"
#include "THCMdefs.H"

#include <sstream>
#include "GlobalDefinitions.H"

// for gethostname in pardebug
#include <unistd.h>

#include "TRIOS_Domain.H"
#include "TRIOS_BlockPreconditioner.H"

using Teuchos::RCP;
using Teuchos::rcp;

typedef Epetra_Vector       VEC;
typedef Epetra_MultiVector  MVEC;
typedef Epetra_Operator     OPER;
typedef Epetra_CrsMatrix    CRSMAT;

Teuchos::RCP<std::ostream> outputFiles(Teuchos::RCP<Epetra_Comm> Comm);
void parDebug(Teuchos::RCP<Epetra_Comm> Comm);
RCP<std::ostream> outFile;
//----------------------------------------------------------------------
// THCM is a singleton, there can be only one instance at a time. 
// As base class Singleton is templated we must instantiate it    
// using the macro defined in Singleton.H:
//----------------------------------------------------------------------
_INSTANTIATE_(THCM);

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
	// -------------------------------------------------------------------
	// Setup THCM parameters:
	// -------------------------------------------------------------------
	RCP<Teuchos::ParameterList> globalParamList =
		rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("thcm_params.xml", globalParamList.ptr());
	Teuchos::ParameterList &thcmList = globalParamList->sublist("THCM");
	thcmList.set("Parameter Name", "Time");
    //-------------------------------------------------------------------
	// Create THCM object
	//-------------------------------------------------------------------	
	RCP<THCM> ocean = rcp(new THCM(thcmList, Comm));
	//-------------------------------------------------------------------
	// Obtain solution vector from THCM
	//-------------------------------------------------------------------	
	RCP<VEC> soln = THCM::Instance().getSolution();
	//-------------------------------------------------------------------
	// Initialize solution vector
	//-------------------------------------------------------------------
	soln->Random();
	soln->Scale(1.0e-3);
	INFO("Initialized solution vector");
	//-------------------------------------------------------------------
	// Obtain Jacobian
	//-------------------------------------------------------------------
	THCM::Instance().evaluate(*soln, Teuchos::null, true);
	RCP<CRSMAT> A = THCM::Instance().getJacobian();

	//-----------------------------------------------------------------------------
	// Block preconditioner 
	//-----------------------------------------------------------------------------
	
	Teuchos::RCP<Teuchos::ParameterList> solverParams =
		Teuchos::rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams.ptr());	
	RCP<TRIOS::Domain> domain = THCM::Instance().GetDomain();
	INFO("====== CONSTRUCT PRECONDITIONER ===========");
	RCP<TRIOS::BlockPreconditioner> blockPrec =
		Teuchos::rcp(new TRIOS::BlockPreconditioner(A, domain, *solverParams));
	INFO("=========== TEST PRECONDITIONER ===========");
	blockPrec->Test();
	
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
