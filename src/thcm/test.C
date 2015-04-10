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
	// Uncomment this when parallel debugging with gdb
	//  - It enters an infinite while loop so you can couple a task to
	//    gdb with the -p option.
	//  - Then with gdb step until you get in the condition of the while loop
	//    and switch the integer to a value that negates the condition (set var i = 7)
	//  - On Ubuntu you need root privileges for this to work.
	// -------------------------------------------------------------------
	// parDebug(Comm);
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
	// Initialize X, MultiRHS and RHS
	// ------------------------------------------------------------------	
	RCP<MVEC> X    = rcp(new MVEC(soln->Map(), 1));
 	RCP<MVEC> MRHS = rcp(new MVEC(soln->Map(), 1));
	// Create non-owning rcp from Epetra_MultiVector MRHS
	// pointing to first Epetra_Vector in MRHS
	RCP<VEC>  RHS = rcp((MRHS->operator()(0)), false);
	double nrm1[1];
	double nrm2[1];
	MRHS->Norm2(nrm1);
	INFO("Norm MRHS before calculation " << nrm1[0]);
	// Calculate rhs and store it in RHS
	THCM::Instance().evaluate(*soln, RHS, false);
	MRHS->Norm2(nrm2);
	INFO("Norm MRHS after calculation " << nrm2[0]);
	//-------------------------------------------------------------------
	// Obtain Jacobian
	//-------------------------------------------------------------------
	THCM::Instance().evaluate(*soln, Teuchos::null, true);
	RCP<CRSMAT> A = THCM::Instance().getJacobian();
	//-------------------------------------------------------------------
	// Rearrange the Jacobian so we can solve it with less iterations
	//-------------------------------------------------------------------
	// Rearranger rearr;
	// rearr.setMatrix(A);
	// rearr.buildOrdering();
	// rearr.setBlockOperator();
	// rearr.fillBlocks();
	// rearr.test();

    //-------------------------------------------------------------------
	// Belos::LinearProblem setup
	//-------------------------------------------------------------------
	RCP<Belos::LinearProblem<double, MVEC, OPER>> problem =
		rcp(new Belos::LinearProblem<double, MVEC, OPER>(A, X, MRHS));
		
	//-----------------------------------------------------------------------------
	// Block preconditioner 
	//-----------------------------------------------------------------------------
	
	Teuchos::RCP<Teuchos::ParameterList> solverParams =
		Teuchos::rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams.ptr());	
	RCP<TRIOS::Domain> domain = THCM::Instance().GetDomain();
	DEBUG(*solverParams);
	RCP<Ifpack_Preconditioner> blockPrec =
		Teuchos::rcp(new TRIOS::BlockPreconditioner(A, domain, *solverParams));
	blockPrec->Compute();
	RCP<Belos::EpetraPrecOp> belosPrec = rcp(new Belos::EpetraPrecOp(blockPrec));
	problem->setRightPrec(belosPrec);
	
	//-------------------------------------------------------------------
	// Belos parameter setup
	//-------------------------------------------------------------------
	Teuchos::RCP<Teuchos::ParameterList> belosList =
		Teuchos::rcp(new Teuchos::ParameterList());
	belosList->set("Block Size", 100);
	belosList->set("Maximum Iterations", 10000);
	belosList->set("Convergence Tolerance", 1e-4);
	belosList->set("Verbosity", Belos::FinalSummary);
	//-------------------------------------------------------------------
	// Belos block GMRES setup
	//-------------------------------------------------------------------
	Belos::BlockGmresSolMgr<double, MVEC, OPER> belosSolver(problem, belosList);
	//--------------------------------------------------------------------
	INFO("Before setProblem()")
		bool set = problem->setProblem();	
	INFO("After setProblem()")
		belosSolver.solve();
	
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

void parDebug(Teuchos::RCP<Epetra_Comm> Comm)
{
	int i = 0;
	int myPID = Comm->MyPID();
	char hostname[256];
	gethostname(hostname, sizeof(hostname));
	if (Comm->MyPID() == 0)
	{
		printf("PID %d on %s, CPU%d is ready for attach\n", getpid(), hostname, myPID);
		fflush(stdout);
		getchar();
	}
	Comm->Barrier();
}
