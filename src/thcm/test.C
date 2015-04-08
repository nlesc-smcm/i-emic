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
//#include "Rearranger.H"

#include "TRIOS_Domain.H"
#include "TRIOS_BlockPreconditioner.H"

using Teuchos::RCP;
using Teuchos::rcp;

typedef Epetra_Vector       VEC;
typedef Epetra_MultiVector  MVEC;
typedef Epetra_Operator     OPER;
typedef Epetra_CrsMatrix    CRSMAT;
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
