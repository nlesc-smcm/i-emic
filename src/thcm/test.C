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

#include "THCM.H"
#include "THCMdefs.H"
#include <sstream>
#include <fstream>

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
	Teuchos::RCP<Epetra_MpiComm> Comm =
		Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	Teuchos::RCP<Epetra_SerialComm> Comm =
		Teuchos::rcp(new Epetra_SerialComm());
#endif
	// -------------------------------------------------------------------
	// Setup THCM parameters:
	// -------------------------------------------------------------------
	Teuchos::RCP<Teuchos::ParameterList> globalParamList =
		Teuchos::rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("thcm_params.xml", globalParamList.ptr());
	Teuchos::ParameterList &thcmList = globalParamList->sublist("THCM");
	thcmList.set("Parameter Name", "Time");
    //-------------------------------------------------------------------
	// Create THCM object
	//-------------------------------------------------------------------	
	Teuchos::RCP<THCM> ocean = Teuchos::rcp(new THCM(thcmList, Comm));
	//-------------------------------------------------------------------
	// Obtain solution vector from THCM
	//-------------------------------------------------------------------	
	Teuchos::RCP<Epetra_Vector> soln = THCM::Instance().getSolution();
	//-------------------------------------------------------------------
	// Initialize solution vector
	//-------------------------------------------------------------------
	soln->Random();
	soln->Scale(1.0e-8);
	INFO("Initialized solution vector");
	//-------------------------------------------------------------------
	Teuchos::RCP<Epetra_Vector> RHS = Teuchos::rcp(new Epetra_Vector(soln->Map()));
	// Calculate rhs and store it in RHS
	THCM::Instance().evaluate(*soln, RHS, false);
	// Obtain Jacobian 
	//-------------------------------------------------------------------
	THCM::Instance().evaluate(*soln, Teuchos::null, true);
	Teuchos::RCP<Epetra_CrsMatrix> A = THCM::Instance().getJacobian();
	double nrm[1];
	soln->Norm2(nrm);
	double nrm2[1];
	RHS->Norm2(nrm2);
	std::cout << "||x||_2    = " << nrm[0] << std::endl;
	std::cout << "||RHS||_2  = " << nrm2[0] << std::endl; 
	std::cout << "||A||_inf  = " << A->NormInf() << std::endl;
	std::cout << "EXITING" << std::endl;
    //------------------------------------------------------------------
	// Finalize MPI
	//------------------------------------------------------------------
	MPI_Finalize();
	return 0;
}
