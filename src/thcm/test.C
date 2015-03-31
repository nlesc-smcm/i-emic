//==================================================================
//TODO: Factorize parameter setup and perhaps more
//==================================================================
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

#include <NOX.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_LinearSystem_Amesos.H>

#include <LOCA.H>
#include <LOCA_Epetra.H>

#include "THCM.H"
#include "OceanGrid.H"
#include "THCMdefs.H"

#include "Utils.H"

#include <sstream>
#include <fstream>
#include "test.H"

#define MAX(x,y) (x > y ? x : y)
#define MIN(x,y) (x < y ? x : y)

//----------------------------------------------------------------------
// Function declarations
//----------------------------------------------------------------------
extern "C"
{ 
	_SUBROUTINE_(write_data)(double*, int*, int*);  // file inout.f
} //extern

void displayStatus(double time, double dt, int niters);
void writeTHCMData(const Epetra_Vector &solution, Teuchos::RCP<Epetra_Comm> Comm);
Teuchos::RCP<std::ostream> outputFiles(Teuchos::RCP<Epetra_Comm> Comm);

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
	updateParametersFromXmlFile("test_params.xml", globalParamList.ptr());
	Teuchos::ParameterList &thcmList = globalParamList->sublist("THCM");
	thcmList.set("Parameter Name", "Time");
    //--------------------------------------------------------------------
	// Parameters -> top level:
	//  -creating RCP pointing to newly allocated ParameterList
	//  -obtaining list by dereferencing the RCP
	//--------------------------------------------------------------------
	Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
		Teuchos::rcp(new Teuchos::ParameterList);
	Teuchos::ParameterList &nlParams = *(nlParamsPtr.get());
	//--------------------------------------------------------------------
	// Set the nonlinear solver method
	//--------------------------------------------------------------------
	nlParams.set("Nonlinear Solver", "Line Search Based");
	nlParams.set("Convergence Tolerance", 1.0e-6);
	//--------------------------------------------------------------------
	// Set the printing parameters in the "Printing" sublist
	//--------------------------------------------------------------------
	Teuchos::ParameterList &printParams = nlParams.sublist("Printing");
	printParams.set("MyPID", Comm->MyPID());
	printParams.set("Output Precision", 3);
	printParams.set("Output Processor", 0);
	bool verbose = true;
	if (verbose)
		printParams.set("Output Information",
						NOX::Utils::OuterIteration +
						NOX::Utils::OuterIterationStatusTest +
						NOX::Utils::InnerIteration +
						//	NOX::Utils::LinearSolverDetails +
						NOX::Utils::Parameters +
						NOX::Utils::Details +
						NOX::Utils::Warning +
						NOX::Utils::Debug +
						NOX::Utils::TestDetails +
						NOX::Utils::Error);
	else
		printParams.set("Output Information",
						NOX::Utils::Warning +
						NOX::Utils::Error);
	
	//--------------------------------------------------------------------
	// Nonlinear solver parameters in the
	// "Line Search" and "Direction" sublists
	//--------------------------------------------------------------------
	Teuchos::ParameterList &searchParams = nlParams.sublist("Line Search");
	searchParams.set("Method", "Backtrack");
	searchParams.set("Max Iters", 10);
	Teuchos::ParameterList &backtrackParams = searchParams.sublist("Backtrack");
	backtrackParams.set("Default Step", 1.0);
	backtrackParams.set("Max Iters", 10);
	backtrackParams.set("Minimum Step", 1e-6);
	backtrackParams.set("Recovery Step", 1e-3);
	
	Teuchos::ParameterList &dirParams = nlParams.sublist("Direction");
	dirParams.set("Method", "Newton");
	Teuchos::ParameterList &newtonParams = dirParams.sublist("Newton");
	newtonParams.set("Forcing Term Method", "Type 2");
	newtonParams.set("Forcing Term Initial Tolerance", 1.0e-2);
	newtonParams.set("Forcing Term Maximum Tolerance", 1.0e-2);
	newtonParams.set("Forcing Term Minimum Tolerance", 1.0e-2);
	newtonParams.set("Rescue Bad Newton Solve", true);

    //-------------------------------------------------------------------
	// Sublist for the linear solver in the Newton method
	//-------------------------------------------------------------------

	Teuchos::ParameterList &lsParams = newtonParams.sublist("Linear Solver");
	lsParams.set("Aztec Solver", "GMRES");
	lsParams.set("Max Iterations", 1000);
	lsParams.set("Tolerance", 1e-3);
	lsParams.set("Output Frequency", 1);
	lsParams.set("Preconditioner", "None");
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
	soln->Scale(1.0e-10);
//	soln->PutScalar(0.001);
	//-------------------------------------------------------------------
	// Obtain Jacobian 
	//-------------------------------------------------------------------
	THCM::Instance().evaluate(*soln, Teuchos::null, true);
	Teuchos::RCP<Epetra_CrsMatrix> A = THCM::Instance().getJacobian();

    //-------------------------------------------------------------------
	// Setup timestep and time
	//-------------------------------------------------------------------

	double dt_min    = 1e-12;    // minimum timestep
	double dt_max    = 0.1;      // max timestep 
	double dt        = 1.0e-08;  // starting timestep 
	double dt_sc     = 1.10;     // scaling for the adaptive timestep
	double dt_sc2    = 1.01;     // scaling for a persistent increase
	double dt_sc3    = 2;        // scaling for a reset
	double t         = 0.0;      // start time
	double t_end     = 2*dt;     // end time
	
	// ------------------------------------------------------------------
	// Setup the solver interface
	// ------------------------------------------------------------------	
	Teuchos::RCP<SolverInterface> interface =
		Teuchos::rcp(new SolverInterface(*soln, dt));
	// ------------------------------------------------------------------	
	Teuchos::RCP<NOX::Epetra::Interface::Required> iReq        = interface;
	Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac        = interface;

	Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
		Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO
					 (printParams, lsParams, iReq, iJac, A, *soln));
	//-------------------------------------------------------------------
	// Use a NOX::Epetra::Vector for the initial guess in the
	// construction of a NOX::Epetra::Group
	//-------------------------------------------------------------------
	NOX::Epetra::Vector noxInitGuess(*soln, NOX::DeepCopy);
	std::cout << "UniqeGIDs: " << soln->Map().UniqueGIDs() << std::endl;
	Teuchos::RCP<NOX::Epetra::Group> grpPtr =
		Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq,
											noxInitGuess, linSys));	

	//-------------------------------------------------------------------
	// Setup some status tests and combine them into combo
	//-------------------------------------------------------------------
	int maxNewtonIters = 100;
	double newtonTolerance = 1.0e-5;
	Teuchos::RCP<NOX::StatusTest::NormF> testNormF =
	 	Teuchos::rcp(new NOX::StatusTest::NormF(newtonTolerance));
	Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters =
	 	Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));	
	Teuchos::RCP<NOX::StatusTest::Combo> combo =
	 	Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
	 											testNormF, testMaxIters));
    //-------------------------------------------------------------------
	// Create the solver
	//-------------------------------------------------------------------
	Teuchos::RCP<NOX::Solver::Generic> solver =
	 	NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
	//-------------------------------------------------------------------
	// Get the domain and create an OceanGrid to use some existing
	// infrastructure for calculating streamfunctions etc.
	//-------------------------------------------------------------------
	Teuchos::RCP<TRIOS::Domain> dom = THCM::Instance().GetDomain();
	Teuchos::RCP<OceanGrid> gridPtr =
		Teuchos::rcp(new OceanGrid(dom));
	// =================================================
	std::ostringstream filename;
	filename << "psidata" << t_end * 2 << ".txt";
	Teuchos::RCP<std::ostream> psidata =
		Teuchos::rcp(new std::ofstream(filename.str().c_str()));
	//-------------------------------------------------------------------
	// Do a few timesteps, solve the linear system in each step
	//-------------------------------------------------------------------
	NOX::StatusTest::StatusType status;

	double merPsiMax, barPsiMax;
	double uMax, vMax, wMax;
	int niters     = 0;
	int niters_min = 3;
	int niters_max = 6;
	int counter    = 0;
	bool pid0 = (Comm->MyPID() == 0);
	while (t < t_end)
	{
		//--------------------------------------------------------------
		// Solve
		//--------------------------------------------------------------
		
		const NOX::Epetra::Vector oldX = // get a copy of the old solution
			dynamic_cast<const NOX::Epetra::Vector &>(grpPtr->getX());
		
		std::cout << "====== Start solve" << std::endl;
		status = solver->solve();
		niters = solver->getNumIterations();

		//--------------------------------------------------------------
		// Calculate new timestep / reset if necessary
		// FACTORIZE THIS PART
		//--------------------------------------------------------------
		if (niters <= niters_min)
		{
			if (pid0)
				std::cout << "======= Increasing dt with a factor "
						  << dt_sc << std::endl;
			dt = MIN(dt * dt_sc, dt_max);
			interface->setTimeStep(dt);
		}
		else if (niters >= niters_max && niters < maxNewtonIters)
		{
			if (pid0)
				std::cout << "======= Decreasing dt with a factor "
						  << dt_sc << std::endl;
			dt = MAX(dt / dt_sc, dt_min);
			interface->setTimeStep(dt);
		}
		else if (niters >= maxNewtonIters)
		{
			dt = MAX(dt / dt_sc3, dt_min);
			if (pid0)
			{
				std::cout << "======= RESETTING " << std::endl;
				std::cout << " New dt: " << dt << std::endl;
			}
			interface->setTimeStep(dt);
			solver->reset(oldX); // Not updated solution oldX
			continue;
		}
		else
		{
			if (pid0)
				std::cout << "======= Increasing dt with a factor "
						  << dt_sc2 << std::endl;
			dt = MIN(dt * dt_sc2, dt_max);
			interface->setTimeStep(dt);
		}
		
		//--------------------------------------------------------------
		// Obtain solution
		//--------------------------------------------------------------
		const Epetra_Vector &newX = 
			(dynamic_cast<const NOX::Epetra::Vector &>(grpPtr->getX()))
			.getEpetraVector();
		//--------------------------------------------------------------
		// -Reset solver with new initial guess.
		// -Update the oldX in the solver interface.
		// -Update the time.
		//--------------------------------------------------------------
		solver->reset(grpPtr->getX());
		interface->setOldX(newX);
		gridPtr->ImportData(newX);
		merPsiMax = gridPtr->psimMax();
		barPsiMax = gridPtr->psibMax();
		uMax      = gridPtr->uMax();
		vMax      = gridPtr->vMax();
		wMax      = gridPtr->wMax();
		++counter;
		t += dt;
		
		if (pid0)
		{
			(*psidata) << counter << " " << niters << " "
					   << t << " " << merPsiMax << '\n';
			displayStatus(t, dt, niters);
			std::cout << "   max(Psi_m) = "
					  << merPsiMax << std::endl;
			std::cout << "   max(Psi_b) = "
					  << barPsiMax << std::endl;
			std::cout << "   max(u)     = "
					  << uMax << std::endl;
			std::cout << "   max(v)     = "
					  << vMax << std::endl;
			std::cout << "   max(w)     = "
					  << wMax << std::endl;
		}
					
		//--------------------------------------------------------------
		// Write the solution in traditional THCM format every 5 steps
		//--------------------------------------------------------------
		if (!(t < t_end) || ((counter % 5) == 0))
		{
			writeTHCMData(newX, Comm);
		}
	}
	THCM::Instance().writeParams();
    //------------------------------------------------------------------
	// Finalize MPI
	//------------------------------------------------------------------
	MPI_Finalize();
	return 0;
}

//***********************************************************************************
// Support functions
//***********************************************************************************

void displayStatus(double time, double dt, int niters)
{
	std::cout << "________________________________________________\n"
			  << "            DAYS  " << time/(1./(2.*365.))  << '\n'
			  << "            YEARS " << time*2               << '\n'
			  << "       Current dt " << dt                   << '\n'
			  << "     # Iterations " << niters               << '\n'
			  << "________________________________________________\n"
			  << std::endl;
}

//***********************************************************************************

void writeTHCMData(const Epetra_Vector &solution, Teuchos::RCP<Epetra_Comm> Comm)
{
	// This function will probably break (?) on a distributed memory system.
	// For now it is convenient.
	// Use some HYMLS functionality to gather the solution in the right way
	Teuchos::RCP<Epetra_MultiVector> fullSol =Utils::Gather(solution, 0);
	int filename = 3;
	int label    = 2;	
	int length   = fullSol->GlobalLength();
	double *solutionArray = new double[length]; 
	if (Comm->MyPID() == 0)
	{
		std::cout << "Writing to fort." << filename
				  << " at label " << label << "." << std::endl;
		(*fullSol)(0)->ExtractCopy(solutionArray); //  (*..)(0) ?? S>>??
		FNAME(write_data)(solutionArray, &filename, &label);
	}
	delete [] solutionArray;
}

//***********************************************************************************

Teuchos::RCP<std::ostream> outputFiles(Teuchos::RCP<Epetra_Comm> Comm)
{
	// Setup output files "fname_#.txt" for P==0 && P==1, other processes
	// will get a blackholestream.
	Teuchos::RCP<std::ostream> outFile;
	if (Comm->MyPID() < 2)
	{
		std::ostringstream infofile;  // setting up a filename
		infofile << "info_" << Comm->MyPID() << ".txt";
		std::cout << "info for P=" << Comm->MyPID() << " is written to "
				  << infofile.str().c_str() << std::endl;
		outFile = 
			Teuchos::rcp(new std::ofstream(infofile.str().c_str()));
	}
	else
		outFile = 
			Teuchos::rcp(new Teuchos::oblackholestream());
	return outFile;
}
