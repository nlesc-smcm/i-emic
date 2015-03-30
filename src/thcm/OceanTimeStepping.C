/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

/************************************************************************
//OCEAN MODEL 
************************************************************************/

#include "Epetra_Time.h"
#include "Epetra_Vector.h"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX.H"
#include "LOCA_Parameter_Vector.H"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <sstream>

#include "HYMLS_MatrixUtils.H"

// global THCM definitions
#include "globdefs.H"

//User's application specific files 
#include "THCM.H"
#include "DefaultParams.H"
#include "TRIOS_DefaultParams.H"
#include "OceanModel.H"

#include "Teuchos_StandardCatchMacros.hpp"
#include "ThetaStepperEvaluator.H"
#include "ImplicitTimeStepper.H"
#include "OceanGrid.H"
#include "Teuchos_Array.hpp"

#include "HYMLS_Tools.H"

extern "C"
{ 
	_SUBROUTINE_(write_data)(double*, int*, int*);  // file inout.f
	
} //extern

// some constants (these should be consistent with
// those in usr.F90, we're to lazy to get them from
// there ;)
const double r0dim = 6370.0e3; //=== Radius of ...
const double udim  = 0.1;      //=== Some velocity scale ...

// a function for nice output of time data
string time_out(double t);

// THCM is a singleton, there can be only one instance at a time. 
// As base class Singleton is templated we must instantiate it    
// using the macro defined in Singleton.H

_INSTANTIATE_(THCM)    // ====== See Singleton.H

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
#endif

	int ierr = 0;   
  
////////////////////////////////////////
// Create the Parallel  Communicator  //
////////////////////////////////////////

#ifdef HAVE_MPI
	Teuchos::RCP<Epetra_MpiComm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	Teuchos::RCP<Epetra_SerialComm> Comm= Teuchos::rcp(new Epetra_SerialComm());
#endif

	//Get process ID and total number of processes
	int MyPID   = Comm->MyPID();
	int NumProc = Comm->NumProc();
	////////////////////////////////////////////////////////////////////////////
	// open static file streams (seperate for each process)                   //
	// Note: when using many processors it is convenient                      //
	// to make most of them dump their output to Teuchos::oblackholestreams   //
	// instead.                                                               //
	////////////////////////////////////////////////////////////////////////////

	if (MyPID==0)
	{
		std::ostringstream infofile;  //=== setting up infofile
		infofile  << "info_"<< MyPID << ".txt";
		std::cout << "info is written to " << infofile.str().c_str() << std::endl;
//		info = Teuchos::rcp(new std::ofstream(infofile.str().c_str()) );
		info = Teuchos::rcp(&std::cout, false); //wrapping std::cout in a non-owning rcp
	}
	else
	{
		info = Teuchos::rcp(new Teuchos::oblackholestream());
	}
  
	(*info) << std::setw(15) << std::setprecision(15);

	Teuchos::RCP<std::ostream> outstream;

	if (MyPID==0)
	{
		outstream = Teuchos::rcp(new std::ofstream("timestepping.out"));
	}
	else
    {
		outstream = Teuchos::rcp(new Teuchos::oblackholestream());
    }
	//=== Set format 
	(*outstream) << std::setw(15) << std::setprecision(15);
	
	//==================== DEBUGGING STUFF ===========================
#ifdef DEBUGGING
	
#if 1 // seperate files, set this manually here
	std::ostringstream debfile;
	debfile << "debug_" << MyPID << ".txt";
	std::cout << "debug output is written to " << debfile.str().c_str() << std::endl;
	Teuchos::RCP<std::ofstream> debug = 
		Teuchos::rcp(new std::ofstream(debfile.str().c_str())); // modified _erik
	(*debug) << std::setw(15) << std::setprecision(15);
#else // dump everything into one file
	debug = info;
	outstream = info;
	std::cout << "debug output is written to " << infofile.str().c_str() << std::endl;
#endif
    //== Note that HYMLS::Tools::InitializeIO_std() redirects std::cout
	HYMLS::Tools::InitializeIO_std(Comm, info, debug);	
#else	
//	HYMLS::Tools::InitializeIO_std(Comm, info); 
#endif 
	DEBUG("*********************************************");
	DEBUG("* Debugging output for process " << MyPID);
	DEBUG("* To prevent this file from being written,  ");
	DEBUG("* omit the -DDEBUGGING flag when compiling. ");
	DEBUG("*********************************************\n\n");
	INFO( *Comm);
	//================================================================

	// the try block starts here to make sure the streams and comms are
	// deleted last
	
	try {	
		////////////////////////////////////////////////////////
		// Setup Parameter Lists                              //
		////////////////////////////////////////////////////////

		INFO("Construct Parameter List...");
  
		// read parameters from files
		// we create one big parameter list and further down we will set some things
		// that are not straight-forward in XML (verbosity etc.)
		
		Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new Teuchos::ParameterList);
		DefaultParams::THCM(*paramList);

		// ====================================================================
		// Note that we now use .ptr() to get a safe wrapper around the
		// raw C++ pointer to the underlying object, instead of using .get().
		// ====================================================================
		
		Teuchos::updateParametersFromXmlFile("thcm_params.xml",paramList.ptr()); 
		DefaultParams::TimeStepping(*paramList);
		Teuchos::updateParametersFromXmlFile("transient_params.xml",paramList.ptr());
		DefaultParams::NOX(*paramList);
				
		TRIOS::DefaultParams::LinearSolver(paramList->sublist("NOX")
										   .sublist("Direction")
										   .sublist("Newton"));
		TRIOS::DefaultParams::BlockPreconditioner(paramList->sublist("NOX")
												  .sublist("Direction")
												  .sublist("Newton")
												  .sublist("Linear Solver"));
		TRIOS::DefaultParams::HYMLS(paramList->sublist("NOX")
									.sublist("Direction")
									.sublist("Newton")
									.sublist("Linear Solver"));
		Teuchos::updateParametersFromXmlFile("solver_params.xml",paramList.ptr());
		
		// Get the THCM sublist
		Teuchos::ParameterList& thcmList = paramList->sublist("THCM");

		// Get the "Solver" parameters sublist to be used with NOX Solvers
		Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

		// list to control output of NOX
		Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");

		//Create the "Direction" sublist for the "Line Search Based" solver
		Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");

		//Create the "Line Search" sublist for the "Line Search Based" solver
		Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");

		//Create the "Direction" sublist for the "Line Search Based" solver
		Teuchos::ParameterList& newtParams = dirParams.sublist("Newton");

		//Create the "Linear Solver" sublist for the "Direction' sublist
		Teuchos::ParameterList& lsParams = newtParams.sublist("Linear Solver");
		
		// Get the "Time Stepping" sublist
		Teuchos::ParameterList& transParams = paramList->sublist("Time Stepping");


		
		
		// ======================================================================
		// .get() returns either the value of "name" or the
		// "default" value in the second argument. In the latter case
		// "name" is set with the default value.
		// ======================================================================
		
		double hdim = thcmList.get("Depth hdim", 5000.0);
		std::string stepperScheme = transParams.get("Scheme","Theta");
		bool gradual_startup = transParams.get("Gradual Start-Up",false);

		std::string startup_param = "None";
		double startup_rate  = 0.0;
		double startup_value, startup_max;
    
		if (gradual_startup)
		{
			startup_param = transParams.get("Start-Up Parameter","Combined Forcing");
			startup_rate  = transParams.get("Start-Up Rate",0.1);
			startup_max   = thcmList.sublist("Starting Parameters").get(startup_param,0.0);
			// =========================================================================
			// So there are two startup parameters, one is called Start-Up Parameter
			// and lives in transient_params.xml, the other one is called startup_param
			// and is added to thcmList for some reason...
			// =========================================================================
		}
		
		///////////////////////////////////////////////////////////
		// Setup the Problem                                     //
		///////////////////////////////////////////////////////////

		transParams.set("Output Stream", outstream);
    
		double t_start    = transParams.get("Start Time", 0.0);
		double t_end      = transParams.get("End Time", 1.0);
		double lteTol     = transParams.get("Error Tolerance", 1.0e-3);
		std::string steadyMon  = transParams.get("Steady-State Monitor", "Norm of F");
		double steadyTol  = transParams.get("Steady-State Tolerance", 0.0);
		std::string stepChoice = transParams.get("Step Size Control", "Constant");

		// step size change factors for "Constant" mode
		double red_fac = transParams.get("Failed Step Reduction Factor", 0.5);
		double inc_fac = transParams.get("Successful Step Increase Factor", 1.5);
    
		// courant-number for CFL check
		double CN = transParams.get("Courant Number", 0.9);
    
		// if the flow stays under cflmin for maxslow steps we increase the time-step
		int    nslow   = 0;
		int    maxslow = transParams.get("Slow Flow Delay", 10);
		double slw     = transParams.get("Slow Flow Condition", 0.4);
		// PQ: What's the condition?
		double cflmin  = slw*std::abs(CN);

		//================= THCM SETUP ===============================
    	// Set up the THCM interface to allow calls to
		// residual (RHS) and Jacobian evaluation routines.
		// THCM is implemented as a Singleton, that means there
		// is only a single instance which should be referenced 
		// using THCM::Instance()

		INFO("Initialize THCM...");

		//== Ocean is a THCM object, providing an interface with
		//== the original FORTRAN code.
		Teuchos::RCP<THCM> ocean = Teuchos::rcp(new THCM(thcmList, Comm));
		
		//== EpetraFactory:
		Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
			Teuchos::rcp(new LOCA::Epetra::Factory);
		
		//== GlobalData:
		Teuchos::RCP<LOCA::GlobalData> globalData =
			LOCA::createGlobalData(paramList, epetraFactory);
		
		thcmList.set("Parameter Name","Time");

		//== Obtain preconditioning type, standard is Ifpack, in
		//== solver_params.xml it is probably "User Defined".
		std::string PrecType = lsParams.get("Preconditioner", "Ifpack");
    
		//== Obtain the Linear Solver preconditioning ParameterList.
		Teuchos::RCP<Teuchos::ParameterList> prec_params = Teuchos::null;
		if (PrecType == "User Defined") 
		{
			prec_params = Teuchos::rcp(&lsParams,false);
		}

		//Create the model interface
		Teuchos::RCP<OceanModel> model = 
			rcp(new OceanModel(thcmList,globalData,prec_params));


		// this vector defines what parameters the problem depends on
		// (or at least which parameters may be varied)
		// and gives initial values for them which overwrite
		// the settings made in usrc.F90::stpnt(). 
		Teuchos::RCP<LOCA::ParameterVector> pVector
			= model->getParameterVector();

		pVector->setValue("Time", t_start);

		//Get the vector from the problem
		Teuchos::RCP<Epetra_Vector> soln = model->getSolution();

		//Initialize solution
		soln->PutScalar(0.0);

		// check for starting solution
		//== In thcm_params.xml it is probably "None"
		std::string StartConfigFile = thcmList.get("Starting Solution File","None");
            
		if (StartConfigFile != "None")
		{
			INFO("Read Start File...");   
			soln = model->ReadConfiguration(StartConfigFile,*pVector);
			try 
			{
				t_start = pVector->getValue("Time");
			} catch (...) {Error("no time info found in start file",__FILE__,__LINE__);}
			transParams.set("Start Time",t_start);
			model->setParameters(*pVector);
		}

//construct Jacobian A and W = A + 1/dt B to have a closer look at them
#if 1
/*
  double dt_ = transParams.get("Step Size",0.1);
  RCP<Epetra_CrsMatrix> A = model->getJacobian();
  model->computeJacobian(*soln,*A);
  MatrixUtils::DumpMatrix(*A, "A.txt");
  model->computeShiftedMatrix(0.0,1.0,*soln,*A);
  MatrixUtils::DumpMatrix(*A, "B.txt");
  model->computeShiftedMatrix(1.0,1.0/dt_,*soln,*A);
  MatrixUtils::DumpMatrix(*A, "W.txt");
  std::ofstream ofs("F.txt");
  ofs << F ;
  ofs.close();
*/
		Epetra_Vector F(*soln);
		model->computeF(*soln, F, OceanModel::Residual);
		double nrm;
		soln->Norm2(&nrm);
		std::cout << "||u|| of starting solution: " << nrm << std::endl;
		F.Norm2(&nrm);
		std::cout << "||f|| of starting solution: " << nrm << std::endl;
#endif		

		INFO("Construct the Implicit Time Stepper...");

		//== Create an ImplicitTimeStepper object using
		//== OceanModel model and the "Time Stepping" and "NOX"
		//== parameters.

		Teuchos::RCP<AbstractTimeStepper> stepper =
			Teuchos::rcp(new ImplicitTimeStepper(model, transParams, nlParams));

		// get a grid-representation of the solution.   
		// For the implicit solver this is used only for
		// post-processing, whereas the explicit schemes
		// do the actual computations on the grid

		Teuchos::RCP<OceanGrid> grid = model->getGrid();

		double t           = t_start;
		bool   increase_dt = false;
		double dt          = transParams.get("Step Size",0.1);

		//== The overloaded constructors for Epetra_Vector seem to have changed :s.
		//== Omitting the Copy bit and just using the plain Epetra_Vector with copy
		//== construction.
		Epetra_Vector soln_old(soln->Map(), Copy);

		double dt_min = transParams.get("Minimum Step Size",(double)(dt/10));
		double dt_max = transParams.get("Maximum Step Size",dt);
		transParams.set("Max Num Steps",(int)((t_end-t_start)/dt));
		int maxsteps  = transParams.get("Max Num Steps",(int)((t_end-t_start)/dt));

		int step = 1;
		double maxEffort = 0.9; // prevent increasing step size if more than
		// 90% of the maximum Newton steps were needed.
		double PsiMaxOld = -1.0;

		//== in timestepping.out gradual_startup is probably 0
		if (gradual_startup)
		{
			pVector->setValue(startup_param,0.0);
			model->setParameters(*pVector);
			startup_value = 0.0;
		}

		INFO("Start Time integration");
		Epetra_Time timer(*Comm);

		model->printSolution(*soln,t);

#if 0
		{
			soln_old = *soln;
			INFO("TRIAL STEP...");
			INFO("dt = " << dt);
			double nrm;
			soln_old.Norm2(&nrm);
			INFO("||x|| of input vector: " << nrm);
			Epetra_Vector F(*soln);
			//== Calling inherited Loca::Epetra::ModelEvaluatorInterface::computeF()
			//== to compute residual
			model->computeF(soln_old, F, OceanModel::Residual);
			F.Norm2(&nrm);
			INFO("||f|| of input vector: " << nrm);
			INFO("do step");
			//== Calling ImplicitTimeStepper::Step() 
			stepper->Step(soln_old, 0, *soln, dt);
			Error("TRIAL STEP ONLY",__FILE__,__LINE__);
		}
#endif

		while (t < t_end)
		{
			if (step > maxsteps)
			{
				(*outstream) << "Maximum number of steps exceeded, stopping...\n";
				break;
			}
  
			(*outstream) << std::endl;
			(*outstream) << "##############################################\n";
			(*outstream) << "start time step " << step << ": t=" << time_out(t) << std::endl;
 
			// make sure that we hit t_end exactly  
			// If the time-step would become too    
			// small, we allow the stepper to       
			// go a little too far (may want to fix 
			// this)
			if (t + dt > t_end)
			{
				dt = std::max(dt_min, t_end-t);
			}

			(*outstream) << "step size dt = " << time_out(dt) << std::endl;
			(*outstream) << "##############################################\n";
			(*outstream) << std::endl;

			soln_old = *soln;
			bool success = stepper->Step(soln_old, t, *soln, dt);
						
			(*outstream) << "##############################################\n";
			(*outstream) << "returned from Step, success = " << success << std::endl;
			(*outstream) << "##############################################\n";
			double nrmu;
			soln->Norm2(&nrmu);
			(*outstream) << "    NORM  ||u|| = " << nrmu << std::endl;
			(*outstream) << std::endl;
			
			increase_dt = success;

			if (!success)
			{
				if (dt == dt_min)
				{
					(*outstream) << "Time stepping failed at minimum step-size, stopping run...\n";
					break;
				}
				else
				{
					(*outstream) << "Step failed!\n";
				}
			}

			(*outstream) << std::endl;
			(*outstream) << "##############################################\n";

			double lte = -1.0;  //== local truncation error
  
			double effortEst = stepper->getEffortEstimate();
			(*outstream) << "Step effort: " << effortEst * 100 << "%\n";
			
			// success only indicates wether the nonlinear solver converged,
			// if we have more criteria we must check them now:
			if (success && (stepper->hasErrorEstimate())) 
			{
				lte = stepper->getErrorEstimate();
				(*outstream) << "Local truncation error tau: " << lte << std::endl;

				if (success && (stepChoice == "Adaptive"))
				{
					if (lte > lteTol)
					{
						if (dt > dt_min)
						{
							success = false;
							(*outstream) << "Step rejected: estimate of LTE (" << lte << ")\n";
							(*outstream) << "               does not satisfy tolerance (" << lteTol << ")\n";
						}
						else
						{
							(*outstream) << "WARNING: step failed to achieve desired accuracy, \n";
							(*outstream) << "         but as the step-size is already at its   \n";
							(*outstream) << "         minimum, we accept it anyway.            \n";
						}
					}
				}
			}    

			// advance one step:
			if (success)
			{
				t += dt;
				pVector->setValue("Time", t);
				model->setParameters(*pVector);
				model->printSolution(*soln,t);
				step++;
				(*outstream) << "Step was successful!\n";
				if (gradual_startup)
				{
					startup_value = std::min(startup_max, startup_value + dt * startup_rate);
					if (startup_value == startup_max) gradual_startup = false;
					pVector->setValue(startup_param, startup_value);
					model->setParameters(*pVector);
				}
			}

			// check stability conditions
    
			// the stepper may have some restraint, for instance
			// explicit steppers require a parabolic stability condition
			// to be satisfied:
			double stabCond1 = stepper->getStabCond();
    
			// we always obey the CFL condition, you can use the Courant number
			// to control/disable this
			double stabCond2 = grid->cflCond();
 
			// adjust time step-size:
			if (dt < cflmin * stabCond2 && effortEst < maxEffort)
			{
				nslow++;
			}
			else
			{
				nslow=0;
			}
    
			if (nslow<maxslow)
			{
				increase_dt = false;
			}
			else
			{
				if (effortEst<maxEffort)
				{
					increase_dt = success;
					nslow = 0;
				}
				else
				{
					(*outstream) << "step size not increased because system to hard to solve\n";
					increase_dt = false;
				}
			}
    
			if (stepChoice=="Constant") // step-size control only based on nonlinear solver convergence
			{
				if (increase_dt)
				{
					if (dt<dt_max)
					{
						(*outstream) << "increasing step size for next step\n";
					}
					dt = std::min(dt_max,inc_fac*dt);
				}
				else if (!success) // reduce step size and try again
				{
					(*outstream) << "reducing step size\n";
					dt = std::max(dt_min, red_fac*dt);
					nslow=0;
				}
			}
			else if (stepChoice=="Adaptive") // adaptivity based on error estimate,
			{
				if (lte!=-1.0)
				{
					if ((lte>lteTol)||(lte<lteTol/4))
					{
						// compute new step size
						double newdt = sqrt(lteTol/(2*lte))*dt;
						newdt = std::max(std::min(newdt,2*dt),dt/2);
						newdt = std::min(newdt,dt_max);
						newdt = std::max(newdt,dt_min);
						DEBVAR(newdt);
						dt=newdt;
					}
				}
				else // same as "Constant". hasErrorEstimate==false may 
				{  // simply mean that this was the first step or a restart
					if (!success) 
					{
						(*outstream) << "reducing step size\n";
						dt = std::max(dt_min,dt/2);
					}
				}
			}
    
			if (stabCond1>0) // otherwise not required/implemented/available
			{
				(*outstream) << "stability condition: dt<"<<stabCond1;
				(*outstream) << " ("<<time_out(stabCond1)<<")\n";
				if (dt>stabCond1)
				{
					(*outstream) << "reducing time step because of stability condition\n";
					//TODO: the 0.9 shouldn't be hard-coded
					dt = 0.9*stabCond1;
				}
				DEBVAR(dt)
					}
			(*outstream) << "CFL condition: dt<"<<stabCond2<<" ("<<time_out(stabCond2)<<")\n";
			if (CN>0)
			{
				DEBVAR(dt);
				DEBVAR(stabCond2);
				if (dt>CN*stabCond2)
				{
					(*outstream) << "restraining time step because of CFL condition\n";
					(*outstream) << "(CFL condition: "<<stabCond2<<", Courant-number: "<<CN<<")\n";
					dt = CN*stabCond2;
				}
				DEBVAR(dt)
					}
			(*outstream) << "##############################################\n";
			if (!success)
			{
				*soln = soln_old;
				stepper->reset();
			}

			// check for convergence to a steady state (when dU/dt->0)
			if (success&&steadyMon=="Norm of F")
			{
				const Epetra_Vector& fn = stepper->getF();
				double nrmF;
				CHECK_ZERO(fn.Norm2(&nrmF));
				(*outstream) << "Norm of F(u): "<<nrmF<<std::endl;
				if (nrmF<steadyTol)
				{
					(*outstream) << "Convergence to a steady-state has been detected, stopping here!"<<std::endl;
					break;
				}
			}

			// check for convergence to a steady state (when dU/dt->0)
			if (success&&steadyMon=="Overturning")
			{
				// compute meridional streamfunction on the 3D 
				// representation of the solution:
				double PsiMax = (grid->psimMax() - grid->psimMin());
				// compute PsiMax in physical units (Sv):
				// we use some constants from thcm,
				// which we defined globally above
				const double transc = r0dim*hdim*udim;

				PsiMax = PsiMax*transc*1e-6;
				(*outstream) << "\nMaximum meridional overturning streamfunction: "<<PsiMax<<" Sv"<<std::endl;
				double eps = std::abs((PsiMax-PsiMaxOld)/(PsiMaxOld*dt));
				if (PsiMaxOld>0) // in first step don't check
				{
					(*outstream) << "convergence monitor: dPsi/dt="<<eps<<std::endl;
					DEBVAR(PsiMax);
					DEBVAR(*grid);
					if (eps<steadyTol)
					{
						(*outstream) << "Convergence to a steady-state has been detected, stopping here!"<<std::endl;
						break;
					}
				}
				PsiMaxOld=PsiMax;      
			}
    
		}//while loop

		double time = timer.ElapsedTime();
		(*info) << "\n*******************************************\n";
		(*info) << "Elapsed time during entire run: "<<(int)(time/60)<<" min\n";
		(*info) << "*******************************************\n\n";

		(*outstream) << std::endl;
		(*outstream) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		(*outstream) << "Time-stepping run finished, store final solution and finish...\n";
		(*outstream) << "Elapsed time during entire run: "<<(int)(time/60)<<" min\n";
		(*outstream) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
		(*outstream) << std::endl;

// TODO: the openDX file is finished automatically by the OceanModel destructor,
//       but the final state is generally not stored!

///////////////////////////////////////////////////////////////

		// write solution in native THCM format
		int filename = 3;
		int label = 2;


		// gather the final solution on the root process (pid 0):
		Teuchos::RCP<Epetra_MultiVector> fullSol =HYMLS::MatrixUtils::Gather(*soln,0);

		model->WriteConfiguration("FinalConfig.txt",*pVector,*soln);

		if (MyPID == 0 )
		{
			int ndim = fullSol->GlobalLength();
			double *solutionbig = new double[ndim];

			(*fullSol)(0)->ExtractCopy(solutionbig);
			INFO(" Time-stepping run finished, store solution...");
			FNAME(write_data)(solutionbig, &filename, &label);
			INFO("done!");
		}

		(*outstream) << "Final Parameters: "<<std::endl;    
		(*outstream) << *paramList << std::endl;


// The HDF5 file is closed when OceanModel is deleted => OceanOutput is 
// deleted, but this happens only at the end of main because of the stepper hook.
// we need to make sure everyone closes the file before calling MPI_Finalize()
// TODO: This is a bug, we probably violated one or more of the 10 Teuchos::RCP-commandments
		model->finishOutput();

	} TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,ierr)
//}catch(...){std::cout << "caught something!\n";}
///////////////////////////////////////////////////////////////
    
//end main
	
		  Comm->Barrier();
	
#ifdef HAVE_MPI
	MPI_Finalize();
#endif
    return ierr;
}

// a function for nice output of time data
std::string time_out(double t)
{
	const double timesc = r0dim/udim; // time scale [s]
	std::string unit;  
	double val;
	// make sure the value is larger than one
	val = t*timesc/(365*24*3600);
	unit = "years";
	if (val<=1.0)
    {
		val*=12;
		unit="months";
		if (val<=1.0)
		{
			val*=(365.0/12.0);
			unit="days";
			if (val<=1.0)
			{
				val*=24;
				unit="h";
				if (val<=1.0)
				{
					val*=60;
					unit="min";
					if (val<=1.0)
					{
						val*=60;
						unit="s";
					}
				}
			}
		}
    }
	std::stringstream s;
	s<<val<<" ["<<unit<<"]";
	return s.str();
}
