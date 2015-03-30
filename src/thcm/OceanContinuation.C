/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
 
/************************************************************************
//OCEAN MODEL 
************************************************************************/

// define this to 
// - read starting solution
// - compute Jacobian A
// - build preconditioner P
// - compute eigenvalues of P\A-I
// - exit the program
//#define ANALYZE_SPECTRUM 1
#include <iostream>

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "LOCA_Epetra.H"
#include "LOCA.H"

#include "HYMLS_Tools.H"
#include "TRIOS_SolverFactory.H"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "HYMLS_MatrixUtils.H"
#include <sstream>

// global THCM definitions
#include "globdefs.H"

//User's application specific files 
#include "THCM.H"
#include "TRIOS_DefaultParams.H"
#include "DefaultParams.H"
#include "OceanModel.H"


#include "Epetra_LinearProblem.h"

// THCM is a singleton, there can be only one instance at a time. 
// As base class Singleton is templated we must instantiate it    
// using the macro defined in Singleton.H
_INSTANTIATE_(THCM)

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // this sets the paths in THCM, should be put somewhere else later
  //F90NAME(m_par,get_path)();

  
  int ierr = 0;   
  
////////////////////////////////////////
// Create the Parallel  Communicator  //
////////////////////////////////////////

#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  Teuchos::RCP<Epetra_SerialComm> Comm=rcp(new Epetra_SerialComm());
#endif

  //Get process ID and total number of processes
  int MyPID = Comm->MyPID();
  int NumProc = Comm->NumProc();

  ////////////////////////////////////////////////////////////////////////////
  // open static file streams (seperate for each process)                   //
  // Note: when using many processors it is convenient                      //
  // to make most of them dump their outupt to Teuchos::oblackholestreams   //
  // instead.                                                               //
  ////////////////////////////////////////////////////////////////////////////
  if (MyPID==0)
    {
    std::ostringstream infofile;  
    infofile << "info_"<<MyPID << ".txt";
    std::cout << "info is written to "<<infofile.str().c_str()<<std::endl;
    info = Teuchos::rcp(new std::ofstream(infofile.str().c_str()) );
    }
  else
    {
    info = Teuchos::rcp(new Teuchos::oblackholestream());
    }

(*info) << std::setw(15) << std::setprecision(15);

Teuchos::RCP<std::ostream> contstream;

contstream = Teuchos::rcp(new std::ofstream("continuation.out") );

(*contstream) << std::setw(15) << std::setprecision(15);

#ifdef DEBUGGING  
#if 1 // seperate files, set this manually here
  std::ostringstream debfile;
  debfile<<"debug_"<<MyPID<<".txt";
  std::cout << "debug output is written to "<<debfile.str().c_str()<<std::endl;
  debug=Teuchos::rcp(new std::ofstream(debfile.str().c_str()));
(*debug) << std::setw(15) << std::setprecision(15);
#else // dump everything into one file
  debug = info;
  contfile = info;
  std::cout << "debug output is written to "<<infofile.str().c_str()<<std::endl;
#endif
HYMLS::Tools::InitializeIO_std(Comm,info,debug);
#else
HYMLS::Tools::InitializeIO_std(Comm,info);
#endif 


  DEBUG("*********************************************")
  DEBUG("* Debugging output for process "<<MyPID)
  DEBUG("* To prevent this file from being written,  ")
  DEBUG("* omit the -DDEBUGGING flag when compiling. ")
  DEBUG("*********************************************")
  (*info) << *Comm << std::endl;
  Comm->Barrier();
  

  try {
  
  ////////////////////////////////////////////////////////
  // Setup Parameter Lists                              //
  ////////////////////////////////////////////////////////
  
    // read parameters from files
    // we create one big parameter list and further down we will set some things
    // that are not straight-forwars in XML (verbosity etc.)
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new 
    Teuchos::ParameterList);

	// modified_14_10_2014_erik
	// argument types in updateParametersFromXmlFile differ between versions
	// solution: .get() -> .ptr()
    DefaultParams::THCM(*paramList);
    Teuchos::updateParametersFromXmlFile("thcm_params.xml",paramList.ptr()); // modified_14_10_2014_erik
    
    DefaultParams::LOCA(*paramList);    
    Teuchos::updateParametersFromXmlFile("loca_params.xml",paramList.ptr()); // modified_14_10_2014_erik

    DefaultParams::NOX(*paramList);
    // create the NOX->Direction->Newton->Linear Solver sublist
    TRIOS::DefaultParams::LinearSolver(paramList->sublist("NOX")
                                                 .sublist("Direction")
                                                 .sublist("Newton"));
    // create the "BlockPreconditioner" sublist
    TRIOS::DefaultParams::BlockPreconditioner(paramList->sublist("NOX")
                                                 .sublist("Direction")
                                                 .sublist("Newton")
                                                 .sublist("Linear Solver"));
    // create the "HYMLS" sublist
    TRIOS::DefaultParams::HYMLS(paramList->sublist("NOX")
                                                 .sublist("Direction")
                                                 .sublist("Newton")
                                                 .sublist("Linear Solver"));
    // override default settings with parameters from user input file
    Teuchos::updateParametersFromXmlFile("solver_params.xml",paramList.ptr()); // modified_14_10_2014_erik

    // extract the final sublists:
            
    // Get the THCM sublist
    Teuchos::ParameterList& thcmList = paramList->sublist("THCM");

    // Get the LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");
    
    // get the Stepper sublist
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    std::string cont_param = stepperList.get("Continuation Parameter","Combined Forcing");
    double start_value = stepperList.get("Initial Value",0.0);

	// Get Anasazi Eigensolver sublist (needs --with-loca-anasazi)
//    Teuchos::ParameterList& aList = stepperList.sublist("Eigensolver");
//    aList.set("Verbosity", Anasazi::Errors+Anasazi::Warnings+Anasazi::FinalSummary);

    // Get the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("MyPID", MyPID);
    nlPrintParams.set("Output Stream",contstream);
    nlPrintParams.set("Error Stream",contstream);
    nlPrintParams.set("Output Process",0);
    nlPrintParams.set("Output Information",
			  NOX::Utils::Details + 
			  NOX::Utils::OuterIteration + 
			  NOX::Utils::InnerIteration +
                          NOX::Utils::OuterIterationStatusTest + 
                          NOX::Utils::LinearSolverDetails + 
                          NOX::Utils::Debug +
			  NOX::Utils::Warning +
			  NOX::Utils::StepperDetails +
			  NOX::Utils::StepperIteration +
			  NOX::Utils::StepperParameters); 


    //Create the "Direction" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");

    //Create the "Line Search" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");

    //Create the "Direction" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& newtParams = dirParams.sublist("Newton");

    //Create the "Linear Solver" sublist for the "Direction' sublist
    Teuchos::ParameterList& lsParams = newtParams.sublist("Linear Solver");

///////////////////////////////////////////////////////////
// Setup the Problem                                     //
///////////////////////////////////////////////////////////
        
    // put the correct starting value for the continuation parameter
    // in the thcm-list
    thcmList.sublist("Starting Parameters").set(cont_param,start_value);
        
    // Set up the THCM interface to allow calls to
    // residual (RHS) and Jacobian evaluation routines.
    // THCM is implemented as a Singleton, that means there
    // is only a single instance which should be referenced 
    // using THCM::Instance()
    Teuchos::RCP<THCM> ocean = Teuchos::rcp(new THCM(thcmList, Comm));

    // these LOCA data structures aren't really used by the OceanModel,
    // but are required for the ModelEvaluatorInterface

    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
    LOCA::createGlobalData(paramList, epetraFactory);

    // the model serves as Preconditioner factory if "User Defined" is selected
    std::string PrecType = lsParams.get("Preconditioner","Ifpack");
    
    Teuchos::RCP<Teuchos::ParameterList> myPrecList=Teuchos::null;
    if (PrecType=="User Defined")
      {
      myPrecList = Teuchos::rcp(&lsParams,false);
      }

    // for some purposes it's good to know which one is the continuation parameter
    // (i.e. backup in regular intervals)
    thcmList.set("Parameter Name",cont_param);                                     

    // this is the LOCA interface (LOCA::Epetra::Interface::TimeDependent) to
    // our EpetraExt ModelEvaluator class 'OceanModel'.
    Teuchos::RCP<OceanModel> model =
      Teuchos::rcp(new OceanModel(thcmList,globalData,myPrecList));

    // this vector defines what parameters the problem depends on
    // (or at least which parameters may be varied)
    // and gives initial values for them which overwrite
    // the settings made in usrc.F90::stpnt(). 
    Teuchos::RCP<LOCA::ParameterVector> pVector = 
                model->getParameterVector();

                                                     
    //Get the vector from the problem
    Teuchos::RCP<Epetra_Vector> soln = model->getSolution();

    //Initialize solution
    soln->PutScalar(0.0);

    // check for starting solution
    std::string StartConfigFile = thcmList.get("Starting Solution File","None");
                      
    if (StartConfigFile!="None")
      {
      THCM::Instance().startTiming("Read Start File");
      soln=model->ReadConfiguration(StartConfigFile,*pVector);
      try 
        {
        start_value = pVector->getValue(cont_param);
        } catch (...) {Error("Bad Continuation Parameter",__FILE__,__LINE__);}
      stepperList.set("Initial Value",start_value);
      model->setParameters(*pVector);
      THCM::Instance().stopTiming("Read Start File");
      }
    //Create the Epetra_RowMatrix for the Jacobian/Preconditioner
    Teuchos::RCP<Epetra_CrsMatrix> A = model->getJacobian();

    // get a pointer to the scaling arrays for the linear system
    Teuchos::RCP<NOX::Epetra::Scaling> scaling = THCM::Instance().getScaling();
#ifdef TEST_SOLVER_ON_THCM_MATRIX
    scaling = null;
#endif

// NOX/LOCA interface setup

    Teuchos::RCP<NOX::Abstract::PrePostOperator> prepost = model;
          
    // register pre- and post operations
    nlParams.sublist("Solver Options").set("User Defined Pre/Post Operator",prepost);
      
    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iReq = model;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = model;
  
    
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys = THCM::Instance().createLinearSystem
                (model, lsParams, nlPrintParams, contstream);                                    
    
    // register the linear system with the THCM session. This
    // is used to implement our own Preconditioner reuse policy
    // via the NOX PrePostOperator interface
    THCM::Instance().setLinSys(linsys);

    // we use the same linear system for the shift-inverted operator
    Teuchos::RCP<NOX::Epetra::LinearSystem> shiftedLinSys = linsys;

#ifdef ANALYZE_SPECTRUM
    INFO("EIGEN-TEST: compute A and B");
    model->computeShiftedMatrix(0.0,1.0,*soln,*A);
    Teuchos::RCP<Epetra_CrsMatrix> B = Teuchos::rcp(new Epetra_CrsMatrix(*A));
    model->computeJacobian(*soln,*A);
    
    // scale operator
    Epetra_LinearProblem linProb(A.get(),soln.get(),soln.get());
  
    INFO("EIGEN-TEST: scale jacobian...");  
    scaling->scaleLinearSystem(linProb);
  
    // create ocean preconditioner
  INFO("EIGEN-TEST: Build Preconditioner...");
  //CHECK_TRUE(linsys->recomputePreconditioner(*sol0,lsParams));
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> azlinsys = Teuchos::rcp_dynamic_cast
        <NOX::Epetra::LinearSystemAztecOO>(linsys);
  CHECK_TRUE(azlinsys->createPreconditioner(*soln,lsParams,false));
  Teuchos::RCP<const Epetra_Operator> P = azlinsys->getPrecOperator();
  Teuchos::ParameterList testList;
  Teuchos::RCP<const Epetra_Operator> A_op = A;
  testList.set("Eigen-Analysis: Operator",A_op);
  testList.set("Eigen-Analysis: Mass-Matrix",B);
  TRIOS::SolverFactory::AnalyzeSpectrum(testList,P);
  Error("Eigen-Analysis done, stopping here!",__FILE__,__LINE__);
#endif

    
#ifdef TEST_SOLVER_ON_THCM_MATRIX
// for testing: read THCM matrix/rhs, solve linear system to compare
test_with_thcm_linsys(linsys,soln,lsParams);
Error("Solver Test: stopping now!",__FILE__,__LINE__);
#endif    

    //Create the loca vector
    NOX::Epetra::Vector locaSoln(soln);
    
    // Create the Group
    Teuchos::RCP<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams,
                                           iReq, locaSoln, linsys, shiftedLinSys,
                                           *pVector));
    grp->computeF();

    double TolNewton = nlParams.get("Convergence Tolerance",1.0e-9);

    // Set up the Solver Convergence tests
    Teuchos::RCP<NOX::StatusTest::NormF> wrms =
       Teuchos::rcp(new NOX::StatusTest::NormF(TolNewton));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters
      = Teuchos::rcp(new NOX::StatusTest::MaxIters(searchParams.get("Max Iters", 10)));
    Teuchos::RCP<NOX::StatusTest::Combo> comboOR = 
       Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR)); 
    comboOR->addStatusTest(wrms);
    comboOR->addStatusTest(maxiters);

    // Create the stepper  
    LOCA::Stepper stepper(globalData, grp, comboOR, paramList);

#ifdef DEBUGGING    
    paramList->print(*debug);
#endif

paramList->print(*info);
//Teuchos::writeParameterListToXmlFile(*paramList,"parameters.xml");

INFO("\n*****************************");
INFO("Start Continuation process...");
INFO("*****************************\n\n");

THCM::Instance().startTiming("entire continuation run");

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();


THCM::Instance().stopTiming("entire continuation run",true);
   
    if (status != LOCA::Abstract::Iterator::Finished) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
        globalData->locaUtils->out()
          << "Stepper failed to converge!" << std::endl;
    }
    globalData->locaUtils->out() << "Continuation status -> \t"<< status << std::endl;

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out()
        << std::endl << "Final Parameters" << std::endl
        << "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RCP<const LOCA::Epetra::Group> finalGroup = 
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(stepper.getSolutionGroup());
    const NOX::Epetra::Vector& finalSolutionNOX = 
      dynamic_cast<const NOX::Epetra::Vector&>(finalGroup->getX());
    const Epetra_Vector& finalSolution = finalSolutionNOX.getEpetraVector(); 


    // write final backup file for restarting
    model->setForceBackup(true);
    // also write THCM files (fort.3, fort.44, fort.15)
    model->setThcmOutput(true);
    model->printSolution(finalSolution, model->getContinuationParameterValue());
    Comm->Barrier(); // make sure the root process has written the files
#if 0    
    MatrixUtils::DumpMatrix(*A, "FinalMatrix.txt");
#endif    
    LOCA::destroyGlobalData(globalData);

  } 
  catch (std::exception& e) {
    std::cerr << e.what() << std::endl; 
  }
  catch (const char*s) {
    std::cerr << s << std::endl; 
  }
  catch (...) {
    std::cerr << "Caught unknown exception!" << std::endl; 
  }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    
//end main
    return ierr;

}

