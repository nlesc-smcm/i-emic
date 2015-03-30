/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
 
/************************************************************************
//OCEAN MODEL 
************************************************************************/
#include "iostream"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

#include "AztecOO.h"

#include "LOCA.H"
#include "LOCA_Epetra.H"

#include "TRIOS_BlockPreconditioner.H"
#include "TRIOS_SolverFactory.H"
#include "TRIOS_Domain.H"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_MultiComm.h"
#include "EpetraExt_MultiMpiComm.h"

#include "NOX_Epetra_Scaling.H"

#include "NOX_Epetra_LinearSystem_Belos.H"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#ifdef HAVE_XDMF
#include "EpetraExt_HDF5.h"
#endif
#include "HYMLS_MatrixUtils.H"
#include "HYMLS_Tools.H"
#include <sstream>

// global THCM definitions
#include "globdefs.H"


#include "TRIOS_Epetra_Interface_xyzt.H"

//User's application specific files 
#include "THCM.H"
#include "DefaultParams.H"
#include "TRIOS_DefaultParams.H"
#include "OceanModel.H"


typedef NOX::Epetra::LinearSystemBelos LinSysType;

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
// Create the Global  Communicator    //
////////////////////////////////////////

#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> CommWorld = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  Teuchos::RCP<Epetra_SerialComm> CommWorld=Teuchos::rcp(new Epetra_SerialComm());
#endif

int MyGlobalPID = CommWorld->MyPID();
int GlobalNumProc=CommWorld->NumProc();

if (argc==1) // print usage and exit
  {
  if (MyGlobalPID==0)
    {
    std::cout <<"USAGE: mpirun -np N ./periodic Nx nt\n";
    std::cout <<"       where                           \n";
    std::cout <<"             Nx: number of procs sharing spatial domain\n";
    std::cout <<"                 (defaults to 1)                       \n";
    std::cout <<"             nt: total number of timesteps             \n";
    std::cout <<"                 (defaults to 1)                       \n";
    }
  MPI_Finalize();
  return 0;
  }

  int spatialProcs=1;
  if (argc>1) { spatialProcs = atoi(argv[1]);}
  int numTimeSteps= 1; // default
  if (argc>2) { numTimeSteps = atoi(argv[2]);}

  Teuchos::RCP<EpetraExt::MultiComm> globalComm =
    Teuchos::rcp(new EpetraExt::MultiMpiComm(MPI_COMM_WORLD, spatialProcs, numTimeSteps));
  Teuchos::RCP<Epetra_Comm> _Comm = 
        Teuchos::rcp(&(globalComm->SubDomainComm()),false);
  Teuchos::RCP<Epetra_MpiComm> Comm = 
        Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(_Comm);

  //Get process ID and total number of processes
  int MyPID = Comm->MyPID();
  int NumProc = Comm->NumProc();
  int MyTimeDomain = globalComm->SubDomainRank();

  ////////////////////////////////////////////////////////////////////////////
  // open static file streams (seperate for each process)                   //
  // Note: when using many processors it is convenient                      //
  // to make most of them dump their outupt to Teuchos::oblackholestreams   //
  // instead.                                                               //
  ////////////////////////////////////////////////////////////////////////////
  if (MyPID==0)
    {
    std::ostringstream infofile;
    infofile << "info_"<<MyPID<<"_"<<MyTimeDomain << ".txt";
    std::cout << "info is written to "<<infofile.str().c_str()<<std::endl;
    info = Teuchos::rcp(new std::ofstream(infofile.str().c_str()) );
    }
  else
    {
    info = Teuchos::rcp(new Teuchos::oblackholestream());
    }

(*info) << std::setw(15) << std::setprecision(15);

Teuchos::RCP<std::ostream> contstream;

if (MyGlobalPID==0)
  {
  contstream = Teuchos::rcp(new std::ofstream("periodic.out") );
  (*contstream) << std::setw(15) << std::setprecision(15);
  }
else
  {
  contstream = Teuchos::rcp(new Teuchos::oblackholestream());
  }

//contstream=info;

#ifdef DEBUGGING
#if 1 // seperate files, set this manually here
  std::ostringstream debfile;
  debfile<<"debug_"<<MyGlobalPID<<".txt";
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
  
//(*info) << *globalComm << std::endl;
//(*info) << *Comm << std::endl;
  Comm->Barrier();
  

  try {
  
  ////////////////////////////////////////////////////////
  // Setup Parameter Lists                              //
  ////////////////////////////////////////////////////////
  
    // read parameters from files
    // we create one big parameter list and further down we will set some things
    // that are not straight-forwars in XML (verbosity etc.)
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new Teuchos::ParameterList);

	// modified_14_10_2014_erik
	// argument types in updateParametersFromXmlFile differ between versions
	// solution: .get() -> .ptr()
    DefaultParams::THCM(*paramList);
    Teuchos::updateParametersFromXmlFile("thcm_params.xml",paramList.ptr()); 	// modified_14_10_2014_erik
    DefaultParams::LOCA(*paramList);
    Teuchos::updateParametersFromXmlFile("loca_params.xml",paramList.ptr());    // modified_14_10_2014_erik
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

    DefaultParams::PeriodicOrbits(*paramList);
    Teuchos::updateParametersFromXmlFile("solver_params.xml",paramList.ptr()); // modified_14_10_2014_erik
            
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
                          //NOX::Utils::Debug +
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

    // this vector defines what parameters the problem depends on
    // (or at least which parameters may be varied)
    // and gives initial values for them which overwrite
    // the settings made in usrc.F90::stpnt(). 
    Teuchos::RCP<LOCA::ParameterVector> pVector;
    pVector=THCM::Instance().getParameterVector();
   std::string PrecType = lsParams.get("Preconditioner","User Defined");

    thcmList.set("Parameter Name",cont_param);                                     

    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);
                
    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, epetraFactory);

    //Create the model interface
    Teuchos::RCP<OceanModel> model = 
      Teuchos::rcp(new OceanModel(thcmList,globalData));

    //Get the vector from the problem
    Teuchos::RCP<Epetra_Vector> soln = model->getSolution();

    //Initialize solution
    soln->PutScalar(0.0);

    // check for starting solution
   std::string StartConfigFile = thcmList.get("Starting Solution File","None");

    Epetra_MultiVector initGuess(soln->Map(), globalComm->NumTimeStepsOnDomain());

  if (StartConfigFile!="None")
    {
    THCM::Instance().startTiming("Read Start File");

    if (StartConfigFile=="Config*.txt")
      {
      for (int i=0; i<globalComm->NumTimeStepsOnDomain(); i++)
        {
        int id = globalComm->FirstTimeStepOnDomain()+i;
        std::stringstream ss;
        ss<<"Config"<<id<<".txt";
       std::string levelStartConfigFile=ss.str();
        soln=model->ReadConfiguration(levelStartConfigFile,*pVector);
        *(initGuess(i)) = *soln;
        }
//      }
//    else if (StartConfigFile=="solution.hdf")
   
//      {
#if 0
//#ifdef HAVE_XDMF
  bool verbose=true;
  bool success;
  Teuchos::RCP<EpetraExt::HDF5> hdf5 = Teuchos::rcp(new EpetraExt::HDF5(*globalComm));
    hdf5->Open("solution.hdf");
    hdf5->Read("vector",initGuess);
    hdf5->Close();
#else
  //Error"Can't read HDF5 file! Recompile with -DHAVE_XDMF.",__FILE__,__LINE__);
#endif      
      }
    else 
      {
      soln=model->ReadConfiguration(StartConfigFile,*pVector);
      for (int i=0; i<globalComm->NumTimeStepsOnDomain(); i++)
        {
        *(initGuess(i)) = *soln;
        }
      }
    try 
      {
      start_value = pVector->getValue(cont_param);
      } catch (...) {Error("Bad Continuation Parameter",__FILE__,__LINE__);}
    stepperList.set("Initial Value",start_value);
    model->setParameters(*pVector);
    THCM::Instance().stopTiming("Read Start File");
    }//restart

    //Create the Epetra_RowMatrix for the Jacobian/Preconditioner
    Teuchos::RCP<Epetra_CrsMatrix> A = model->getJacobian();
    
    /*    Epetra_Vector F(*soln);
    model->computeF(*soln,F,OceanModel::Residual);
    model->computeShiftedMatrix(0.0,-1.0,*soln,*A);
    INFO("B-MATRIX:");
    MatrixUtils::DumpMatrix(*A,"MassMatrix.txt");
    model->computeJacobian(*soln,*A);
    INFO("A-MATRIX:");
    MatrixUtils::DumpMatrix(*A,"StiffnessMatrix.txt");
    F.PutScalar(0.0);
    model->setSigma(F,-1.0);
    model->computeJacobian(*soln,*A);
    INFO("A-B MATRIX:");
    MatrixUtils::DumpMatrix(*A,"AminusB.txt");
*/

    Teuchos::RCP<Teuchos::ParameterList> precLSParams = Teuchos::rcp(new Teuchos::ParameterList(lsParams));
  Teuchos::RCP<Teuchos::ParameterList> precPrintParams =
       Teuchos::rcp(new Teuchos::ParameterList(nlPrintParams));

    // get output behaviour right
    // note: the 'outer streams' are in fact on the intermediate level of the diagonal
    // blocks here. Since xyztPrecond doesn't care about Aztec's bad output habits,
    // both the Jacobian and diagonal block solves will be printed to contstream.

  precPrintParams->set("Output Stream", info);
  precPrintParams->set("Error Stream", info);

    Teuchos::RCP<Teuchos::FancyOStream> fancyOuter = Teuchos::rcp
        (new Teuchos::FancyOStream(contstream), false);
    Teuchos::RCP<Teuchos::FancyOStream> fancyInner = Teuchos::rcp
        (new Teuchos::FancyOStream(info), false);
    precLSParams->sublist("Block Preconditioner").set("Outer Output Stream",fancyOuter);
    precLSParams->sublist("Block Preconditioner").set("Outer Error Stream",fancyOuter);
    precLSParams->sublist("Block Preconditioner").set("Inner Output Stream",fancyInner);
    precLSParams->sublist("Block Preconditioner").set("Inner Output Stream",fancyInner);
       
  double T = thcmList.get("Time Period",1.0);
  double dt = T/numTimeSteps;
  
  Teuchos::RCP<const Epetra_MultiVector> nullSpace = 
        model->getNullSpace();
        
  precLSParams->set("Null Space",nullSpace);
  
  bool innerSolver = precLSParams->get("XYZT Inner Solve",true);
  precLSParams->set("Use Preconditioner as Solver",!innerSolver);
  
  Teuchos::RCP<TRIOS::Epetra::Interface::xyzt> ixyzt = 
              Teuchos::rcp(new TRIOS::Epetra::Interface::xyzt(model,
                               initGuess, A,
                               globalComm,
                               *soln,dt,
                               precPrintParams.get(),
                               precLSParams.get()));

  // this object defines left and right scaling for a single time level
  Teuchos::RCP<NOX::Epetra::Scaling> scaling3D = 
    THCM::Instance().getScaling();
  // this object defines left and right scaling for the whole 4D matrix
  Teuchos::RCP<NOX::Epetra::Scaling> scaling4D = 
    THCM::Instance().getXyztScaling(ixyzt->getSolution().Map(), globalComm);

  Teuchos::RCP<EpetraExt::BlockCrsMatrix> Axyzt =
     Teuchos::rcp(&(ixyzt->getJacobian()),false);

  Epetra_Vector& solnxyzt = ixyzt->getSolution();
  
  Teuchos::RCP<Epetra_Operator> Mxyzt =
     Teuchos::rcp(&(ixyzt->getPreconditioner()),false);

Teuchos::RCP<TRIOS::Epetra::xyztPrec> blockprec = 
  Teuchos::rcp_dynamic_cast<TRIOS::Epetra::xyztPrec>(Mxyzt);
  if (blockprec!=Teuchos::null)
    {
    INFO("set scaling for diagonal blocks...");
    // this is needed by the Circulant solver
    blockprec->resetScaling(scaling3D);
    }

  if (PrecType=="User Defined")
    {
   std::string UserDefPrecType=precLSParams->get("User Defined Preconditioner","Block Preconditioner");
   std::string xyztPrecType=precLSParams->get("XYZTPreconditioner","BlockDiagonal");
    
    if (UserDefPrecType!="Block Preconditioner")
      {
      Error("Only 'Block Preconditioner' is allowed as 'User Defined Preconditioner' in 4D",
        __FILE__,__LINE__);
      }
      
    if (xyztPrecType=="Circulant")
      {
      Error("'User Defined' Preconditioner doesn't work with 'Circulant' XYZTPreconditioner yet!",
        __FILE__,__LINE__);
      }
    
    TRIOS::BlockPreconditioner *dummyPrec=NULL;
    if (blockprec!=Teuchos::null)
      {
      CHECK_ZERO(blockprec->setUserDefinedPreconditioner(dummyPrec));
      }
    }


  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = ixyzt;

  // Create the Linear System
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = ixyzt;

  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec =
    Teuchos::rcp(&(ixyzt->getPreconditioner()),false);

  Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
    
  DEBUG("Create global Linear System...");

/*
if (scaling4D==null)
  {
  INFO("Warning: no scaling is used, this is typically a bad idea with THCM!");
  INFO("file "<<__FILE__<<", line "<<__LINE__<<")");
  }
*/

if (Mxyzt==Teuchos::null)
  {
  INFO("Create linear system without user-defined preconditioner");
  //Create the linear systems without Ocean preconditioner
  linsys = Teuchos::rcp(new LinSysType(nlPrintParams,
                   lsParams, iReq,iJac, Axyzt, solnxyzt, scaling4D));
  }
else 
  {
  INFO("Create linear system with user-defined preconditioner");
  lsParams.set("Preconditioner","User Defined");
  //Create the linear systems with Ocean preconditioner
  linsys = Teuchos::rcp(new LinSysType(nlPrintParams,
                   lsParams, iJac, Axyzt, iPrec, Mxyzt, solnxyzt, scaling4D));
  }

    // register the linear system with the THCM session. This
    // is used to implement our own Preconditioner reuse policy
    // via the NOX PrePostOperator interface
//    THCM::Instance().setLinSys(linsys);

    // we use the same linear system for the shift-inverted operator
    //Teuchos::RCP<NOX::Epetra::LinearSystem> shiftedLinSys = linsys;

    //Create the loca vector
    NOX::Epetra::Vector locaSoln(solnxyzt);
    
    // Create the Group
  DEBUG("Create LOCA group...");
    Teuchos::RCP<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams,
                                           iReq, locaSoln, linsys,
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
  DEBUG("Create LOCA stepper...");
    LOCA::Stepper stepper(globalData, grp, comboOR, paramList);


paramList->print(*info);
//Teuchos::writeTeuchos::ParameterListToXmlFile(*paramList,"parameters.xml");

INFO("\n*****************************");
INFO("Start Continuation process...");
INFO("*****************************\n\n");

THCM::Instance().startTiming("entire continuation run");


  // Initialize time integration parameters
  int maxTimeSteps = 1; // No longer need a time integration loop
  int timeStep = 0;
  double time = 0.;

  // Time integration loop
//  while(timeStep < maxTimeSteps) {

    timeStep++;
    time += dt;

    globalData->locaUtils->out()
      << "Time Step: " << timeStep << ",\tTime: " << time << std::endl; 

//    NOX::StatusTest::StatusType status = solver.solve();
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status == LOCA::Abstract::Iterator::Finished)
      globalData->locaUtils->out() << "All tests passed" << std::endl; 
    else
       globalData->locaUtils->out() << "Stepper failed to converge!" << std::endl; 


    // Get the Epetra_Vector with the final solution from the solver
    const LOCA::Epetra::Group& finalGroup =
      dynamic_cast<const LOCA::Epetra::Group&>(*(stepper.getSolutionGroup()));

    const Epetra_Vector& finalSolution =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // End Nonlinear Solver **************************************

  //} // end time step while loop

  // Output the parameter list
  if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
      globalData->locaUtils->out()
        << std::endl << "Final Parameters" << std::endl
        << "****************" << std::endl; 
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl; 
    }
                                                         

model->setForceBackup(true);
//model->setThcmOutput(true); //I think this doesn't work properly, yet
ixyzt->printSolution(finalSolution,model->getContinuationParameterValue());

#ifdef HAVE_XDMF
// also store it in a single vector:
  Teuchos::RCP<EpetraExt::HDF5> hdf5 = Teuchos::rcp(new 
  EpetraExt::HDF5(finalSolution.Comm()));
  bool new_file=true;
  if (new_file)
    {
    hdf5->Create("solution.hdf");
    }
  else
    {
    hdf5->Open("solution.hdf");
    }

  hdf5->Write("vector",finalSolution);
  hdf5->Close();
#endif

#if 0
// construct the mean of January minus the mean of July (see Marotzke paper)
*contstream << "Compute and save january minus july solution..."<<std::endl; 
try {
  solnxyzt = finalSolution;
  EpetraExt::BlockVector& all_year_sol =
        dynamic_cast<EpetraExt::BlockVector&>(solnxyzt);
    Epetra_Vector january_sol(*soln);
    Epetra_Vector july_sol(*soln);
    Epetra_Vector partial_sum(*soln);
    Epetra_Vector seasonal_sol(*soln);//january - july
    
    january_sol.PutScalar(0.0);
    july_sol.PutScalar(0.0);
    seasonal_sol.PutScalar(0.0);
    partial_sum.PutScalar(0.0);
    
    int first_month = globalComm->FirstTimeStepOnDomain();
    int last_month = first_month + globalComm->NumTimeSteps();
    
    int january = 0; int july = 6;
    if (first_month<=january && last_month>=january)
      {
      CHECK_ZERO(all_year_sol.ExtractBlockValues(january_sol,january));
      partial_sum = january_sol;
      }
    if (first_month<=july && last_month>=july)
      {
      CHECK_ZERO(all_year_sol.ExtractBlockValues(july_sol,july));
      CHECK_ZERO(partial_sum.Update(-1.0,july_sol,1.0));
      }
    double *partial_result, *result;
    CHECK_ZERO(partial_sum.ExtractView(&partial_result));
    CHECK_ZERO(seasonal_sol.ExtractView(&result));
    int ndim = seasonal_sol.MyLength();
    CHECK_ZERO(globalComm->SumAll(partial_result,result,ndim));    

  if (first_month==january)
    {
    // TODO: this interferes with our backup and normal THCM output,
    //       we have to work out how to do this seperately.
    model->dataForPrintSolution(0,0,0);
    model->setForceBackup(true);
    model->setThcmOutput(true);
    model->printSolution(seasonal_sol,0.0);
    }

    } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,ierr)

#endif

THCM::Instance().stopTiming("entire continuation run",true);
   
    Comm->Barrier(); // make sure the root process has written the files
    
    LOCA::destroyGlobalData(globalData);

  } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,ierr)
/*  
catch (std::exception& e) {
    std::cerr << e.what() << std::endl; 
  }
  catch (const char*s) {
    std::cerr << s << std::endl; 
  }
  catch (...) {
    std::cerr << "Caught unknown exception!" << std::endl; 
  }
*/
#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    
//end main
    return ierr;

}

