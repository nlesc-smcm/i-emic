#include "DefaultParams.H"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"


void DefaultParams::LOCA(Teuchos::ParameterList& list)
  {
  
  Teuchos::ParameterList& locaList = list.sublist("LOCA");


    /* We do not do bifurcation tracking yet, but here is the place 
       the parameters for it:                                       */
       
    //
    locaList.sublist("Bifurcation").set("Method","None");
    
    /* Parameters for the predictor method */
    
    // "Constant", "Secant" or "Tangent"
    locaList.sublist("Predictor").set("Method","Tangent");

    // "Constant" or "Tangent" (required if you choose "Secant" above)
    locaList.sublist("Predictor").sublist("First Step Predictor").set("Method","Tangent");
    
    /* Parameters influencing step size control in the continuation */

    // possible choices are "Constant" or "Adaptive"
    locaList.sublist("Step Size").set("Method","Adaptive");
    locaList.sublist("Step Size").set("Initial Step Size",5.0e-3);
    locaList.sublist("Step Size").set("Min Step Size",1.0e-5);
    locaList.sublist("Step Size").set("Max Step Size",5.0e-2);
    locaList.sublist("Step Size").set("Failed Step Reduction Factor",0.5);

    // This affects the "Adaptive" mode only 
    locaList.sublist("Step Size").set("Aggressiveness",0.5);

    // This affects the "Constant" mode only 
    locaList.sublist("Step Size").set("Successful Step Increase Factor",1.26);
    
    /* Parameters for the parameter advancing during the continuation */
    
    // "Natural" or "Arc Length" 
    locaList.sublist("Stepper").set("Continuation Method","Arc Length");

    // the "Continuation Parameter" has to be set by the user to one of the values
    // defined in class THCM (THCM.C), e.g. "Combined Forcing" or "Horizontal Ekman-Number"
    locaList.sublist("Stepper").set("Continuation Parameter","(Undefined)");

    locaList.sublist("Stepper").set("Initial Value",1.0e-3);
    locaList.sublist("Stepper").set("Min Value",0.0);
    locaList.sublist("Stepper").set("Max Value",1.0);

    locaList.sublist("Stepper").set("Max Steps",500);

    locaList.sublist("Stepper").set("Compute Eigenvalues",0);

    // I did not get the eigenvalue solver to run properly with our block solver,
    // so these settings are not reference values
    locaList.sublist("Stepper").sublist("Eigensolver").set("Method","Anasazi");
    // "Jacobian Inverse", "Cayley"
    
    //Cayley: (A-tau B)(A-sigma B)^{-1}  [I think...]
    locaList.sublist("Stepper").sublist("Eigensolver").set("Operator","Cayley");    
    // tau if you choose the Cayley operator
    locaList.sublist("Stepper").sublist("Eigensolver").set("Cayley Zero",-0.1);
    // sigma if you choose the Cayley operator
    locaList.sublist("Stepper").sublist("Eigensolver").set("Cayley Pole",0.0);
        
    // settings for the outer Krylov iteration 
        
    // block size for Block Krylov method 
    locaList.sublist("Stepper").sublist("Eigensolver").set("Block Size",1);
    locaList.sublist("Stepper").sublist("Eigensolver").set("Convergence Tolerance",1.0e-4);
    locaList.sublist("Stepper").sublist("Eigensolver").set("Max Restarts",3);
    //length of Krylov sequence
    locaList.sublist("Stepper").sublist("Eigensolver").set("Num Blocks",30);
    locaList.sublist("Stepper").sublist("Eigensolver").set("Num Eigenvalues",5);
    locaList.sublist("Stepper").sublist("Eigensolver").set("Num Blocks",30);
    // convergence check after this many steps
    locaList.sublist("Stepper").sublist("Eigensolver").set("Step Size",5);
    locaList.sublist("Stepper").sublist("Eigensolver").set("Linear Solver Tolerance",1.0e-5);
      
    // "Bordering" (default), "Householder" (prec reuse not working) 
    locaList.sublist("Stepper").set("Bordered Solver Method","Bordering");
    
    // there's a large number of advanced params here - arclength scaling and such - but
    // I leave these to LOCA to set default values.
      
  }



void DefaultParams::NOX(Teuchos::ParameterList& list)
  {
  
  Teuchos::ParameterList& noxList = list.sublist("NOX");
  
  // "Line Search Based" or "Trust Region Based" (not tested)
  noxList.set("Nonlinear Solver","Line Search Based");
  // convergence tolerance for the Newton process
  noxList.set("Convergence Tolerance",1.0e-6);

  /* Line Search parameters  */

  // "Full Step" is standard Newton's, "Backtrack" is a simple line search
  // (more robust)
  noxList.sublist("Line Search").set("Method","Backtrack");
  
  noxList.sublist("Line Search").sublist("Backtrack").set("Default Step",1.0);
  noxList.sublist("Line Search").sublist("Backtrack").set("Max Iters",10);
  noxList.sublist("Line Search").sublist("Backtrack").set("Minimum Step",0.01);
  noxList.sublist("Line Search").sublist("Backtrack").set("Recovery Step",0.001);
  noxList.sublist("Line Search").sublist("Backtrack").set("Reduction Factor",0.5);
      
  // maximum number of Newton iterations 
  noxList.sublist("Line Search").set("Max Iters",10);
      
  // we leave the "Backtrack" sublist to NOX for default values
    
  /* Parameters for the actual Newton solver */  
  noxList.sublist("Direction").set("Method","Newton");
       
        // Method to determine convergence tolerance for linear solver 
        // "Constant", "Type 1" (fails!) or "Type 2"
        noxList.sublist("Direction").sublist("Newton").set("Forcing Term Method","Type 2");
        noxList.sublist("Direction").sublist("Newton").set("Forcing Term Initial Tolerance",1.0e-2);
        noxList.sublist("Direction").sublist("Newton").set("Forcing Term Maximum Tolerance",1.0e-2);
        noxList.sublist("Direction").sublist("Newton").set("Forcing Term Minimum Tolerance",1.0e-5);        
        
        // Use Newton step even if linear solve failed (default 1)
        // If you say "0" here the continuation run _stops_ when  
        // the linear solver fails to achieve the requested tol   
        noxList.sublist("Direction").sublist("Newton").set("Rescue Bad Newton Solve",true);

  }
        

void DefaultParams::THCM(Teuchos::ParameterList& list)
  {

//                                                            
// This test case is an idealized single-hemispheric basin    
// at 0.5 degree resolution (128x128x16), idealized restoring 
// T and S forcing, wind from data. The solution was obtained 
// using the 'old' linear equation of state (alphaT=1.0).     
// compile the code with -DOLD_EOS for both fortran and C++   
// files.                                                     
// Ah has the value 10^5, it is still hard-coded in usr.F90.  
//                                                            

  Teuchos::ParameterList& thcmList = list.sublist("THCM");

    // a descriptive name to identify the settings 
    thcmList.set("Problem Description","not specified");

    // we don't really specify a default problem here, just some settings that
    // require reasonable defaults

    // read starting guess from backup file XYZ. set this to "None" to start from 0 
    // in case this is given, both the parameter values and the solution vector are 
    // read from a file as generated by the "Backup Interval" option                
    thcmList.set("Starting Solution File","None");

    // store approximation at regular intervals (set to -X to disable)   
    // (the value refers to the continuation parameter or dim.less time) 
    thcmList.set("Backup Interval",0.1);

    // this can be used to write a series file in Xdmf (-X to disable).  
    // This works for both time-integration and continuation             
    // (the value refers to the continuation parameter or dim.less time) 
    thcmList.set("Output Frequency",-1.0);
    
    // Type of scaling applied to the linear systems.                            
    // We currently support "None", "THCM" and "THCM/TS"                         
    // "None" is not really recomended.                                          
    // "THCM" is the block scaling used in THCM                                  
    // "THCM/TS" is like THCM with additional row scaling of the T/S equations   
    thcmList.set("Scaling","THCM");
    
    // this can be used to have a new preconditioner computed after k time- or    
    // continuation steps, or after k Newton solves. You probably want to set    
    // the "Max Age of Prec" option                                              
    // in solver_params.xml to -2 and the "Reuse Policy" to "Reuse". The reason  
    // why we implement our own reuse strategy is that NOX tends to rebuild the  
    // Preconditioner in the middle of a Newton-solve.                           
    // If you do not want to reuse preconditioner, use the according options in  
    // solver_params.xml ("Preconditioner Reuse Policy" = "Rebuild").            
    // If you want to reuse the preconditioner indefinitely or control the age   
    // by the according option in solver_params.xml, set this one to -2.         
    thcmList.set("Max Age Of Prec",-2);

    // this can be - "None" (leave this to NOX, see solver_params.xml)           
    //             - "Newton Solve" (recompute after k Newton solves)            
    //             - "Step" (recompute after k time-/continuation steps)         
    // here k is set by "Max Age Of Prec" above.                                 
    thcmList.set("Preconditioner Reuse Policy","None");

    // this is only done if you compile with -DHAVE_XDMF (requires -lhdf5) and
    // set "Output Interval" to a postive value.
    
      // currently only version 2.0 is implemented 
      thcmList.sublist("Xdmf Output").set("Xdmf Version","2.0");
      // filename for the output (.xmf and .h5 are appended) 
      thcmList.sublist("Xdmf Output").set("File Name","thcm");
      // write coordinates in real x/y/z positions rather than phi/theta/z 
      thcmList.sublist("Xdmf Output").set("Store X/Y/Z Mesh",true);

      // If you choose "Store X/Y/Z Mesh", this determines the radius of the 
      // earth for visualization. The ocean has depth 1.                     
      thcmList.sublist("Xdmf Output").set("Earth-Radius for Visualization",10.0);

      // select which components are to be stored 
      thcmList.sublist("Xdmf Output").set("Store Velocity",true);
      thcmList.sublist("Xdmf Output").set("Store Pressure",true);
      thcmList.sublist("Xdmf Output").set("Store Temperature",true);
      thcmList.sublist("Xdmf Output").set("Store Salinity",true);
    
  }

void DefaultParams::TimeStepping(Teuchos::ParameterList& list)
  {

  Teuchos::ParameterList& timeList = list.sublist("Time Stepping");
  
    // choose the time-stepping scheme.      
    // presently we support                  
    //                                       
    // explicit schemes:                     
    //            "Forward Euler"            
    //                                       
    // implicit schemes:                     
    //            "Theta", "Adaptive Theta"  
    //                                       
    // the latter chooses theta smartly (?)  
    // as it was implemented in THCM         
    timeList.set("Scheme","Adaptive Theta");
    
    // only relevant for "Theta" scheme 
    timeList.set("Theta","0.5");

    // time is measured in dimensionless units 
    // 1 unit is about 2 years (R/u0)          
    // 0.04 would be (roughly) a month,        
    // 0.0014 is a day                         

    timeList.set("Start Time",0.0);

    timeList.set("End Time",1000.0);

    timeList.set("Step Size",0.12);

    // currently the step counter is not incremented when a step fails, 
    timeList.set("Max Num Steps",5000);

    // the predictor can currently be                         
    // "Constant": start Newton with the previous slution u_n 
    // "Secant":   start with u_n + dt*(u_n-u_{n-1})/dt_n     
    // In the first time-step or after a failed step, it is   
    // always "Constant".                                     
    timeList.set("Predictor","Secant");
    
    // make sure the stepper satisfies the CFL condition. 
    // for implicit time-stepping, you can choose CN>1    
    // to disable CFL checks, set CN<0.                   
    timeList.set("Courant Number",50.0);

    // this indicates when the time-stepper considers the flow slow 
    // (i.e. when it is at factor*CN*cfl-number)                    
    timeList.set("Slow Flow Condition",0.4);
    
    // after nslw steps of slow flow we allow the step-size to be increased 
    timeList.set("Slow Flow Delay",3);
    
    // slowly increase a parameter at the beginning of the simulation 
    timeList.set("Gradual Start-Up",false);
    // which parameter to turn on slowly? 
    timeList.set("Start-Up Parameter","Wind Forcing");
    // by how much should the parameter be increased per unit time? 
    timeList.set("Start-Up Rate",0.05);
       
    // "Constant" or "Adaptive"                                                 
    // if "Constant", the time step is changed only based on success or         
    // failure of the Newton solve. "Adaptive" adjusts the time step based on   
    // an estimate of the local truncation error.                               
    // In case there is a stability condition, it is applied after step-size    
    // control based on error estimation                                        
    timeList.set("Step Size Control","Constant");
    
    // this causes the program to terminate if the step-size becomes too small 
    timeList.set("Minimum Step Size",1.0e-4);
    
    // puts an absolute bound on the step size           
    timeList.set("Maximum Step Size",10.0);
    
    // this is only for the "Constant" method above 
    timeList.set("Failed Step Reduction Factor",0.5);

    // this is only for the "Constant" method above 
    timeList.set("Successful Step Increase Factor",1.25);
    
    // this is for the "Adaptive" mode only                    
    timeList.set("Error Tolerance",1.0);

    // stop the time-integration when a steady-state is detected. the following 
    // criteria can be used:                                                    
    // "None": no convergence monitoring                                       
    // "Norm of F": ||F(u)||_1<tol                                              
    // "Overturning": monitor the meridional streamfunction, i.e.               
    //                the maximum of the y-z streamfunction (v is averaged over x) 
    timeList.set("Steady-State Monitor","Overturning");

    // criterion to stop if the monitor chosen above goes below a given  
    // tolerance, indicating that a steady state has been found and      
    // further time integration is of little use.                        
    timeList.set("Steady-State Tolerance",-1.0);
     
}

void DefaultParams::PeriodicOrbits(Teuchos::ParameterList& list)
  {
  // enable monthly forcing
  list.sublist("THCM").set("Time Dependent Forcing",true);
  // one year
  list.sublist("THCM").set("Time Period",0.4951);
  //
  list.sublist("THCM").sublist("Starting Parameters").set("Seasonal Forcing",1.0);

   // Global Preconditioner for computing periodic orbits  
   // can be either "Circulant" (our own method) or one of 
   // "BlockDiagonal", "Sequential", "Parallel" etc        
   // (see LOCA_Epetra_xyztPrec class)                     
   // This option is only relevant for the 'periodic'      
   // driver                                               
   list.sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver")
         .set("XYZTPreconditioner","BlockDiagonal");

   // do inner GMRESR (1) or just apply prec to blocks (0) 
   list.sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver")
         .set("XYZT Inner Solve",true);
   
   // (time-)periodic problem (for XYZT preconditioner) 
   list.sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver")
         .set("Periodic",true);

  
  }
