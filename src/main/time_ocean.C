//=======================================================================
// Main continuation of the ocean model
//=======================================================================

#include "RunDefinitions.H"
#include "OceanTheta.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

using JDQZsolver = JDQZ<JDQZInterface<Teuchos::RCP<Ocean>,
                                      ComplexVector<Epetra_Vector> > >;

//------------------------------------------------------------------
void runOceanModel(RCP<Epetra_Comm> Comm);
void writeData(Teuchos::RCP<OceanTheta> model, bool describe,
               int step, double dt, double time, int niters);

//------------------------------------------------------------------
int main(int argc, char **argv)
{
	// Initialize the environment:
	//  - MPI
	//  - output files
	//  - returns Trilinos' communicator Epetra_Comm
	RCP<Epetra_Comm> Comm = initializeEnvironment(argc, argv);

	// run the ocean model
	runOceanModel(Comm);

    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	MPI_Finalize();
}

//------------------------------------------------------------------
void runOceanModel(RCP<Epetra_Comm> Comm)
{
	TIMER_START("Total time...");

	//------------------------------------------------------------------
	// Check if outFile is specified
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");

	// Create parameter object for Ocean
	RCP<Teuchos::ParameterList> oceanParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
    oceanParams->setName("Ocean parameters");

	// Create parallelized OceanTheta object
	RCP<OceanTheta> oceanTheta = Teuchos::rcp(new OceanTheta(Comm, oceanParams));

    // Create parameter object for time stepping
	RCP<Teuchos::ParameterList> timeParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("timestepper_params.xml", timeParams.ptr());
    timeParams->setName("time stepper parameters");

    // get parameters from xml
    double dt     = timeParams->get("initial time step size", 1.0e-03);
    double mindt  = timeParams->get("minimum step size", 1.0e-8);
    double maxdt  =  timeParams->get("maximum step size", 1.0);
    double iscale = timeParams->get("increase step size", 2.0);
    double dscale = timeParams->get("decrease step size", 2.0);

    int minK = timeParams->get("minimum desired Newton iterations", 3);
    int maxK = timeParams->get("maximum desired Newton iterations", 3);

    double tend   = timeParams->get("end time (in y)", 1.0);
    double theta  = timeParams->get("theta", 1.0);
    int nsteps    = timeParams->get("number of time steps", 10);
    int output    = timeParams->get("HDF5 output frequency", 1);
    double Ntol   = timeParams->get("Newton tolerance", 1e-6);
    int Niters    = timeParams->get("maximum Newton iterations", 10);

    // conversion from nondimensional time to years
    double years  = 2.01991; // r0dim / udim / 3600 / 24 / 365
    double time   = 0.0;     // keep track of time
    double normF  = 0.0;     // keep track of residual norm
    double normdx = 0.0;     // keep track of update infnorm
    int    step   = 0;       // keep track of time steps

    Teuchos::RCP<Epetra_Vector> F    = oceanTheta->getRHS('V');
    Teuchos::RCP<Epetra_Vector> x    = oceanTheta->getState('V');
    Teuchos::RCP<Epetra_Vector> dx   = oceanTheta->getSolution('V');
    Teuchos::RCP<Epetra_Vector> xdot = oceanTheta->getState('C');
    Teuchos::RCP<Epetra_Vector> xold = oceanTheta->getState('C');
    dx->PutScalar(0.0);
    xdot->PutScalar(0.0);
    xold->PutScalar(0.0);
    
    oceanTheta->initializeSolver();
    writeData(oceanTheta, true, 0, 0, 0, 0); // write data to tdata.txt
    oceanTheta->setTheta(theta);

    bool test_step = ( nsteps < 0 ) ? true : step < nsteps;
    while ( ( time < tend ) ||
            ( test_step ) )
    {
        INFO("----------------------------------------------------------");
        INFO("Timestepping:    t = " << time * years << " y");        

        oceanTheta->preProcess();        
        oceanTheta->store();  // save current state to oldstate

        // compute time discretization -G:
        // -G(x) =  -1 * (B d/dt x + theta*F(x) + (theta-1) * F(x_old))
        oceanTheta->setTimestep(dt);
        oceanTheta->computeRHS();
        normF = Utils::norm(F);
        
        // Timestep adjustments
        // while (normF > 10)
        // {
        //     dt = std::max(dt / dscale, mindt);
        //     oceanTheta->setTimestep(dt);
        //     oceanTheta->computeRHS();
        //     normF = Utils::norm(F);
        //     INFO("             ||F|| = " << normF);
        // }

        INFO("              step = " << step);
        INFO("              Newton solve, ||F|| = " << normF);

        // // Secant prediction (not sure if needed)
        // x->Update(dt, *xdot, 1.0);
        // *xold = *x;
        
        // Newton solve
        int k = 0;
        for (; k != Niters; ++k)
        { 
            // create jacobian of time discretization
            oceanTheta->computeJacobian();

            // solve J(x) dx = -G(x)
            oceanTheta->solve(F);

            // correct for pressure modes (superfluous due to prec)
            oceanTheta->pressureProjection(dx);

            // update state
            x->Update(1.0, *dx, 1.0);

            // compute time discretization -G(x)
            oceanTheta->computeRHS();
            
            normF  = Utils::norm(F);
            normdx = Utils::normInf(dx);
            INFO("                            ||F||2   = " << normF);
            INFO("                           ||dx||inf = " << normdx);
            if ( normdx < Ntol )
                break;
            else if ( normdx > 1e2)
            {
                k = Niters;
                break;
            }
        }
        
        if (k == Niters)
        {
            WARNING("Newton did not converge! ||F|| = "
                    << normF << "\nRestoring model",
                    __FILE__, __LINE__);
            dt = std::max(dt / dscale, mindt);
            oceanTheta->restore();
            continue;
        }
        
        oceanTheta->postProcess();
        
        step++;
        time += dt;
        
        std::stringstream outFile;
        outFile << "ocean_time_" << std::setprecision(8)
                << time * years << ".h5";
        
        if ((output > 0) && ((step % output) == 0))
            oceanTheta->saveStateToFile(outFile.str());

        writeData(oceanTheta, false, step, dt, time*years, k);

        // Timestep adjustments
        if (k < minK)
            dt = std::min(dt * iscale, maxdt);
        else if (k > maxK)
            dt = std::max(dt / dscale, mindt);

        test_step = ( nsteps < 0 ) ? true : step < nsteps;

        // create tangent
        xdot->Update(1.0/dt, *x, -1.0/dt, *xold, 0.0);
        //    INFO(" ||xdot|| = " << Utils::norm(xdot));
    }
    
    // print the profile
    if (Comm->MyPID() == 0)
        printProfile(profile);
    
	//------------------------------------------------------------------
	TIMER_STOP("Total time...");
}

//-----------------------------------------------------------------------------
void writeData(Teuchos::RCP<OceanTheta> model, bool describe,
               int step, double dt, double time, int niters)
{
    // Write continuation data
    std::ostringstream cdatastring;

    if (describe) // write description of entries
    {
        cdatastring << std::setw(_FIELDWIDTH_) 
                    << "time"
                    << std::setw(_FIELDWIDTH_)
                    << "step" 
                    << std::setw(_FIELDWIDTH_) 
                    << "dt" 
                    << std::setw(_FIELDWIDTH_)
                    << "|x|"
                    << std::setw(_FIELDWIDTH_) 
                    << "NR"
                    << model->writeData(describe);
        
        WRITETDATA(cdatastring.str());
    }
    else
    {
        cdatastring.str("");
        cdatastring.clear();
        cdatastring.precision(_PRECISION_);

        cdatastring << std::scientific
                    << std::setw(_FIELDWIDTH_) <<  time
                    << std::setw(_FIELDWIDTH_) <<  step
                    << std::setw(_FIELDWIDTH_) <<  dt  
                    << std::setw(_FIELDWIDTH_)
                    <<  Utils::norm(model->getState('V'))
                    << std::setw(_FIELDWIDTH_) <<  niters
                    << model->writeData();    

        WRITETDATA(cdatastring.str());
    }
}

