//=======================================================================
// Main continuation of the ocean model
//=======================================================================

#include "RunDefinitions.H"
#include "OceanTheta.H"
#include "Theta.H"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

using JDQZsolver = JDQZ<JDQZInterface<Teuchos::RCP<Ocean>,
                                      ComplexVector<Epetra_Vector> > >;

//------------------------------------------------------------------
void runOceanModel(RCP<Epetra_Comm> Comm);
void writeData(Teuchos::RCP<Theta<Ocean> > model, bool describe,
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

    double comb =
        oceanParams->sublist("THCM").
        sublist("Starting Parameters").get("Combined Forcing", 99.0);

    if (std::abs(comb) < 1e-7)
    {
        WARNING("Nothing will happen without any forcing: par(comb) = "
                << comb, __FILE__, __LINE__);
    }

	// Create parallelized OceanTheta object
    Teuchos::RCP<Theta<Ocean> > oceanTheta =
        Teuchos::rcp(new Theta<Ocean>(Comm, oceanParams));

    // Create parameter object for time stepping
	RCP<Teuchos::ParameterList> timeParams = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("timestepper_params.xml", timeParams.ptr());
    timeParams->setName("time stepper parameters");

    // get parameters from xml
    double dt     = timeParams->get("initial time step size", 1.0e-03);
    double mindt  = timeParams->get("minimum step size", 1.0e-8);
    double maxdt  = timeParams->get("maximum step size", 1.0);
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
    double time   = 0.0;     // keep track of (dimensional)time
    double normF  = 0.0;     // keep track of residual norm
    double normdx = 0.0;     // keep track of update infnorm
    int    step   = 0;       // keep track of time steps

    Teuchos::RCP<Epetra_Vector> F;
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
    while ( ( time < tend ) && ( test_step ) )
    {
        INFO("----------------------------------------------------------");
        INFO("Timestepping:    t = " << time << " y");

        oceanTheta->preProcess();
        oceanTheta->store();  // save current state to oldstate

 
        // Timestep adjustments
        // while (normF > 10)
        // {
        //     dt = std::max(dt / dscale, mindt);
        //     oceanTheta->setTimestep(dt);
        //     oceanTheta->computeRHS();
        //     normF = Utils::norm(F);
        //     INFO("             ||F|| = " << normF);
        // }

        // Secant prediction (not sure if needed)
        // x->Update(dt, *xdot, 1.0);
        *xold = *x;

        // Newton solve
        int k = 0;
        for (; k != Niters; ++k)
        {
            oceanTheta->setTimestep(dt);

            // compute time discretization -F:
            // -F(x) =  -1 * (B d/dt x / theta + F(x) + (theta-1) / theta * F(x_old))

            oceanTheta->computeRHS();

            // create jacobian of time discretization
            oceanTheta->computeJacobian();

            // solve J(x) dx = -G(x)
            F = oceanTheta->getRHS('V');
            normF = Utils::norm(F);                   
            F->Scale(-1.0);
            oceanTheta->solve(F);

            // correct for pressure modes (superfluous due to prec)
            // oceanTheta->pressureProjection(dx);

            // update state
            x->Update(1.0, *dx, 1.0);

            // compute time discretization -G(x)
            oceanTheta->computeRHS();

            normF  = Utils::norm(F);
            normdx = Utils::normInf(dx);

            INFO("               Newton iter           = " << k);
            INFO("                           ||F||2    = " << normF);
            INFO("                           ||dx||inf = " << normdx);

            if ( (normdx < Ntol ) && (normF < Ntol) )
                break;
            // else if ( normdx > 1e2)
            // {
            //     k = Niters;
            //     break;
            // }
        }

        if (k == Niters)
        {
            WARNING("Newton did not converge! ||F|| = "
                    << normF << "\nRestoring model",
                    __FILE__, __LINE__);
            INFO("    adjusting time step.. old dt = " << dt);
            dt = std::max(dt / dscale, mindt);
            INFO("    adjusting time step.. new dt = " << dt);

            if (dt == mindt)
            {
                INFO("min timestep reached, exiting...");
                return;
            }
            oceanTheta->restore();
            continue;
        }

        step++;
        time += dt * years;

        INFO("           step = " << step);
        INFO("           time = " << time);
        INFO("           Newton solve, ||F||2    = " << normF);
        INFO("           Newton solve  ||dx||inf = " << normdx);

        oceanTheta->postProcess();

        std::stringstream outFile;
        outFile << "ocean_time_" << std::setprecision(8)
                << time << ".h5";

        if ((output > 0) && ((step % output) == 0))
            oceanTheta->saveStateToFile(outFile.str());

        writeData(oceanTheta, false, step, dt, time, k);

        // Timestep adjustments
        if (k < minK)
            dt = std::min(dt * iscale, maxdt);
        else if (k > maxK)
            dt = std::max(dt / dscale, mindt);

        test_step = ( nsteps < 0 ) ? true : step < nsteps;

        // create tangent
        xdot->Update(1.0/dt, *x, -1.0/dt, *xold, 0.0);
        // INFO(" ||xdot|| = " << Utils::norm(xdot));
    }

    // print the profile
    if (Comm->MyPID() == 0)
        printProfile(profile);

	//------------------------------------------------------------------
	TIMER_STOP("Total time...");
}

//-----------------------------------------------------------------------------
void writeData(Teuchos::RCP<Theta<Ocean> > model, bool describe,
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

