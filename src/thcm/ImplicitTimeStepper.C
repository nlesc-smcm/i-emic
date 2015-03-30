/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

#include "Teuchos_RCP.hpp"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_LinearSystem.H"
#include "NOX.H"

#include "Ifpack_Preconditioner.h"

#include "TRIOS_SolverFactory.H"

#include "THCM.H"
#include "ImplicitTimeStepper.H"
#include "ThetaStepperEvaluator.H"
#include "OceanModel.H"

//! this class needs a serious makeover or possibly should be replaced by Rythmos
    
///////////////////////////////////////////////
// Construct from OceanModel and NOX solver //
/////////////////////////////////////////////
ImplicitTimeStepper::ImplicitTimeStepper(Teuchos::RCP<OceanModel> model,
										 Teuchos::ParameterList& stepperParams,
										 Teuchos::ParameterList& nlParams)
	:
	d_model(model),
	d_paramList(stepperParams),
	d_nlParams(nlParams)
{
	DEBUG("ImplicitTimeStepper constructor done");
	
	//== info is a filestream implemented in Filestreams.C
	d_out       = d_paramList.get("Output Stream", info);   
	d_predictor = d_paramList.get("Predictor", "Constant");
	d_scheme    = d_paramList.get("Scheme", "Theta");
	
	if (d_scheme == "Theta")
    {
		d_theta = d_paramList.get("Theta", 1.0); // backward Euler by default
    }
	else if (d_scheme != "Adaptive Theta")
    {
		//== Using __FILE__ and __LINE__ preprocessor macros. 
		Error("only 'Theta' and 'Adaptive Theta' Schemes are implemented!", __FILE__, __LINE__);
    }
	
	// initial step size
	double dt     = d_paramList.get("Step Size" , 1.0);
	double t      = d_paramList.get("Start Time", 0.0);
	d_haveNormEst = true; // Currently we have only theta-schemes that
	                      // support local truncation error (LTE) estimates.

    // create model evaluator for theta time-stepping 
	d_thetaModel = Teuchos::rcp(new ThetaStepperEvaluator(d_model, d_model->getJacobian(), t, d_theta, dt));
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = d_thetaModel;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = d_thetaModel;

    if (iReq == Teuchos::null || iJac == Teuchos::null)
	{
		Error("null pointer detected!", __FILE__, __LINE__);
	}

    // register pre/post operator (TODO: do we need it here?) //== ??
    Teuchos::RCP<NOX::Abstract::PrePostOperator> iPrePost = d_model;

    //Create the "Solver Options" sublist
    Teuchos::ParameterList& optParams = d_nlParams.sublist("Solver Options");

    // we need this so that the vmix_fix flag is handled correctly in THCM //== ??
    optParams.set("User Defined Pre/Post Operator", iPrePost);

	//== document this?
    double TolNewton = d_nlParams.get("Convergence Tolerance", 1.0e-9);
    Teuchos::ParameterList& searchParams = d_nlParams.sublist("Line Search");
    int MaxIt = searchParams.get("Max Iters", 10); 

    // Set up the Solver Convergence tests. //==??
    Teuchos::RCP<NOX::StatusTest::NormF> wrms =
		Teuchos::rcp(new NOX::StatusTest::NormF(TolNewton));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters
		= Teuchos::rcp(new NOX::StatusTest::MaxIters(MaxIt));
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
		Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(wrms);
    combo->addStatusTest(maxiters);

    d_convTests = combo;
#ifdef DEBUGGING //== NEWDEBUG 
	d_convTests->print(*debug);
#endif

    Teuchos::ParameterList& nlPrintParams = d_nlParams.sublist("Printing");
    nlPrintParams.set("MyPID", d_model->get_x_map()->Comm().MyPID());
    nlPrintParams.set("Output Stream", d_out);
    nlPrintParams.set("Error Stream" , d_out);
    nlPrintParams.set("Output Process", 0);
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

	//== Secant prediction
	d_uprime = Teuchos::rcp(new Epetra_Vector(*(d_model->get_x_map())));
	d_ustart = Teuchos::rcp(new Epetra_Vector(*(d_model->get_x_map())));

	//== t-derivative
	d_f_nm1 = Teuchos::rcp(new Epetra_Vector(*(d_model->get_x_map())));
	d_f_n   = Teuchos::rcp(new Epetra_Vector(*(d_model->get_x_map())));
	d_f_np1 = Teuchos::rcp(new Epetra_Vector(*(d_model->get_x_map())));

	DEBUG("Construct Linear System...");
	DEBUG("Construct Linear System...");
	Teuchos::ParameterList& lsParams = d_nlParams.
        sublist("Direction").
        sublist("Newton").
        sublist("Linear Solver");
           
	// note that we want the jacobian/rhs to be evaluated by the ThetaStepperEvaluator,
	// not the original OceanModel
	d_linSys  = THCM::Instance().createLinearSystem(d_model, lsParams, nlPrintParams, d_out,
													d_thetaModel, d_thetaModel);
	NOX::Epetra::Vector noxSoln(*d_ustart);

	DEBUG("ImplicitTimeStepper: create NOX Group");
    // Create the Group
    d_curGroup = Teuchos::rcp(new NOX::Epetra::Group(nlPrintParams, iReq, noxSoln, d_linSys));

	NOX::Solver::Factory nox_factory;
	DEBUG("ImplicitTimeStepper: create NOX Solver");
	d_nlSolver = nox_factory.buildSolver(d_curGroup, d_convTests, Teuchos::rcp(&d_nlParams, false));
	DEBUG("solver constructed successfully");

	DEBUG(" call reset()...");
	this->reset();
    
	(*d_out) << "Stepper Parameters: \n";
	(*d_out) << d_paramList;
	d_restarted = true;
	DEBUG("ImplicitTimeStepper constructor done");
}

////////////////////////////////////////////
// Destructor //////////////////////////////
////////////////////////////////////////////
ImplicitTimeStepper::~ImplicitTimeStepper()
{
}

//////////////////////////////////////////////////////////////
// reset stepper after failed step or at start of new run   //
//////////////////////////////////////////////////////////////
void ImplicitTimeStepper::reset()
{
    DEBUG("ImplicitTimeStepper: reset");
    d_restarted = true;
    d_lteEst    = -1.0;
    d_effortEst =  1.0;
    d_dt_nm1    =  0.0;
    CHECK_ZERO(d_uprime->PutScalar(0.0));
    CHECK_ZERO(d_f_nm1->PutScalar(0.0));
    CHECK_ZERO(d_f_n->PutScalar(0.0));
    DEBUG("ImplicitTimeStepper: reset done");
}

//////////////////////////////////////////////////////////////

bool ImplicitTimeStepper::hasErrorEstimate(void) const
{
    return (d_haveNormEst && (d_lteEst!=-1.0));
}
    
///////////////////////////////////////////////////////////////////////////
// advance one step, Bu_np1 = Bu_n + dt*(theta*F(u_n+1)+(1-theta)F(u_n)) //
///////////////////////////////////////////////////////////////////////////
  
// after a step fails, the user has to call reset(). In that case we need to recompute
// f(u_n) before solving for u_{n+1}, and we are not able to provide an error stimate 
// for u_{n+1}.
bool ImplicitTimeStepper::Step(const Epetra_Vector&   u_n, double t_n, 
							         Epetra_Vector& u_np1, double dt)
{
	if (d_restarted)  // first step or restart after failed step:
    {                 // need to evaluate f_n = F(u_n)
		d_thetaModel->reset(t_n, u_n);
		*d_f_n = d_thetaModel->get_f();
		// destroy preconditioner, regardless of what our strategy may be. This may mean
		// that the new precond is reused for a shorter time.
		d_linSys->destroyPreconditioner();
    }
  
	bool success = true; //indicates wether Newton converged or not
  
	double theta_old = d_theta;
  
	if (d_scheme == "Adaptive Theta")
    {
		double kappa; // kappa is used to select theta somewhere between 0.5 and 1.0, 
		// depending on the time step-size dt
		// TODO: check the correct scaling of dt, it should be 0<ds<=1 here!
		kappa   = std::min(1.0, 0.5 / dt);
		d_theta = 0.5 + kappa*dt;
		(*d_out) << d_scheme << ": " << "Stepper: using theta = " << d_theta << std::endl;
    }
	
	if (d_theta != theta_old)
	{
		d_thetaModel->set_theta(d_theta);
	}
	d_thetaModel->set_dt(dt);
    
	// predictor: we use the secant method or no predictor at all (constant)
	if (d_predictor == "Constant")
    {
		(*d_out) << "Constant Predictor..." << std::endl;
		CHECK_ZERO(d_ustart->Update(1.0, u_n, 0.0));
    }
	else if (d_predictor == "Secant")
    {
		(*d_out) << "Secant Predictor..." << std::endl;
		//== this is the inherited Epetra_MultiVector::Update()
		CHECK_ZERO(d_ustart->Update(1.0, u_n, dt, *d_uprime, 0.0));
    }
	else
    {
		// using preprocessor macros __FILE__ and __LINE__.
		Error("Invalid Predictor for implicit time-integration ", __FILE__, __LINE__);
    }

	// d_nlSolver->reset(d_curGroup, d_convTests, d_nlParams);
	DEBUG("reset the solver...");
	NOX::Epetra::Vector noxvec(*d_ustart);
	d_nlSolver->reset(noxvec);

	// Compute next point on continuation curve
	DEBUG("solve using nox...");
	(*d_out) << "Newton Corrector..." << std::endl;
	
	NOX::StatusTest::StatusType solverStatus = d_nlSolver->solve();  

	// Check solver status
	if (solverStatus == NOX::StatusTest::Failed) 
    {
		DEBUG("WARNING: Nonlinear Solve failed!");
		success = false;
		d_lteEst=-1;
		// skip the rest of the step. The caller has to call reset()
		// and try with a different step-size
		return success;
    }

	// Copy solution out of solver //==?? dynamic casts... why? is there a cleaner way?
	(*d_curGroup) = dynamic_cast<const NOX::Epetra::Group&>(d_nlSolver->getSolutionGroup());
	
	const NOX::Epetra::Vector& noxSoln = dynamic_cast<const NOX::Epetra::Vector&>(d_curGroup->getX());
	u_np1 = noxSoln.getEpetraVector();
			
	int numIters                         = d_nlSolver->getNumIterations();
	Teuchos::ParameterList& searchParams = d_nlParams.sublist("Line Search");
	int maxIters                         = searchParams.get("Max Iters", 10);
  
	// the user can then ask how 'difficult' this step was:
	d_effortEst = ((double) numIters) / ((double) maxIters);
  
	// we can normalize the pressure here, but this
	// is not the same as in the THCM implementation
	
	// TODO: do we need this???
	// THCM::Instance().normalizePressure(u_np1);

	// compute secant to approximate u' for the next predictor
	if (d_predictor == "Secant")
    {
		//== Epetra_MultiVector::Update() //==?? CHECK_ZERO()
		d_uprime->Update(1.0/dt, u_np1, -1.0/dt, u_n, 0.0);
    }
	
	// prepare evaluator for next step
	//== I DON'T SEEM TO UNDERSTAND THIS RESET
	//d_thetaModel->reset(t_n + dt, u_np1);
	
	// get f(u_new) from the theta evaluator
	*d_f_np1 = d_thetaModel->get_f();
  
	double nrm;
	d_f_np1->Norm2(&nrm);
	(*d_out) << "t = " << t_n + dt << ": NORM OF F = " << nrm << std::endl;
  
	if (!d_restarted) // after a reset() we can't do this estimate
    {
		d_lteEst = estimateError(d_theta,*d_f_nm1,*d_f_n,*d_f_np1, d_dt_nm1, dt);
		//(*d_out) << "estimated truncation error: " << d_lteEst << std::endl;
    }

	d_dt_nm1    = dt;
	d_restarted = false;

	// avoid deleting the Teuchos::rcp d_f_nm1:
	Teuchos::RCP<Epetra_Vector> tmp = d_f_nm1;
  
	d_f_nm1  = d_f_n;
	d_f_n    = d_f_np1;
	d_f_np1  = tmp;
	d_f_np1->PutScalar(0.0);

	d_curPrecAge++;
	return success;	
}
// Step

// this error estimate was taken from the THCM code (file time.f version 6.0)
// rh0/1/2 are the value of F at u_{n-1}, u_n and u_{n+1} respectively.
double ImplicitTimeStepper::estimateError(double theta,
										  const Epetra_Vector& rh0,
										  const Epetra_Vector& rh1, 
										  const Epetra_Vector& rh2,
										  double dt_old, double dt) const
{
	Epetra_Vector D(rh0.Map());

	DEBUG("estimating local truncation error...");
	DEBVAR(theta);
	DEBVAR(dt_old);
	DEBVAR(dt);

	Epetra_Vector D1 = rh1;

	double alpha = dt_old/dt;
	double denom = dt_old + alpha*alpha*dt;

	// approximate t-derivative of u at time level n (2nd order accurate)
	//      f_n-1 - (1-a^2) f_n - a^2 f_n+1
	// f' = ---------------------------     + O(dt^2),  a = (dt_n-1) / dt_n
	//             dt_n-1 + a^2 dt_n

	CHECK_ZERO(D1.Update(1.0 / denom, rh0, -(alpha * alpha) / denom, rh2,
						 (alpha * alpha - 1.0) / denom));

	// approximate tt-derivative of u at time level n (2nd order accurate)
	//
	//         f_n-1 - (1+a) f_n + a f_n+1
	// f'' = -------------------------------
	//       (1/2) ((dt_n-1)^2 + a (dt_n)^2)
	//
	
	Epetra_Vector D2 = rh1;
	denom = dt_old*dt_old/2 + alpha*dt*dt/2;
	CHECK_ZERO(D2.Update(1.0/denom,rh0,alpha/denom,rh2,-(1.0+alpha)/denom));

	// estimate LTE
	double c1 = dt*dt*(0.5-theta);
	double c2 = dt*dt*dt*(1.0/6.0 - theta/2.0);
	CHECK_ZERO(D.Update(c1,D1,c2,D2,0.0));

	//DEBVAR(D);

	double error;
	CHECK_ZERO(D.Norm1(&error));
	DEBVAR(error);
	return error;
}

