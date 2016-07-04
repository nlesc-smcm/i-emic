//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef TOPO_H
#define TOPO_H

#include "TopoDecl.H"
#include "GlobalDefinitions.H"

#include <vector>
#include <sstream>
#include <math.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

//=====================================================================
using Teuchos::RCP;
using Teuchos::rcp;

//==================================================================
// constructor
template<typename Model, typename ParameterList>
Topo<Model, ParameterList>::
Topo(Model model, ParameterList pars)
	:
	model_              (model),
	pars_               (pars),
	nMasks_             (pars->get("Number of mask files", 0)),
	delta_              (pars->get("Delta", 0.0)),
	storeEverything_    (pars->get("Store everything", false)),
	solverInitialized_  (false)
{
	INFO("Topo constructor...");

	double nrm1, nrm2;

	// Set the parameter and coefficients
	setPar(delta_);

	model_->computeRHS();
	nrm1 = model_->getRHS('V')->norm();
	INFO("   norm RHS (1) = " << nrm1);

	loadMaskFileNames();
	loadMasks();
	
	model_->setLandMask(landMasks_[0], true);
	
	// We initialize our vectors with the ones in the model.
	stateView_ = model_->getState('V');
	solView_   = model_->getSolution('V');
	rhs_       = model_->getRHS('C');	

	// Setup indices k,a and b
	k_ = static_cast<int>(std::floor(delta_));
	a_ = std::vector<int>(nMasks_);
	b_ = std::vector<int>(nMasks_);		
	a_[0] = 0;
	b_[0] = 1;
	for (int k = 0; k != nMasks_-1; ++k)
	{
		a_[k+1] = a_[k] + 2 * ((k + 1) % 2);
		b_[k+1] = b_[k] + 2 * ( k % 2 );
	}	

	model_->computeRHS();
	nrm2 = model_->getRHS('V')->norm();
	INFO("   norm RHS (2) = " << nrm2);
	
	if (std::abs(nrm1 - nrm2) > 1e-10)
		ERROR("Norms should not differ much!", __FILE__, __LINE__);

	matA_  = rcp(new Matrix(*model_->getJacobian()));
	matB_  = rcp(new Matrix(*model_->getJacobian()));

	initPrecs_ = std::vector<bool>(nMasks_, false);
	
	INFO("Topo constructor... done");
}

//==================================================================
// constructor
template<typename Model, typename ParameterList>
Topo<Model, ParameterList>::
~Topo()
{
	INFO("Topo destructor");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::loadMaskFileNames()
{
	INFO(" Topo: loading mask files...");
	
	landMaskFileNames_ = std::vector<std::string>(nMasks_);
	std::stringstream fileID;
	for (int i = 0; i != nMasks_; ++i)
	{
		fileID << "Mask file " << i;		
		landMaskFileNames_[i] = pars_->get(fileID.str().c_str(), "");
		assert(landMaskFileNames_[i] != "");
		INFO("  " << fileID.str() << ": " << landMaskFileNames_[i]);
		fileID.str("");
		fileID.clear();
	}
	
	INFO(" Topo: loading mask files... done");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::loadMasks()
{
	INFO(" Topo: loading mask arrays...");

	landMasks_ = std::vector<LandMask>(nMasks_);
	
	for (int i = 0; i != nMasks_; ++i)
	{
		INFO("   loading " << landMaskFileNames_[i]);
		landMasks_[i] = model_->getLandMask(landMaskFileNames_[i]);
	}
	
	INFO(" Topo: loading mask arrays... done");
}

//==================================================================
template<typename Model, typename ParameterList>
double Topo<Model, ParameterList>::getPar()
{
	return delta_;
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::setPar(double par)
{
	if (par < 0 || par > nMasks_ - 1)
	{
		std::stringstream err;
		err << "Topo: invalid parameter, par = " << par; 
		ERROR(err.str(), __FILE__, __LINE__);
	}

	// set parameter
	delta_ = par;

	// calculate index and coefficients
	k_     = static_cast<int>(std::floor(delta_));
	facA_  = pow(cos(M_PI * delta_ / 2), 2);
	facB_  = pow(sin(M_PI * delta_ / 2), 2);

	// ALTERNATIVE
	// double dhat = delta_ - k_;
	// facA_ = (k_ % 2 == 0) ? 1 - dhat : dhat;
	// facB_ = (k_ % 2 == 0) ? dhat : 1 - dhat; 	

	// correct for small values
	facA_ = (std::abs(facA_) < 1e-16) ? 0 : facA_;
	facB_ = (std::abs(facB_) < 1e-16) ? 0 : facB_;
	
	INFO("  <><>  delta = " << delta_ <<
		 " k = " << k_ <<
		 " facA = " << facA_ <<
		 " facB = " << facB_  << "  <><>  ");

	if (std::abs(facA_ + facB_ - 1.0) > 1e-7)
	{
		std::stringstream errmsg;
		errmsg << "facA and facB do not add up to 1.0";
		errmsg << " facA + facB = " << facA_ + facB_;
		ERROR(errmsg.str(), __FILE__, __LINE__);
	}
}

//==================================================================
template<typename Model, typename ParameterList>
typename Topo<Model, ParameterList>::VectorPtr
Topo<Model, ParameterList>::getState(char mode)
{
	if (mode == 'V')
	{
		return stateView_;
	}
	else if (mode == 'C')
	{
		VectorPtr stateCopy = model_->getState('C');
		// copy construction
		*stateCopy = *stateView_;
		return stateCopy;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}
}

//==================================================================
template<typename Model, typename ParameterList>
typename Topo<Model, ParameterList>::VectorPtr
Topo<Model, ParameterList>::getRHS(char mode)
{
	if (mode == 'V')
	{
		return rhs_;
	}
	else if (mode == 'C')
	{
		VectorPtr rhsCopy = model_->getRHS('C');
		// copy construction
		*rhsCopy = *rhs_;
		return rhsCopy;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//==================================================================
template<typename Model, typename ParameterList>
typename Topo<Model, ParameterList>::VectorPtr
Topo<Model, ParameterList>::getSolution(char mode)
{
	if (mode == 'V')
	{
		return solView_;
	}
	else if (mode == 'C')
	{
		VectorPtr solCopy = model_->getSolution('C');
		// copy construction
		*solCopy = *solView_;
		return solCopy;
	}
	else
	{
		WARNING("Invalid mode", __FILE__, __LINE__);
		return Teuchos::null;
	}	
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::computeRHS()
{
 	// Set landmask (A) in the model
	model_->setLandMask(landMasks_[a_[k_]]);	
	
	// Evaluate and compute RHS (A)
	model_->computeRHS();

	// Get a copy of RHS (A)
	VectorPtr rhsA = model_->getRHS('C');

	INFO("    ||F(A)|| = " <<  rhsA->norm() <<
		 " facA " << facA_ << " : " << facA_ * rhsA->norm());

	// Set landmask (B) in the model
	model_->setLandMask(landMasks_[b_[k_]]);

	// Evaluate and compute RHS (B)
	model_->computeRHS();
	
	// Let our rhs be a view
	*rhs_ = *model_->getRHS('V');

	INFO("    ||F(B)|| = " <<  rhs_->norm() <<
		 " facB " << facB_ << " : " << facB_ * rhs_->norm());

	// Add facB * RHS(B) to facA * RHS(A)
	rhs_->update(facA_, *rhsA, facB_);

	INFO("    ||F(A) + F(B)|| = " << rhs_->norm());
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::computeJacobian()
{
	// Get Jacobians
	model_->setLandMask(landMasks_[a_[k_]]);

	model_->computeJacobian();
	*matA_ = *model_->getJacobian();

	model_->setLandMask(landMasks_[b_[k_]]);

	model_->computeJacobian();
	*matB_ = *model_->getJacobian();

	combMat_.A      = matA_;
	combMat_.B      = matB_;
	combMat_.PA     = PreconPtr();
	combMat_.PB     = PreconPtr();
	combMat_.facA   = facA_;
	combMat_.facB   = facB_;		
	combMat_.isPrec = false;
}

//==================================================================
template<typename Model, typename ParameterList>
double Topo<Model, ParameterList>::computeResidual(VectorPtr b)
{
	Vector resA = *getSolution('C');
	Vector resB = *getSolution('C');
	resA.zero();
	resB.zero();
	
	combMat_.A->Apply(*solView_->getOceanVector(),
					  *resA.getOceanVector());

	combMat_.B->Apply(*solView_->getOceanVector(),
					  *resB.getOceanVector());
	resA.scale(facA_);
	resB.scale(facB_);
	
	resA.update(1.0, *b, -1.0);
	INFO("   TOPO residual ||b-Ax|| / ||b||: " << resA.norm() / b->norm());
	INFO("                 ||b-Bx|| / ||b||: " << resB.norm() / b->norm());
	resB.update(1.0, resA, -1.0);
	INFO("             ||b-(A+B)x|| / ||b||: " << resB.norm() / b->norm());
	return resB.norm();	
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::buildPreconditioner()
{
	combPrec_.A      = MatrixPtr();
	combPrec_.B      = MatrixPtr();
	combPrec_.facA   = facA_;
	combPrec_.facB   = facB_;		

	//------------------------------------------------------------------
	// Preconditioner A
	//------------------------------------------------------------------
	
	model_->setLandMask(landMasks_[a_[k_]]);
	model_->computeJacobian();
	
	if ( !initPrecs_[a_[k_]] )
	{
		initPrecs_[a_[k_]] = true;
		model_->buildPreconditioner(true);
		precA_ = model_->getPreconPtr();
	}	

	precA_->Compute();
	
	//------------------------------------------------------------------
	// Preconditioner B
	//------------------------------------------------------------------

	model_->setLandMask(landMasks_[b_[k_]]);
	model_->computeJacobian();
	
	if ( !initPrecs_[b_[k_]] )
	{
		initPrecs_[b_[k_]] = true;
		model_->buildPreconditioner(true);
		precB_ = model_->getPreconPtr();
	}

	precB_->Compute();		

	//------------------------------------------------------------------
	
	combPrec_.PA = precA_;
	combPrec_.PB = precB_;
	combPrec_.isPrec = true;
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::initializeSolver()
{
	solverParams_ = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams_.ptr());
	
	// Belos LinearProblem setup
	Teuchos::RCP<Combined_Operator<Model> > combMatPtr =
		Teuchos::rcp(&combMat_, false);
	
	problem_ =
		rcp(new Belos::LinearProblem
			<double, Epetra_MultiVector, Combined_Operator<Model> >
			(combMatPtr,
			 solView_->getOceanVector(),
			 rhs_->getOceanVector()) );

	buildPreconditioner();

	Teuchos::RCP<Combined_Operator<Model> > combPrecPtr =
		Teuchos::rcp(&combPrec_, false);
	
	problem_->setRightPrec(combPrecPtr);

	// A few FGMRES parameters are made available in solver_params.xml:
	int gmresIters  = solverParams_->get("FGMRES iterations", 500);
	double gmresTol = solverParams_->get("FGMRES tolerance", 1e-8);

	int NumGlobalElements = stateView_->getOceanVector()->GlobalLength();
	int maxrestarts       = 0;
	int blocksize         = 1; // number of vectors in rhs
	int maxiters          = NumGlobalElements / blocksize - 1;

	// Create Belos parameterlist
	belosParamList_ = rcp(new Teuchos::ParameterList());
	belosParamList_->set("Block Size", blocksize);
	belosParamList_->set("Flexible Gmres", true);
	belosParamList_->set("Adaptive Block Size", true);
	belosParamList_->set("Num Blocks", gmresIters);
	belosParamList_->set("Maximum Restarts", maxrestarts);
	belosParamList_->set("Orthogonalization","DGKS");
	belosParamList_->set("Output Frequency", 100);
	belosParamList_->set("Verbosity", Belos::TimingDetails +
						 Belos::Errors +
						 Belos::Warnings +
						 Belos::StatusTestDetails );
	belosParamList_->set("Maximum Iterations", maxiters); 
	belosParamList_->set("Convergence Tolerance", gmresTol); 
	belosParamList_->set("Explicit Residual Test", false); 
	belosParamList_->set("Implicit Residual Scaling", "Norm of RHS");

	// Belos block FGMRES setup
	belosSolver_ =
		rcp(new Belos::BlockGmresSolMgr
			<double, Epetra_MultiVector, Combined_Operator<Model> >
			(problem_, belosParamList_));

	solverInitialized_ = true;
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::solve(VectorPtr b)
{
	// Check whether solver is initialized, if not perform the
	// initialization here
	if (!solverInitialized_)
		initializeSolver();	

	// Initialize to zero
	solView_->putScalar(0.0);

	// int dominantMask = (facA_ >= facB_) ? a_[k_] : b_[k_];
	// model_->setLandMask(landMasks_[dominantMask]);
	// model_->computeJacobian();
	// model_->solve(b);
	
	bool set;
	set = problem_->setProblem(solView_->getOceanVector(), b->getOceanVector());
	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
							   "*** Belos::LinearProblem failed to setup");

	TIMER_START("  TOPO:  solve...");
	INFO("  TOPO:  solve...");
	int    iters;
	double tol;

	try
	{
		belosSolver_->solve(); 	// Solve
	}
	catch (std::exception const &e)
	{
		INFO("Ocean: exception caught: " << e.what());
	}

	INFO("  TOPO:  solve... done");
	TIMER_STOP("  TOPO:  solve...");

	iters = belosSolver_->getNumIters();
	tol   = belosSolver_->achievedTol();
	INFO("  TOPO: FGMRES, i = " << iters << ", ||r|| = " << tol);
	
	TRACK_ITERATIONS("  TOPO: FGMRES iterations...", iters);

	computeResidual(b);
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::applyMatrix(Vector const &v, Vector &out)
{}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::applyPrecon(Vector const &v, Vector &out)
{}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::preProcess()
{
	model_->preProcess();
	buildPreconditioner();
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::postProcess()
{
	// For plotting we use the landmask nearest to delta, so we round.
	int r = static_cast<int>(std::round(delta_));
	model_->setLandMask(landMasks_[r], true);

	// If we are near a landmask we reset new land values
	if (std::abs(facA_) < 1e-3)
	{
		model_->applyLandMask(landMasks_[b_[k_]], 0.0);
	}
	else if (std::abs(facB_) < 1e-3)
	{
		model_->applyLandMask(landMasks_[a_[k_]], 0.0);
	}
	
	model_->postProcess();

	if (storeEverything_)
	{
		std::stringstream fname1, fname2;
		fname1 << "state_topo_" << std::setprecision(4)
			  << delta_;
		model_->copyFiles(fname1.str());
		fname2 << "mask_" << std::setprecision(4)
			   <<  r;
		model_->copyMask(fname2.str());
	}
}

#endif
