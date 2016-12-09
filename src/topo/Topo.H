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
	model_               (model),
	pars_                (pars),
	nMasks_              (pars->get("Number of mask files", 0)),
	delta_               (pars->get("Delta", 0.0)),
	startMask_           (pars->get("Starting mask", 0)),
	stopTol_             (pars->get("Stopping tolerance homotopy", 5.0)),
	storeEverything_     (pars->get("Store everything", false)),
	solverInitialized_   (false),
	usePredictor_        (pars->get("Use predictor", true)),
	useCorrector_        (pars->get("Use corrector", true)),
	maxNRit_             (pars->get("Max Newton steps", 20)),
	tolNR_               (pars->get("Tolerance in Newton iteration", 1e-6)),
	increaseBT_          (pars->get("Test norm increase factor", 2.0)),
	disablePostprocess_  (pars->get("Disable postprocessing", false)),
	recompPreconditioner_(true)
{
	INFO("Topo constructor...");

	double nrm1, nrm2;

	// Setup indices k,a and b
	k_ = startMask_;
	a_ = std::vector<int>(nMasks_);
	b_ = std::vector<int>(nMasks_);	
	a_[0] = 0;
	b_[0] = 1;
	for (int k = 0; k != nMasks_-1; ++k)
	{
		a_[k+1] = a_[k] + 1;
		b_[k+1] = b_[k] + 1;
	}

	// Set the parameters and coefficients
	modelPar_ = model_->getPar();
	setPar(delta_);	

	model_->computeRHS();
	nrm1 = model_->getRHS('V')->norm();

	loadMaskFileNames(); // load filenames
	loadMasks();         // load actual data

	// Set the starting mask A
	model_->setLandMask(landMasks_[a_[k_]], true);
	
	// We initialize our vectors with model vectors.
	stateView_ = model_->getState('V');
	solView_   = model_->getSolution('V');
	stateDot_  = model_->getState('C');
	stateDot_->zero();
	rhs_       = model_->getRHS('C');
	rhsView_   = model_->getRHS('V');
	dFdPar_    = model_->getRHS('C');
	
	model_->computeRHS();
	nrm2 = model_->getRHS('V')->norm();

	// Check whether ocean has been setup with the same a_[k_]-th mask
	if (std::abs(nrm1 - nrm2) > 1e-10)
		ERROR("Norms should not differ much! Choose correct initial landmask!"
			  , __FILE__, __LINE__);

	// A boolean array to indicate whether we have initialized the
	//  associated preconditioner
	initPrecs_ = std::vector<bool>(nMasks_, false);

	// More initializations that we might want to repeat
	initialize();
	
	INFO("Topo constructor... done");
}

//==================================================================
// destructor
template<typename Model, typename ParameterList>
Topo<Model, ParameterList>::
~Topo()
{
	INFO("Topo destructor");
}

//==================================================================
// initialize
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::initialize()
{
	INFO("Topo: initialize...");
	matA_  = rcp(new Matrix(*model_->getJacobian()));
	matB_  = rcp(new Matrix(*model_->getJacobian()));

	// Set landmask (B) in the model
	model_->setLandMask(landMasks_[b_[k_]]);

	// Obtain M, binary diagonal matrix to select transient parts,
	// implemented as a vector
	vecM_ = model_->getM('C');
	
	// Zero out land points
	model_->applyLandMask(vecM_->getOceanVector(), landMasks_[b_[k_]], 0.0);
	
	INFO("Topo: initialize... done");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::predictor()
{
	if (usePredictor_)
	{
		INFO("Topo: predictor...");
		if (!model_->isCurrentMask(landMasks_[b_[k_]].label))
			model_->setLandMask(landMasks_[b_[k_]]);
		
		model_->computeRHS();		
		INFO("  nrm before applyLandmask: " << rhsView_->norm());

		INFO("  k = " << k_ << " a[k] = " << a_[k_] << " b[k] = " << b_[k_]);
		INFO("  using differences between " <<
			 landMasks_[b_[k_]].label << " and " << landMasks_[a_[k_]].label);
			
		model_->applyLandMask(stateView_->getOceanVector(),
							  landMasks_[a_[k_]],
							  landMasks_[b_[k_]]);
	
		model_->computeRHS();

		INFO("  nrm after  applyLandmask: " << rhsView_->norm());
		INFO("Topo: predictor... done");
	}	
	stateCopy_ = model_->getState('C');
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
void Topo<Model, ParameterList>::setMaskIndex(int index)
{
	if ((index < 0) || (index > nMasks_))
		ERROR("Wrong mask index", __FILE__, __LINE__);

	INFO("Topo: setting mask idx: " << index);
	k_ = index;
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::setPar(double par)
{
	if (par < 0 || par > nMasks_)
	{
		std::stringstream err;
		err << "Topo: invalid parameter, par = " << par; 
		ERROR(err.str(), __FILE__, __LINE__);
	}

	// set parameter
	delta_ = par;

	// calculate coefficients
	// facA_ = 1 - delta_;
	// facB_ = delta_;

	// set coefficients
	
	facA_  = pow(cos(M_PI * delta_ / 2), 2);
	facB_  = pow(sin(M_PI * delta_ / 2), 2);

	// correct for small values	
	facA_ = (std::abs(facA_) < 1e-16) ? 0 : facA_;
	facB_ = (std::abs(facB_) < 1e-16) ? 0 : facB_;
	
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
	// Get a copy of RHS (A)
	VectorPtr x = model_->getState('C');

	// Overwrite state with initial
	x->update(-1.0, *stateCopy_, 1.0);
	INFO("            norm x-xA = " << x->norm());

	// dangerous
	CHECK_ZERO(x->getOceanVector()->Multiply(1.0, *(x->getOceanVector()),
											 *(vecM_->getOceanVector()), 0.0));
	
	// Set landmask (B) in the model
	if (!model_->isCurrentMask(landMasks_[b_[k_]].label))
		model_->setLandMask(landMasks_[b_[k_]]);

	// Evaluate and compute RHS (B)
	model_->computeRHS();
	
	// Get the view	
	*rhs_ = *model_->getRHS('V');

	normfB_ = rhs_->norm();
	INFO("              norm fB = " << normfB_);

	// Scale the RHS, skipping land, w and p rows.

	Vector scaleVec(*vecM_);	
		
	int NumMyElements  = scaleVec.getOceanVector()->Map().NumMyElements();
	
	for (int i = 0; i != NumMyElements; ++i)
		if  (std::abs((*scaleVec.getOceanVector())[i]) < 1e-12)
			(*scaleVec.getOceanVector())[i] = 1.0;
		else
			(*scaleVec.getOceanVector())[i] = facB_;

	CHECK_ZERO(rhs_->getOceanVector()->Multiply(1.0, *(rhs_->getOceanVector()),
												*(scaleVec.getOceanVector()), 0.0));

	rhs_->update(facA_, *x, 1.0);

	INFO("                 facA = " << facA_);
	INFO("                 facB = " << facB_);
	INFO("                gamma = " << gamma_);
	INFO("facA*(x-xA) + facB*fB = " << rhs_->norm());
}

//==================================================================
template<typename Model, typename ParameterList>
bool Topo<Model, ParameterList>::monitor()
{
	return (normfB_ < stopTol_) || (delta_ > 1.0);
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::computeDFDPar()
{
	INFO(" TOPO: Calculating dFdpar...");
 	// Set landmask (A) in the model
	model_->setLandMask(landMasks_[a_[k_]]);	
	
	// Evaluate and compute RHS (A)
	model_->computeRHS();

	// Get a copy of RHS (A)
	VectorPtr rhsA = model_->getRHS('C');

	// Set landmask (B) in the model
	model_->setLandMask(landMasks_[b_[k_]]);

	// Evaluate and compute RHS (B)
	model_->computeRHS();

	// Get it
	*dFdPar_ = *model_->getRHS('C');

	// Calculate pi*cos(theta)*sin(theta)(F(B) - F(A))
	double factor = M_PI * sin(M_PI * delta_ / 2)
		* cos(M_PI * delta_ / 2);

	dFdPar_->update(-factor, *rhsA, factor);

	//done
	INFO(" TOPO: Calculating dFdpar... done " << dFdPar_->norm());
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::computeJacobian()
{
	INFO("Topo: compute Jacobian...");

	if (!model_->isCurrentMask(landMasks_[b_[k_]].label))
		model_->setLandMask(landMasks_[b_[k_]]);
	
	model_->computeJacobian();
	matB_ = model_->getJacobian();	

	Vector leftScale(*vecM_);

	// to maintain correct sparsity pattern
	// double factor = (std::abs(facB_) < 1e-16) ? 1e-16 : facB_;
	
	// leftScale.scale(facB_);

	// Set zeros to 1, let scaling have no effect at land and w,p points.
	// Apply scaling to other points
	int NumMyElements = leftScale.getOceanVector()->Map().NumMyElements();
	for (int i = 0; i != NumMyElements; ++i)
		if (std::abs((*leftScale.getOceanVector())[i]) < 1e-12)
			(*leftScale.getOceanVector())[i] = 1.0;
		else
			(*leftScale.getOceanVector())[i] = facB_;
	
	//Row scale matrix
	CHECK_ZERO( matB_->LeftScale(*(leftScale.getOceanVector()) ));
	
	// Obtain diagonal of matrix
	VectorPtr tmp = model_->getSolution('C');
	tmp->putScalar(0.0);
	matB_->ExtractDiagonalCopy(*(tmp->getOceanVector()));
	tmp->update(facA_, *vecM_, 1.0);
	matB_->ReplaceDiagonalValues(*(tmp->getOceanVector()));

	matB_->SetLabel("TopoCont Jacobian");
	
	combMat_.A      = matA_;
	combMat_.B      = matB_;
	combMat_.PA     = PreconPtr();
	combMat_.PB     = PreconPtr();
	combMat_.facA   = facA_;
	combMat_.facB   = facB_;		
	combMat_.isPrec = false;
	combMat_.precType = 0;

	// We don't have to update the preconditioners
	// but we do need to let them know the correct factors!
	combPrec_.facA   = facA_;
	combPrec_.facB   = facB_;
	INFO("Topo: compute Jacobian... done");
}

//==================================================================
// DEPRECATED
template<typename Model, typename ParameterList>
double Topo<Model, ParameterList>::computeResidual(VectorPtr b)
{
	return 0.0;
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::buildPreconditioner()
{
	if (!recompPreconditioner_)
		return;

	INFO("Topo: build preconditioner...");
	
	recompPreconditioner_ = false;
	
	combPrec_.facA   = facA_;
	combPrec_.facB   = facB_;

	// If the landmask is new we need to rebuild the preconditioner
	if ( !initPrecs_[b_[k_]] )
	{
		initPrecs_[b_[k_]] = true;
		model_->buildPreconditioner(true);
		precB_ = model_->getPreconPtr();
	}

	// Compute preconditioner
	precB_->Compute();

	combPrec_.A        = matA_;
	combPrec_.PA       = PreconPtr(); // DEPRECATED
	combPrec_.B        = matB_;	
	combPrec_.PB       = precB_;
	combPrec_.isPrec   = true;
	combPrec_.precType = pars_->get("Preconditioning type", 'H');

	INFO("Topo: build preconditioner... done");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::initializeSolver()
{
	if (solverInitialized_)
		return;
	
	INFO("Topo: initialize solver...");
	solverParams_ = rcp(new Teuchos::ParameterList);
	updateParametersFromXmlFile("solver_params.xml", solverParams_.ptr());
	
	// Belos LinearProblem setup
	Teuchos::RCP<Combined_Operator<Model> > combMatPtr =
		Teuchos::rcp(&combMat_, false);
	
	problem_ =
		rcp(new Belos::LinearProblem
			<double, Epetra_MultiVector, Combined_Operator<Model> >
			(combMatPtr, solView_->getOceanVector(), rhs_->getOceanVector()));
			
	buildPreconditioner();

	Teuchos::RCP<Combined_Operator<Model> > combPrecPtr =
		Teuchos::rcp(&combPrec_, false);
	
	problem_->setRightPrec(combPrecPtr);

	// A few FGMRES parameters are made available in solver_params.xml:
	int gmresIters  = solverParams_->get("Topo FGMRES iterations", 400);
	double gmresTol = solverParams_->get("Topo FGMRES tolerance", 1e-2);
	int maxrestarts = solverParams_->get("Topo FGMRES restarts", 2);
	int output      = solverParams_->get("Topo FGMRES output", 20);

	int NumGlobalElements = stateView_->getOceanVector()->GlobalLength();
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
	belosParamList_->set("Output Frequency", output);
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
	INFO("Topo: initialize solver... done");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::solve(VectorPtr b)
{
	
	initializeSolver();    // Initialize solver
	buildPreconditioner(); // Build preconditioner			

	TIMER_START("  TOPO:  solve...");
	INFO("  TOPO:  solve...");

	int    innerIters;
	double innerTol, outerTol;
	bool   set;
	
	Vector resB = *getSolution('C');	

	solView_->putScalar(0.0);

	set = problem_->setProblem(solView_->getOceanVector(), b->getOceanVector());

	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
							   "*** Belos::LinearProblem failed to setup");
	try
	{
		belosSolver_->solve(); 	// Solve
	}
	catch (std::exception const &e)
	{
		INFO("Ocean: exception caught: " << e.what());
	}

	innerIters = belosSolver_->getNumIters();
	innerTol   = belosSolver_->achievedTol();

	INFO("  TOPO: FGMRES, inner = " <<
		 innerIters << ", ||r|| = " << innerTol);
	
	INFO("  TOPO: FGMRES  ||x|| = " << solView_->norm());
	
	combMat_.B->Apply(*solView_->getOceanVector(),
					  *resB.getOceanVector());

	resB.update(1.0, *b, -facB_);
	resB.update(-facA_, *solView_, 1.0);
		
	outerTol = resB.norm();
	INFO("  TOPO: FGMRES, actual ||r|| = " << outerTol);
		
	TRACK_ITERATIONS("  TOPO: FGMRES iterations...", innerIters);		
	
	INFO("  TOPO:  solve... done");
	TIMER_STOP("  TOPO:  solve...");
}

//==================================================================
template<typename Model, typename ParameterList>
int Topo<Model, ParameterList>::corrector()
{
	INFO("Topo: corrector...");

	if (!useCorrector_)
	{
		INFO("Topo: corrector disabled! ");
		return 0;
	}
	
	model_->setLandMask(landMasks_[b_[k_]]);	
	model_->computeJacobian();
	rhsView_->putScalar(0.0);
	solView_->putScalar(0.0);
	model_->computeRHS();
	
	// force initializations VERY IMPORTANT!
	model_->buildPreconditioner(true);
	model_->initializeSolver();

	double oldNorm = normfB_;
	// double iniNorm = normfB_;

	int	iter = 0;
	int maxBTiter = 10;

	while ((iter < maxNRit_) && (normfB_ > tolNR_))
	{
		iter++;
		rhsView_->scale(-1.0);
		model_->solve(rhsView_);
		model_->printResidual(rhsView_);
		stateView_->update(1.0, *solView_, 1.0);
		model_->computeRHS();
 		normfB_ = rhsView_->norm();

		INFO("  NR " << iter);
		INFO("  NR norm fB = " << normfB_);
		INFO("    testnorm = " << increaseBT_ * oldNorm);
		INFO("        tol  = " << tolNR_);

		// relaxed backtracking
		double r = -1.0; // reduction
		int iterBT = 0;
		while (normfB_ > increaseBT_ * oldNorm && iterBT < maxBTiter)
		{
			iterBT++;
			r /= 2.0;
			INFO("                   r  = " << r);
			stateView_->update(r, *solView_, 1.0);
			model_->computeRHS();
			normfB_ = rhsView_->norm();
			INFO("   back track norm fB = " << normfB_);
		}
		
		oldNorm = normfB_;
		model_->computeJacobian();
	}

	if (iter >= maxNRit_)
	{
		INFO("Topo: corrector failed! ");
		return 1;
	}
			
	INFO("Topo: corrector... done");
	return 0;
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
	solverInitialized_    = false;  
	recompPreconditioner_ = true; 
	INFO("Topo  pre-processing: enabling computation of preconditioner.");
	model_->preProcess();
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::postProcess()
{
	if (!disablePostprocess_)
	{
		// For plotting we use the landmask nearest to delta, so we round.
		int r = static_cast<int>(std::round(delta_));
		model_->setLandMask(landMasks_[k_+r], true);
		model_->postProcess();

		if (storeEverything_)
		{
			std::stringstream fname1, fname2;
			fname1 << "state_topo_" << std::setprecision(12)
				   << k_ + delta_;
			model_->copyFiles(fname1.str());
			fname2 << "mask_" << std::setprecision(4)
				   << k_ + r;
			model_->copyMask(fname2.str());
		}
	}
	else
		INFO("Topo: postprocessing disabled!");
}

#endif
