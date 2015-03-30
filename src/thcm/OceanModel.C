/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#include <iostream>
#include <iomanip>
#include <stdio.h>

/* Trilinos */

// Teuchos
#include "Teuchos_oblackholestream.hpp"

// Epetra
#include "Epetra_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"

#ifdef STORE_MATRICES
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#endif

// LOCA
#include "LOCA_Parameter_Vector.H"

/* my own packages */

// HYMLS
#include "HYMLS_Tools.H"  // added_14_10_2014_erik
#include "HYMLS_MatrixUtils.H"

// TRIOS
#include "TRIOS_Domain.H"
#include "TRIOS_BlockPreconditioner.H"
#include "TRIOS_HYMLSolver.H"

//trilinos_thcm
#include "globdefs.H"
#include "THCM.H"
#include "OceanModel.H"
#include "OceanGrid.H"

#include "OceanOutputXdmf.H"

// thcm

extern "C" {
	_SUBROUTINE_(get_grid_data)(int* nn,int* mm,int* ll,double* xx,double* yy,double* zz);
	_SUBROUTINE_(write_data)(double*, int*, int*);  // file inout.F90
	_SUBROUTINE_(append_data)(double*, int*, int*);  // file inout.F90
	_SUBROUTINE_(writematrhs)(double*);  // file matetc.F90
	_MODULE_SUBROUTINE_(m_global,compute_flux)(double*);
}//extern

OceanModelEvaluator::OceanModelEvaluator(ParameterList& plist):
	paramList(plist),pVector(null)
{
	DEBUG("Create Ocean Model...");
    
	DEBVAR(paramList);
  
#ifdef STORE_MATRICES
	store_step_jac = 0;
	store_step_rhs = 0;
#endif

	std::string probdesc = paramList.get("Problem Description","Unnamed");
	this->SetLabel(("OceanModel ("+probdesc+")").c_str());

// continuation parameter (may be "Time" in transient mode)
// This is typically set by the main program
	cont_param = paramList.get("Parameter Name","Undefined");
	// sign s of exponent if you want to change e in param = 10^{s*e}
	cont_s   = paramList.get("Continuation in Exponent",1.0);
	// if you use "Exponent" as cont. param, the following parameter
	// is given the value 10^{s*e}:
	exp_cont_param = paramList.get("Exp. Cont. Parameter","Undefined");
  
	pVector = THCM::Instance().getParameterVector();
  
	// update parameter names and initial values for ModelEvaluator
	Epetra_SerialComm scomm;
	int npar = _NPAR_+_NPAR_TRILI+1;//+1 because we also have time.
	// time is par(0) and some more are added in THCM.C 
	// (_NPAR_TRILI)
	p_map   = Teuchos::rcp(new Epetra_Map(npar,0,scomm));
	p_init  = Teuchos::rcp(new Epetra_Vector(*p_map));
	p_names = Teuchos::rcp(new Teuchos::Array<std::string>(npar));
	double start_value;
	for (int i=0;i<pVector->length();i++)
    {
		(*p_names)[i] = pVector->getLabel(i);
		(*p_init)[i] = pVector->getValue(i);
    }
	try {
		start_value = pVector->getValue(cont_param);
    } catch(...) {
        Error("specified Parameter not found in ParameterVector!",__FILE__,__LINE__);
	}

	prec_reuse_policy = paramList.get("Preconditioner Reuse Policy", "None");
	// if this option is not set, leave it to LOCA:
	max_prec_age = paramList.get("Max Age Of Prec", -2);
	prec_age=0;

	backup_interval = paramList.get("Backup Interval",-1.0);
	last_backup = start_value-1e-12;// we subtract eps because of the way 

	// backuping is treated in XYZT mode. 
	// otherwise the initial solution would
	// be bckuped in that case.
	output_interval=paramList.get("Output Frequency",-1.0);
	last_output=last_backup-1.1*output_interval;
	ParameterList& xdmfParams = paramList.sublist("Xdmf Output");
	Teuchos::RCP<TRIOS::Domain> domain = THCM::Instance().GetDomain();
	bool print = (output_interval>0);
	xdmfWriter = Teuchos::rcp(new OceanOutputXdmf(domain,xdmfParams,print));
	gridPtr = xdmfWriter->getGrid();

	// pressure correction not implemented correctly
	pres_corr = "never";
    
}

OceanModelEvaluator::~OceanModelEvaluator()
{ 
	DEBUG("Destroy Ocean Model...");
}

Teuchos::RCP<Epetra_Vector> OceanModelEvaluator::ReadConfiguration(std::string filename ,LOCA::ParameterVector& pVec)
{

	Teuchos::RCP<std::istream> in;
	in = Teuchos::rcp(new std::ifstream(filename.c_str()) );
	std::string s1,s2,s3;
	(*in) >> s1;
	DEBVAR(s1);
	if (s1!="LOCA::ParameterVector")
    {
		Error("Error reading start config",__FILE__,__LINE__);
    }

	// read THCM Parameter vector
	int npar;
	(*in) >> s1 >> s2 >> npar >>s3;
	DEBVAR(npar)

		int j;
	std::string key;
	double value;
	for (int i=0;i<npar;i++)
    {
		read_parameter_entry(in,key,value);
		if (pVec.isParameter(key))
		{
			pVec.setValue(key, value);
		} 
		else
		{
			pVec.addParameter(key,value);
		}
    }

	// read current solution
	Teuchos::RCP<Epetra_Vector> dsoln = THCM::Instance().getSolution();
	Teuchos::RCP<Epetra_Map>     dmap = THCM::Instance().GetDomain()->GetSolveMap();

	Teuchos::RCP<Epetra_MultiVector> gsoln = HYMLS::MatrixUtils::Gather(*dsoln,0);

	if (THCM::Instance().GetComm()->MyPID()==0)
    {
		(*in) >> s1;
		if (s1!="Epetra::Vector")
		{
			INFO("Bad Vector label: should be Epetra::Vector, found "<<s1<<std::endl);
			Error("Error reading start config",__FILE__,__LINE__);
		}
		(*in) >> s1 >> s2 >> s3;
		if (s1+s2+s3!="MyPIDGIDValue")
		{
			Error("Error reading start config",__FILE__,__LINE__);
		}
		int pid,gid;
		double val;
		for (int i=0;i<gsoln->GlobalLength();i++)
		{
			(*in) >> pid >> gid >> val;
			(*gsoln)[0][gid]=val;
		}
    }
    
	Teuchos::RCP<Epetra_MultiVector> tmp = HYMLS::MatrixUtils::Scatter(*gsoln,*dmap);
	*dsoln = *((*gsoln)(0));
  
	try {
		last_backup=pVec.getValue(cont_param);
    } catch (...) {
		Error("Missing continuation parameter in starting file!",__FILE__,__LINE__);
    }
	return dsoln;
}
    
void OceanModelEvaluator::read_parameter_entry(RCP<std::istream> in,std::string& key, double& value)
{
	int j;
	std::string tmp;
	(*in) >> j;
	(*in) >> key;
	while (1)
    {
		(*in) >> tmp;
		if (tmp=="=")
		{
			(*in) >> value;
			break;
		}
		else
		{
			key = key + " "+tmp;
		}
    }
	DEBVAR(j);
	DEBVAR(key);
	DEBVAR(value);
}
    

void OceanModelEvaluator::WriteConfiguration(std::string filename , const LOCA::ParameterVector& 
											 pVector, const Epetra_Vector& soln)
{
	Teuchos::RCP<std::ostream> out;
	if (THCM::Instance().GetComm()->MyPID()==0)
    {
		out = Teuchos::rcp(new std::ofstream(filename.c_str()) );
    }
	else
    { // dummy stream
		out = Teuchos::rcp(new Teuchos::oblackholestream());
    }
	(*out) << std::setw(15) << std::setprecision(15);
	out->setf(std::ios::scientific); 
	(*out) << pVector;
	(*out) << *(HYMLS::MatrixUtils::Gather(soln,0));
}

// compute and store streamfunction in 'fort.7'
void OceanModelEvaluator::Monitor(double conParam)
{
	// data in grid-object is assumed to be current solution

	// some constants (these should be consistent with
	// those in usr.F90)
	const double r0dim = 6370.0e3;
	const double udim = 0.1;

	double hdim = THCM::Instance().hDim();
	double transc = r0dim*hdim*udim*1e-6;


	double psimmin = transc*gridPtr->psimMin();
	double psimmax = transc*gridPtr->psimMax();
	double psibmin = transc*gridPtr->psibMin();
	double psibmax = transc*gridPtr->psibMax();
  
	int itp = 0; // bifurcation point? Can't say up to now!
	int icp = THCM::Instance().par2int(cont_param); // continuation parameter
	double xl = conParam;
  
	if (THCM::Instance().GetComm()->MyPID()==0)
    {
		std::string filename="fort.7";
		// this is for the periodic orbit problem, where multiple instances
		// of THCM may be running on different parts of the global comm:
#ifdef HAVE_MPI
		Epetra_MpiComm globcomm(MPI_COMM_WORLD);
		int gpid = globcomm.MyPID();
		int gnp = globcomm.NumProc();
		if (gnp!= THCM::Instance().GetComm()->NumProc())
		{
			std::stringstream ss;
			ss << filename <<"."<<gpid;
			filename = ss.str();
		}
#endif    
		//TODO: find a smart way of handling append or no append
		std::ofstream fort7(filename.c_str(),std::ios::app);
      
		fort7 << itp << " ";
		fort7 << icp << " ";
		fort7 << xl  << " ";
		fort7 << psimmin << " ";
		fort7 << psimmax << " ";
		fort7 << psibmin << " ";
		fort7 << psibmax << std::endl;
		fort7.close();
    }

}

//////////// EpetraExt::ModelEvaluator interface //////////////////////

// get the map of our variable vector [u,v,w,p,T,S]'
Teuchos::RCP<const Epetra_Map> OceanModelEvaluator::get_x_map() const
{
    return THCM::Instance().GetDomain()->GetSolveMap();
}

///////////////////////////////////////////////////////////////////////

// get the map of our 'model response' F(u)
Teuchos::RCP<const Epetra_Map> OceanModelEvaluator::get_f_map() const
{
    return THCM::Instance().GetDomain()->GetSolveMap();
}
  
///////////////////////////////////////////////////////////////////////

// get initial guess (all zeros)
Teuchos::RCP<const Epetra_Vector> OceanModelEvaluator::get_x_init() const
{
    return THCM::Instance().getSolution();
}
    
///////////////////////////////////////////////////////////////////////

// create the Jacobian
Teuchos::RCP<Epetra_Operator> OceanModelEvaluator::create_W() const
{
    return THCM::Instance().getJacobian();
}
    
///////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::InArgs OceanModelEvaluator::createInArgs() const
{
    EpetraExt::ModelEvaluator::InArgsSetup inArgs;
    inArgs.setModelEvalDescription("Ocean Model");
    inArgs.setSupports(IN_ARG_x,    true);
    inArgs.setSupports(IN_ARG_x_dot,true); 
    inArgs.setSupports(IN_ARG_alpha,true); 
    inArgs.setSupports(IN_ARG_beta, true);
    inArgs.setSupports(IN_ARG_t,    true);
    inArgs.set_Np(1); // note: there are actually _NPAR_+2 parameters,
                      // but we store them in a single Epetra_Vector
    return inArgs;            
}

///////////////////////////////////////////////////////////////////////
  
EpetraExt::ModelEvaluator::OutArgs OceanModelEvaluator::createOutArgs() const
{
    EpetraExt::ModelEvaluator::  OutArgsSetup outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(OUT_ARG_f,true);
    outArgs.setSupports(OUT_ARG_W,true);
    // TODO: is this correc? I think I just copied it and left it there...
    outArgs.set_W_properties(
		DerivativeProperties(DERIV_LINEARITY_NONCONST
							 ,DERIV_RANK_FULL,true // supportsAdjoint
			));
    return outArgs;
}
    
///////////////////////////////////////////////////////////////////////

void OceanModelEvaluator::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
	using Teuchos::dyn_cast;
	using Teuchos::rcp_dynamic_cast;
	
	THCM::Instance().startTiming("Compute F and/or Jacobian");
	//
	// Get the input arguments
	//
	const Epetra_Vector &x    = *inArgs.get_x();
	const Epetra_Vector &xdot = *inArgs.get_x_dot();
	double alpha              =  inArgs.get_alpha();
	double beta               =  inArgs.get_beta();
	// there are two ways of setting the time t:  
	// * using inArgs.set_t(),                    
	// * using parameter p0.                      
	// if inArgs.get_t() has the default value 0, 
	// we get it from the parameter vector.       
	double t = inArgs.get_t();
	DEBVAR(t);

	Teuchos::RCP<const Epetra_Vector> p_values = inArgs.get_p(0);
	
	if (t == 0.0)
		t = (*p_values)[0];

	// note: p_values[0] is time: we get it from InArgs instead
  
	// if we do continuation in the exponent, we have to put a different  
	// value into THCM, namely par(exp_cont_par)*10^{s*par(exp)}.         

	double factor =  1.0;
	int     index = -1;    // no param should be multiplied by this factor
	if (cont_param == "Exponent")
    {
		index = THCM::Instance().par2int(exp_cont_param);
		// this has to be set by the user in thcm_params.xml as "Exp. Comp. Parameter"
		if (index < 0)
			Error("Invalid Exp. Cont. Parameter",__FILE__,__LINE__);
		int exp_idx   = THCM::Instance().par2int("Exponent");
		double cont_e = (*p_values)[exp_idx];
		factor = std::pow(10.0, cont_s*cont_e);
    }  

	for (int i = 1; i < p_values->MyLength(); i++)
    {
		std::string label = (*p_names)[i];
		double value      = (*p_values)[i];
		if (index == i)
			value *= factor;
		THCM::Instance().setParameter(label, value);
    }                                
  
	//
	// Get the output arguments
	//
	Teuchos::RCP<Epetra_Vector>    f_out = outArgs.get_f();
	Teuchos::RCP<Epetra_Operator>  W_out = outArgs.get_W();
//  if (showGetInvalidArg_) 
//    {
//    Epetra_Vector *g_out = outArgs.get_g(0).get();
//    }

	// note that setting "Time" in THCM will switch to monthly data instead of 
	// the default (annual mean) as soon as t>0.
	THCM::Instance().setParameter("Time",t);

	bool want_A = ((W_out!=null)&&(beta!=0.0));
	// compute Jacobian and/or RHS
	bool result = THCM::Instance().evaluate(x, f_out, want_A);
  
	if (!result) Error("Error evaluating model!",__FILE__,__LINE__);
  
	//
	// Compute the functions
	//

	if (f_out != null) 
    {
		double n4;
		f_out->NormInf(&n4);
		double *v_values;
		v_values = new double[f_out->MyLength()];

		Teuchos::RCP<Epetra_Vector> vv = THCM::Instance().getSolution();

		FILE  *myfile;
		FILE  *matB;
		myfile = fopen("some.txt"   , "a");
		matB   = fopen("matrixb.txt", "w");
		
		const Epetra_Vector& B = THCM::Instance().DiagB();
		// add diagonal matrix B times xdot

		for (int i = 0; i < f_out->MyLength() ;i++) 
		{			
			(*f_out)[i] = B[i] * xdot[i] + (*f_out)[i];
			(*vv)[i]    = B[i] * xdot[i];
			v_values[i] = B[i] * xdot[i];
			fprintf(matB, "%d    %e    %e  \n", i, B[i], xdot[i]);         
		}
		
		
		delete [] v_values;
		fclose(myfile);
		fclose(matB);

#ifdef HAVE_MPI
		Epetra_MpiComm globcomm(MPI_COMM_WORLD);
		int islemci = globcomm.MyPID();
		int sayi = globcomm.NumProc();
		double n1,n2,n3,n5;
		f_out->NormInf(&n3);
		B.Norm2(&n1);
		xdot.NormInf(&n2);
		vv->NormInf(&n5);

		printf(" number of processor is %d of %d  norm of fout is %e and old fout is %e norm of v is %e  \n",
			   islemci, sayi, n3, n4, n5);
#endif
    }
	
	if (W_out != null) 
    {
		// after the evaluate call above, THCM contains a Teuchos::RCP to the matrix A, which 
		// we hereby extract:
		Teuchos::RCP<Epetra_CrsMatrix> W = THCM::Instance().getJacobian();
   
		printf(" norm of w is %f  \n", W->NormOne());
    
		DEBUG("construct THCM Matrix alpha*B + beta*A");
		DEBVAR(alpha);
		DEBVAR(beta);

		// scale to get beta*A.
		CHECK_ZERO(W->Scale(beta));
  
		// get diagonal of beta*A in a vector
		Teuchos::RCP<Epetra_Vector> diag = Teuchos::rcp(new Epetra_Vector(*get_x_map()));
		CHECK_ZERO(W->ExtractDiagonalCopy(*diag));
  
		// add -alpha*B to diagonal. Note that the sign is reversed as compared to the LOCA interface
		const Epetra_Vector& B = THCM::Instance().DiagB();
		CHECK_ZERO(diag->Update(alpha,B,1.0));
  
		// replace diagonal in Jacobian:
		CHECK_ZERO(W->ReplaceDiagonalValues(*diag));

		// pass it on to the OutArgs:
		Teuchos::RCP<Epetra_CrsMatrix> W_out_crs = 
			rcp_dynamic_cast<Epetra_CrsMatrix>(W_out,true);
		if (W_out_crs.get()!=W.get())
		{
			*W_out_crs = *W;
		}

#ifdef STORE_MATRICES
		std::stringstream ss;
		ss << "Jac_"<<store_step_jac<<".txt";
		HYMLS::MatrixUtils::Dump(*W_out_crs,ss.str());
		store_step_jac++;
#endif

		// compute new scaling
		// TODO: This function uses THCM's internal matrix, which now
		// does not contain beta*A+alpha*B, but only A. Is that a problem?
		// If not, it has already been called in the THCM::evaluate function.
		//THCM::Instance().RecomputeScaling();
	}

	THCM::Instance().stopTiming("Compute F and/or Jacobian",false);
}

///////////////////////////////////////////////////////////////////////


// implementation of NOX::Abstract::PrePostOperator
// (functions that should be called before and after each nonlinear solve
// and nonlinear solver iteration, respectively)

// executed at the start of a call to iterate()
void OceanModelEvaluator::runPreIterate(const NOX::Solver::Generic& solver)
{
	DEBUG("NOX pre-iteration function called");
	if (pres_corr=="pre-iter")
    {
		Error("Pressure Correction not implemented",__FILE__,__LINE__);
    }
}
    
// executed at the end of a call to iterate()
void OceanModelEvaluator::runPostIterate(const NOX::Solver::Generic& solver)
{
	DEBUG("NOX post-iteration function called");
	if (pres_corr=="post-iter")
    {
		Error("Pressure Correction not implemented",__FILE__,__LINE__);
    }
}
        
// executed at the start of a call to solve()
void OceanModelEvaluator::runPreSolve(const NOX::Solver::Generic& solver)
{
	DEBUG("NOX pre-solve function called");
	if (pres_corr=="pre-solve")
    {
		Error("Pressure Correction not implemented",__FILE__,__LINE__);
    }
}
            
// executed at the end of a call to solve()
void OceanModelEvaluator::runPostSolve(const NOX::Solver::Generic& solver)
{
	DEBUG("NOX post-solve function called");
	// invalidate preconditioner after Newton-solve
	if (prec_reuse_policy=="Newton Solve")
    {
		prec_age++;
		if (prec_age>=max_prec_age)
		{
			THCM::Instance().invalidatePreconditioner();
			prec_age=0;
		}
    }
	if (pres_corr=="post-solve")
    {
		Error("Pressure Correction not implemented",__FILE__,__LINE__);
    }
}                              

///////////////////////////////////////////////////////////////////////

// enhanced interface: OceanModel (for NOX/LOCA)

OceanModel::OceanModel(ParameterList& plist, const Teuchos::RCP<LOCA::GlobalData>& globalData,
					   Teuchos::RCP<ParameterList> lsParams)
	: OceanModelEvaluator(plist),
	  LOCA::Epetra::ModelEvaluatorInterface(globalData, rcp(this, false)),
	  backup_filename("IntermediateConfig.txt"), force_backup(false),
	  thcm_output(true), thcm_label(2)
{  
	if (lsParams!=null)
    {

#ifdef NEWDEBUG 
		lsParams->sublist("Block Preconditioner").set("Verbosity",10);
#endif

		Teuchos::RCP<TRIOS::Domain> domainPtr = THCM::Instance().GetDomain();
		Teuchos::RCP<Epetra_CrsMatrix> jacPtr = THCM::Instance().getJacobian();

		std::string prec_type = lsParams->get("User Defined Preconditioner","Block Preconditioner");

		if (prec_type == "Block Preconditioner")
		{
			precPtr = Teuchos::rcp(new TRIOS::BlockPreconditioner(jacPtr,domainPtr,*lsParams));
		}
		else if (prec_type == "HYMLS")
		{
			precPtr = Teuchos::rcp(new TRIOS::HYMLSolver(jacPtr, domainPtr,*lsParams));
		}
		else
		{
			Error("unkown 'User Defined Preconditioner': '"+prec_type+"'.",__FILE__,__LINE__);
		}
		std::cout << *lsParams << std::endl;
    }  
}

///////////////////////////////////////////////////////////////////////

// compute preconditioner, which can then be retrieved by getPreconditioner()
bool OceanModel::computePreconditioner(const Epetra_Vector& x,
									   Epetra_Operator& Prec,
									   ParameterList* p)
{
	DEBUG("enter OceanModel::computePreconditioner");
	if (precPtr == null)
    {
		// no preconditioner parameters passed to constructor
		Error("No Preconditioner available!",__FILE__,__LINE__);
    }
	bool result=precPtr->Initialize();
	if (result==0)
    {
		result=precPtr->Compute();
    }
	DEBUG("leave OceanModel::computePreconditioner");
	return result;
}

///////////////////////////////////////////////////////////////////////

// for XYZT output
void OceanModel::dataForPrintSolution(const int conStep, const int timeStep,
									  const int totalTimeSteps)
{
	// figure out 'what time it is'
	double T = paramList.get("Time Period",0.0);
	double dt = T/totalTimeSteps;
	double t = timeStep*dt;
	THCM::Instance().setParameter("Time",t);
	std::stringstream ss;
	ss << "Config"<<timeStep<<".txt";
	backup_filename = ss.str();
	thcm_label = 2+timeStep;
	INFO("Data for printSolution:");
	INFO(conStep<<" "<<timeStep<<" "<<totalTimeSteps);
	INFO(" backup file \""<<backup_filename<<"\"");
}

//! return the preconditioner operator. Will only be non-null    
//! if you passed a preclist to the constructor. Before using    
//! the preconditioner, computePreconditioner() should be called.
Teuchos::RCP<Epetra_Operator> OceanModel::getPreconditioner()
{
    return Teuchos::rcp_dynamic_cast<Epetra_Operator>(precPtr);
}
  

////////////////
// Call user's own print routine for vector-parameter pair
void OceanModel::printSolution(const Epetra_Vector& x,
							   double conParam) 
{
	bool xmf_out = false; // false: only import solution to grid
  
	if (output_interval>=0)
    {
		if (conParam-last_output>output_interval)
		{
			THCM::Instance().startTiming("Store Solution (XDMF)");
			INFO("Writing Xdmf File...");
			last_output = conParam;
			xmf_out = true; //true: store Xdmf file (if available)
		}
    }

    // either just import solution to grid or
    // also write HDF5/XML files
    xdmfWriter->Store(x,conParam,xmf_out);
    if (xmf_out) THCM::Instance().stopTiming("Store Solution (XDMF)",false);
			
	if ((backup_interval >= 0) || force_backup)
    {      
		if ((conParam - last_backup > backup_interval) ||
			force_backup ||
			(conParam == last_backup) )
		{
			//three cases where we write a backup:
			// - some time has passed since last backup (backup_interval)
			// - the user forces us to (i.e. at the end of a continuation, force_backup)
			// - last_backup indicates that this function is being called repeatedly,   
			// - which probably means that we're in LOCA XYZT mode
			THCM::Instance().startTiming("Store Solution (Backup)");
			INFO("Writing Backup at param value "<<conParam<<"...");
			WriteConfiguration(backup_filename,*pVector,x);
			last_backup = conParam;
			THCM::Instance().stopTiming("Store Solution (Backup)",false);

			if (thcm_output)
			{
				// write solution in native THCM format

				int filename = 3; // write to 'fort.3'
				int label    = thcm_label;

				// gather the solution on the root process (pid 0):
				Teuchos::RCP<Epetra_MultiVector> fullSol = HYMLS::MatrixUtils::Gather(x,0);

				if (THCM::Instance().GetComm()->MyPID() == 0)
				{
					int ndim = fullSol->GlobalLength();
					double *solutionbig;

					(*fullSol)(0)->ExtractView(&solutionbig);
					(*info) << " Store solution in THCM format..." << std::endl;
					if (label == 2)
					{
						INFO("writing data, label = " << label);
						FNAME(write_data)(solutionbig, &filename, &label);
					}
					else 
					{
						INFO("appending data, label = "<<label);
						FNAME(append_data)(solutionbig, &filename, &label);
					}
					//fort.15...
					F90NAME(m_global,compute_flux)(solutionbig);
					(*info) << "done!"<<std::endl;
				}
			}
      
		}
	}
    
	// in every step, we compute the meridional and barotropic streamfunctions
	// and store their maximum in fort.7 (in the old THCM format)
	// At this point, 'grid' contains the complete solution in 3D array format
	this->Monitor(conParam);

	// invalidate preconditioner after Time-/Continuation step
	if (prec_reuse_policy=="Step")
    {
		prec_age++;
		if (prec_age>=max_prec_age)
		{
			THCM::Instance().invalidatePreconditioner();
			prec_age=0;
		}
    }
  
}


///////////////////////////////////////////////////////////////////////
