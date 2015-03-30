/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef OCEANMODEL_H
#define OCEANMODEL_H

#include "THCM.H"
#include "EpetraExt_ModelEvaluator.h"
#include "LOCA_Epetra_ModelEvaluatorInterface.H"
#include "NOX.H"

#include "Epetra_Vector.h" 
#include "LOCA_Epetra.H"
#include <fstream>
#include <iostream>
#include <sstream>

#include "Teuchos_Array.hpp"

//Forward Declaration

class OceanOutputXdmf;

namespace LOCA
{
	class ParameterVector;
}

class OceanGrid;


//! class to interface between Trilinos and THCM

/*! This class implements the EpetraExt ModelEvaluator class 
  (which is also compatible with the Thyra package).
    
  The class interfaces THCM via the Singleton class THCM (THCM.H).
    
  The ModelEvaluator interface provides Jacobian A, Residual f,            
  matrix alpha*B + beta*A etc, through the evalModel(...) function.        
  You can also use class OceanModel (see below), which implements          
  the full LOCA interface called "TimeDependent" and also the              
  NOX::Interface::Preconditioner for integrating our block preconditioner. 
                                                                             
  The ModelEvaluator class uses the following formulation of the model:    
  F(xdot,x,t) := B xdot + f(x,t) = 0, W = alpha*B + beta*df/dx             
*/
class OceanModelEvaluator :
	public EpetraExt::ModelEvaluator, 
	public NOX::Abstract::PrePostOperator,
	public Epetra_Object
{

public:

	//! constructor: plist is the "THCM" sublist
	OceanModelEvaluator(Teuchos::ParameterList& plist);

	//! destructor
	~OceanModelEvaluator();

//! \name EpetraExt/Thyra interface                              

//@{

	//! get the map of the variable vector [u,v,w,p,T,S]'
	Teuchos::RCP<const Epetra_Map> get_x_map() const;

	//! get the map of the result vector from a function evaluation f(x) (same as the x_map)
	Teuchos::RCP<const Epetra_Map> get_f_map() const;

	//! starting guess
	Teuchos::RCP<const Epetra_Vector> get_x_init() const;

	//! get the Jacobian
	Teuchos::RCP<Epetra_Operator> create_W() const;

	//! get the names of parameters
	Teuchos::RCP<const Teuchos::Array<std::string> >
	get_p_names(int l) const {return p_names;}
  
	//! get initial values for parameters
	Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const
		{return p_init;}

	//! get map for parameters
	Teuchos::RCP<const Epetra_Map> get_p_map(int l) const
		{return p_map;}

	//!
	InArgs createInArgs() const;

	//!
	OutArgs createOutArgs() const;

	//! compute F, A, alpha*B+beta*A, ...
	//! See Thyra_ModelEvaluator for info on how to use this
	void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

//@}


	//!\name pre- and post operations during nonlinear solve

	/*! we may use these during time integration for things like I/O, 
      fixing the pressure, recomputing the preconditioner etc.
      
      Note: since the argument is 'const' I don't really see how these functions can be
	  of any use, how can we modify the current approximation etc.?
	*/
  
	//@{

	//! User defined method that will be executed at the start of a call to NOX::Solver::Generic::iterate().
	void runPreIterate(const NOX::Solver::Generic& solver);

	//! User defined method that will be executed at the end of a call to NOX::Solver::Generic::iterate().
	void runPostIterate(const NOX::Solver::Generic& solver);

	//! User defined method that will be executed at the start of a call to NOX::Solver::Generic::solve().
	void runPreSolve(const NOX::Solver::Generic& solver);

	//! User defined method that will be executed at the end of a call to NOX::Solver::Generic::solve().
	void runPostSolve(const NOX::Solver::Generic& solver);

	//@}



	//! Return the parameters from the application.
	Teuchos::RCP<LOCA::ParameterVector> getParameterVector()  
		{ 
			return THCM::Instance().getParameterVector();
		}

	//! Return the Jacobian from the application.
	Teuchos::RCP<Epetra_CrsMatrix> getJacobian()  
		{ 
			return THCM::Instance().getJacobian();
		}

	//! Return the current solution from the application.
	Teuchos::RCP<Epetra_Vector> getSolution()
		{
			return THCM::Instance().getSolution();
		}

	//! Return the null space (two singular vectors of p)
	Teuchos::RCP<const Epetra_MultiVector> getNullSpace()
		{
			return THCM::Instance().getNullSpace();
		}

	//! Return scaling object (may be null)
	Teuchos::RCP<NOX::Epetra::Scaling> getScaling()
		{
			return THCM::Instance().getScaling();
		}


	//! apply pressure correction
	//! 'never'       
	//! 'pre-iter'    
	//! 'post-iter'   
	//! 'pre-solve'   
	//! 'post-solve'  
	//!               
	//! NOTE: this is NOT IMPLEMENTED!
	void setPressureCorrection(std::string when)
		{
			if ((when=="pre-iter")
				||(when=="post-iter")
				||(when=="pre-solve")
				||(when=="post-solve"))
			{
				pres_corr=when;    
			}
			else
			{
				Error("Invalid PresureCorrection Scheme",__FILE__,__LINE__);
			}
		}

	//! get pointer to grid-representation, which always contains the solution 
	//! from the last call to 'printSolution'
	Teuchos::RCP<OceanGrid> getGrid() { return gridPtr; }
      
	//! read parameters and solution from an ASCII file
	Teuchos::RCP<Epetra_Vector> ReadConfiguration(std::string filename, LOCA::ParameterVector& params);

	//! write parameters and solution to an ASCII file
	void WriteConfiguration(std::string filename, const LOCA::ParameterVector& pvec,
							const Epetra_Vector& soln);

protected:  
  
	//! this is a grid-representation of the solution, which we share with the XdmfWriter
	//! Even if Xdmf is not available, it can be used to compute statistics like the max.
	//! streamfunction etc.

	Teuchos::RCP<OceanGrid> gridPtr;
    
	//! parameters obtained from the continuation process
	Teuchos::RCP<LOCA::ParameterVector> pVector;
  
	//! parameter names (this is redundant, but needed for the ModelEvaluator)
	Teuchos::RCP<Teuchos::Array<std::string> > p_names;

	//! parameter values (this is redundant, but needed for the ModelEvaluator)
	Teuchos::RCP<Epetra_Vector> p_init;
  
	//! simple replicated map for the p_values
	Teuchos::RCP<Epetra_Map> p_map;
  
#ifdef STORE_MATRICES
	//! a debugging facility: the first 'STORE_MATRICES'
	// Jacobians are written to numbered files
	mutable int store_step_jac, store_step_rhs;
#endif

	//! global shared parameter list
	Teuchos::ParameterList& paramList;

	//! continuation parameter name
	std::string cont_param;

	//! actual continuation parameter name if cont_param=="Exponent" (see cont_s)
	std::string exp_cont_param;
  
	//! If cont_param=="Exponent", we do                          
	//! continuation in e, p0*10^(s*e), instead of cont_param.    
	//! This is s in the formula above, set by the option         
	//!  "Continuation in Exponent"   (double)                    
	//! in thcm_params.xml. p0 is the initial value given to      
	//! par, which is the "Exp. Cont. Parameter" (exp_cont_param) 
	//! from thcm_params.xml. By default s=1.0 (continuation in   
	//! exp_cont_param). Use s=-1.0 to continue backward,         
	//! for instance in the horizontal friction                   
	//! ("Horizontal Ekman-Number").                              
	double cont_s;

	//! preconditioner reuse policy ("None","Newton Solve","Step")
	std::string prec_reuse_policy;

	//! monitor the age of the preconditioner in time/cont. steps
	int prec_age, max_prec_age;

	//! index of reference cell (rowintcon in THCM, here it is global)
	int refcell_ind;
  
	//! can be set by the user to to implement time stepping
  
	//! backup interval
	double backup_interval;

	//! when was the last backup written?
	double last_backup;
  
	//! time series output interval
	double output_interval;

	//! when was the last series output written?
	double last_output;

	//! Object for writing Xdmf data files at regular intervals
	Teuchos::RCP<OceanOutputXdmf> xdmfWriter;
  
	//! apply pressure correction during or after Newton?
	std::string pres_corr;

	//! compute meridional and barotropic streamfunction based on solution in 'grid',
	//! and store the THCM formatted file 'fort.7'. conParam is the current value of the
	//! continuation parameter (or time), whose name is stored in cont_param
	void Monitor(double conParam);

private:

	//! read parameter entry from ascii stream
	void read_parameter_entry(Teuchos::RCP<std::istream>, std::string&, double&);
    
};

//! this class is an adaptor for LOCA: it inherits everything      
//! from OceanModelEvaluator and ModelEvaluatorInterface to form a 
//! complete Interface::TimeDependent.                             
//! This interface uses the following formulation of the model:    
//! B dx/dt = f(x,t), W = alpha*df/dx + beta*B                     
//! In addition, the 'Preconditioner' interface is implemented for NOX. The 
//! preconditioner supplied by this class is the block preconditioner imple-
//! mented in OceanPreconditioner.H/C. To obtain the prec operator, use     
//! getPreconditioner(), passing in the "Block Preconditioner" sublist.     
    
class OceanModel : 
	public OceanModelEvaluator,
	public LOCA::Epetra::ModelEvaluatorInterface
{
public:
	
	//! constructor:
	//! if preclist!=null, a preconditioner will be          
	//! available through computePreconditioner() and getPreconditioner().
	//! plist is the "THCM" sublist, preclist is the                      
	//! "Linear Solver"->"Block Preconditioner" sublist.                  
    OceanModel(Teuchos::ParameterList& plist, const Teuchos::RCP<LOCA::GlobalData>& globalData,
			   Teuchos::RCP<Teuchos::ParameterList> preclist=Teuchos::null);

//! \name NOX/LOCA interface 

//!@{
  
/*!
  \brief Provides data to application for output files.

  This routine is called from Interface::xyzt::printSolution() just
  before the call to Interface::Required::printSolution(x,param),
  and gives the application some indices that can be used for
  creating a unique name/index for the output files.
*/
	void dataForPrintSolution(const int conStep, const int timeStep,
							  const int totalTimeSteps);
  
  
	//! Call user's own print routine for vector-parameter pair
	void printSolution(const Epetra_Vector& x_,
					   double conParam);

	//!@}
                           
	//! close output files (happens automatically during destruction but
	//! it may have to be done by hand before calling MPI_Finalize)
	//! This means that the solution is no longer written in printSolution().
	void finishOutput() {xdmfWriter = Teuchos::null; output_interval=-1;}

//!\name NOX Preconditioner interface

//!@{

//! Compute preconditioner \f$M\f$.

//! note: this may cause an error if you didn't give a preclist 
//!       to the constructor.                                   
	virtual bool computePreconditioner(const Epetra_Vector& x,
									   Epetra_Operator& Prec,
									   Teuchos::ParameterList* p = NULL);                                                                                              
                                                                                              
//!@}

	//! return the preconditioner operator. Will only be non-null    
	//! if you passed a preclist to the constructor. Before using    
	//! the preconditioner, computePreconditioner() should be called.
	Teuchos::RCP<Epetra_Operator> getPreconditioner();
  
	//! always backup, regardless of parameter value (mainly used for saving finl config)
	void setForceBackup(bool value){force_backup=value;}

	//! also write THCM's native fort.* files (mainly used for saving finl config)
	void setThcmOutput(bool value){thcm_output=value;}
  
	//! get name of continuation parameter
	std::string getContinuationParameterName(){return cont_param;}

	//! get value of continuation parameter
	double getContinuationParameterValue() 
		{
			double value=0.0;
			try {
				value=pVector->getValue(cont_param);
			} catch (...) {std::cerr<<"Warning: caught an exception ("<<__FILE__<<", line "<<__LINE__<<")\n";}
			return value;
		}

protected:

	//! filename of backup file (by default it's "IntermediateConfig.txt", this only
	//! changes in XYZT mode)
	std::string backup_filename;
  
	//! this is usually false, but if it's true a backup file is always written
	bool force_backup;
  
	//! this is usually false, but if it's true, the 'fort.3' and 'fort.44' files 
	//! are written in addition to the standard Trilinos-THCM output format.      
	bool thcm_output;
  
	//! append to fort.3 or create new one. 2: create, >2: append. This label is stored in   fort.3
	int thcm_label;

	//! custom preconditioner, constructed when you first call computePreconditioner
	Teuchos::RCP<Ifpack_Preconditioner> precPtr;
      
};

#endif
