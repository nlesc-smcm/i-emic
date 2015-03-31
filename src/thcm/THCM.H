/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef THCM_H
#define THCM_H

#include "Singleton.H"
#include "Epetra_Object.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Utils.H"

namespace TRIOS
{
	class Domain;
}

class Epetra_Comm;
class Epetra_Map;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_BlockMap;
class Epetra_MultiVector;
  
namespace EpetraExt
{
	class MultiComm;
	class BlockVector;
}
  
//!                                                              
//! This class provides the interface between THCM and Trilinos. 
//! It uses a domain-object (TRIOS::Domain) to decompose            
//! the domain into subdomains and generate two Epetra_Maps.     
//!  The 'Assembly' or 'Local'                                   
//! map includes a layer of ghost-nodes between subdomains, i.e. 
//!                                                              
//! \verbatim                                                    
//!  _______________........                                     
//! |o o o : o* o* | * * * :                                     
//! |o o o : o* o* | * * * :                                     
//! |______:_______|........                                     
//!                                                              
//!                                                              
//! \endverbatim                                                 
//!                                                              
//! whereas the Solve-map is a bijection from local              
//! to global indices. The Jacobian and rhs/sol vectors are      
//! duplicated in memory, the 'Assembly' version used internally 
//! to interface with THCM and the 'Solve' version used for the  
//! Trilinos solvers.                                            
//! THCM operates on the extended (assembly) subdomain, so that  
//! it can be mostly oblivious of the fact that the program is   
//! running in parallel. A minimum of pre- and post processing   
//! on the 'global' domain is necessary, though.                 
//!                                                              
//! rhs is the result of a function evaluation, rhs = f(u), if   
//! the model is written as Bdu/dt + f(u) = 0. The Jacobian      
//!  is A=df/du. To avoid confusion it is                        
//! recommended to use class OceanModel to 'evaluate the model'  
//!                                                              
//! This class is implemented as a 'Singleton', which means that 
//! there is only one instance allowed at a time. This is reason-
//! able because the THCM fortran data structures can only       
//! support one 'instance' of THCM at a time. Other classes can  
//! access this object (once it has been constructed) by a call  
//! to THCM::Instance().                                         
//! When attempting to create another instance, the old one will 
//! be deleted and a warning is issued.                          
//!                                                              
class THCM : 
	public Singleton<THCM>, public Epetra_Object
{

public:

	//! Constructor  
	THCM(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Comm> comm);
  
	//! Destructor

	/*! \note: currently the Fortran data structures are not 
	  completely deallocated, so it is dangerous to 
	  delete and re-allocate the THCM object during
	  a single run (this would create a memory hole 
	  or cause Fortran errors!)
	*/  
	virtual ~THCM();

	//! compute the rhs vector and/or the jacobian.                        
  
	/*! The rhs is computed and returned in *rhsVector if it is not null.  
      Note that the sign of the rhs is reversed as compared to THCM.     
                                                                         
      If computeJac=true the Jacobian is computed and can be obtained    
      by calling getJacobian(). The Jacobian in THCM is A-sigma*B, but   
      we keep igma set to 0. Use diagB to access the B matrix.
	*/
	bool evaluate (const Epetra_Vector& solnVector,
				   Teuchos::RCP<Epetra_Vector> rhsVector, bool computeJac=false);
  
	//! only recompute the diagonal matrix B
  
	/*! the matrix B is used by THCM to 'switch off' some equations.
      It is diagonal and the diagonal entries can be accessed by  
      the function diagB().                                       
      We use the opposite sign as compared to THCM.
      The entries of B are (see assemble.f::fillcolB and          
      usrc.F90::matrix): 
	  Ro for u,v    
	  0  for w,p    
	  1 for T,S     
	*/
	void evaluateB();

	//! Returns initial guess (Global/Solve form)
	Teuchos::RCP<Epetra_Vector> getSolution();
  
	//! returns the Jacobian matrix (Global/Solve form)
	Teuchos::RCP<Epetra_CrsMatrix> getJacobian();
  
	//! get null space
	Teuchos::RCP<const Epetra_MultiVector> getNullSpace();
    
	//! get row scaling vector
	Teuchos::RCP<Epetra_Vector> getRowScaling() {return row_scaling;}

	//! get column scaling vector
	Teuchos::RCP<Epetra_Vector> getColScaling() {return col_scaling;}

	//! start a timer for a given part of the program
	void startTiming(std::string id);

	//! stop timer for a given part of the program
	void stopTiming(std::string id,bool print=false);
  
	//! print timer status for all timed program parts
	void printTiming(std::ostream&);
    
	//! set the pressure in one specific point to zero by subtracting
	//! that point's value from all pressure points.
	void normalizePressure(Epetra_Vector& soln) const;
  
	//! Set a bifurcation parameter in the application physics
	bool setParameter(std::string label, double value);

	//! Get a bifurcation parameter from the application physics
	bool getParameter(std::string label, double& value);

	//! Write all THCM params in fort.7
	bool writeParams();

	//! returns the domain decomposition object
	Teuchos::RCP<TRIOS::Domain> GetDomain() {return domain;}
  
	//! get the diagonal matrix B
	Epetra_Vector& DiagB(){return *diagB;}

	//! return the communicator object
	Teuchos::RCP<Epetra_Comm> GetComm(){return Comm;}
  
	//! get the SRES parameter (non-restoring salt condition)
	bool getSRES() const {return sres;}

	//! ite/its=0 T/S forcing from data (Levitus)
	int ite,its;
  
	//! use internal temperature and salinity forcing
	bool internal_forcing;
  
	//! wind forcing method (0: data (trenberth), 1: zonally averaged, 2: idealized)
	int iza;

	//! get the global index of the row with the integral condition (if sres==0)
	int getRowIntCon() const {return rowintcon_;}
  
	//! set the THCM flag 'vmix_fix' to 0 or 1
	//! note: vmix_fix doesn't have to be set if you have vmix_flag=1 in mix_imp.f 
	//! (recommended). If you REALLY want vmix_flag=2, call this function in the NOX
	//! pre-solve function of class OceanModel.
	void fixMixing(int value);

	//! \name get physical global domain bounds
	//@{
  
	//!
	inline double xMin() const {return xmin;}
	//!
	inline double xMax() const {return xmax;}
	//!
	inline double yMin() const {return ymin;}
	//!
	inline double yMax() const {return ymax;}
	//!
	inline double hDim() const {return hdim;}
  
	//@}

	//! \name these are used by the block preconditioner:
	//!@
	inline double AlphaT() const {return alphaT;}
	inline double AlphaS() const {return alphaS;}
	inline double rhoMixing() const {return rho_mixing;}
	//@}

	//! convert parameter name to integer (i.e. "Combined Forcing" => 19)
	int par2int(std::string par);

	// convert parameter index to std::string (i.e. 19 => "Combined Forcing")
	std::string int2par(int ind);

private:

	//! object for domain-decomposition:
	Teuchos::RCP<TRIOS::Domain> domain;

	//!\name maps and import/exprt objects for distributed data structures
	//! (obtained from the domain-object)
    
	//!@{
  
	//! overlapping map for 'THCM objects': 
	Teuchos::RCP<Epetra_Map> AssemblyMap;

	//! non-overlapping map for Trilinos objects: 
	Teuchos::RCP<Epetra_Map> StandardMap;  
  
	//! non-overlapping load-balanced map for Trilinos objects 
	//! (may contain non-rectangular subdomains)               
	Teuchos::RCP<Epetra_Map> SolveMap;
  
	//!@}
  
	//! used only to define vector format (i.e. map), I think
	Teuchos::RCP<Epetra_Vector> initialSolution;
	//! used to import current approximation to THCM:
	Teuchos::RCP<Epetra_Vector> localSol;
	//! used to export computed RHS vector to Trilinos:
	Teuchos::RCP<Epetra_Vector> localRhs;
  
	//! the diagonal matrix B stored as a vector
	Teuchos::RCP<Epetra_Vector> localDiagB, diagB;
  
	//! Jacobian in globally assembled and load-balanced
	Teuchos::RCP<Epetra_CrsMatrix> Jac;

	//! Jacobian based on standard subdomains
	Teuchos::RCP<Epetra_CrsMatrix> localJac;
  
	//! static graph which defines the maximal pattern of non-zero entries.
	//! the matrix Jac is based on this graph, so depending on forcing and 
	//! convective adjustment, it may contain some zeros
	Teuchos::RCP<Epetra_CrsGraph> MatrixGraph;

	//! same crs graph based on standard map
	Teuchos::RCP<Epetra_CrsGraph> localMatrixGraph;
  
	//! global communicator for 4D things
	Teuchos::RCP<EpetraExt::MultiComm> xyzt_comm;
  
	//! type of scaling to be done
	std::string scaling_type;
  
	//! vectors storing the THCM scaling
	Teuchos::RCP<Epetra_Vector> row_scaling, col_scaling,local_row_scaling,local_col_scaling;

	//! vectors storing the THCM scaling (4D case)
	Teuchos::RCP<EpetraExt::BlockVector> xyzt_row_scaling, xyzt_col_scaling;

	//! additional row scaling for T and S (for scaling type THCM/TS)
	Teuchos::RCP<Epetra_Vector> row_scaling_TS;

	//! (MPI) communicator
	Teuchos::RCP<Epetra_Comm> Comm; 
  
	//! nullspace (p-vectors)
	Teuchos::RCP<Epetra_MultiVector> nullSpace;
  
	//! global shared parameter list
	Teuchos::ParameterList& paramList;
  
	//! timer status list
	Teuchos::ParameterList timerList;

	//! \name the matrix in THCM (CSR-format)
   
	/*! Memory for the Jacobian on the subdomain is allocated by the
	  C++ code. In THCM, pointers are set to these locations so that  
	  the Fortran code fills these arrays directly. Afterwards the    
	  matrix is copied to the 'Jac' member, ignoring 'ghost' rows.    
	*/

	//!@{
  
	//! row pointer
	int* begA; 
	//! column indices
	int* jcoA; 
	//! values
	double* coA; 
	//! values of the diagonal matrix B
	double* coB;
	//!@}
  
	//! global grid dimensions
	int n,m,l,la; 

	//! physical domain bounds (global long./lat./depth)
	double xmin,xmax,ymin,ymax,hdim;
  
	//! periodic domain in x-direction?
	bool periodic;

	//! inhomogeneous mixing?
	int ih;

	//! mixing?
	int vmix_GLB;
  
	//! parameters for linear equation of state
	double alphaT, alphaS;
  
	//! taper
	int tap;
  
	//! use fred's formulation of conv. adjustment
	bool rho_mixing;
  
	//! factors sigma for the Jacobian A+sigma*B
	double sigmaUVTS,sigmaWP;
  
	//! sres=0: non-restoring salt forcing => integral condition in A and f
	int sres;

	//! tres=0: non-restoring temperature forcing => integral condition in A and f
	int tres;
  
	//! which row is replaced by integral condition (global index of last row)
	int rowintcon_;
  
	//! vector with coefficients for integral condition (if sres=0)
	Teuchos::RCP<Epetra_Vector> intcond_coeff;
  
	//! pressure points where equation is replaced by P=0 (-1 means none)
	int rowPfix1, rowPfix2;

	//! asks THCM to recompute scaling vectors
	void RecomputeScaling(void);

	//! distribute land array after global initialization
	Teuchos::RCP<Epetra_IntVector> distributeLandMask(Teuchos::RCP<Epetra_IntVector> landm_glob);
  
	//! implement integral condition for S in Jacobian and B-matrix
	void intcond_S(Epetra_CrsMatrix& A, Epetra_Vector& B);

	//! implement dirichlet values P=0 in cells rowPfix1/2 (if >=0)
	void fixPressurePoints(Epetra_CrsMatrix& A, Epetra_Vector& B);

	//! this subroutine defines the maximal matrix graph (pattern of nonzeros in the jacobian
	//! if convective adjustment occurs in all cells).
	Teuchos::RCP<Epetra_CrsGraph> CreateMaximalGraph();

    
	//! private function used by CreateMaximalGraph()
	void insert_graph_entry(int* indices, int& pos, 
							int i, int j, int k, int var,
							int N, int M, int L) const;

	//! read monthly levitus fields, distribute, pass back to local THCM
	//! note that wind is read in by each process and interpolated, but 
	//! for some reason that doesn't work for T and S
	void SetupMonthlyForcing();

	//! setup load-balancing (optionally reads weights or creates them
	//! using thcm_utils, which only looks at land mask at this stage)
	//! This alters the behaviour of the TRIOS::Domain class, which will
	//! provide a different SolveMap after the call.                  
	void SetupLoadBalancing();

};

// this is called from forcing.F90 to compute some global integral if SRES=0
/*
  extern "C" {
  void thcm_forcing_integral(double* qfun2, int* landm, double* fsint);
  };
*/
#endif
