/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef TRIOS_SADDLEPOINT_H
#define TRIOS_SADDLEPOINT_H

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

// typedef'd Teuchos pointers

class AztecOO;

namespace TRIOS {
	
#ifdef HAVE_ZOLTAN
	class Repart;
#endif

//////////////////////////////////////////////////////////////////////////////////////////

//! general Saddlepoint matrix class

/*! implements the block matrix
\begin{verbatim}

      | A11 A12 |
Spp = | A21  0  |

\end{verbatim

*/
class SaddlepointMatrix : public Epetra_Operator
  {
  
 protected:
   
   //! basic constructor without blocks
   SaddlepointMatrix(Teuchos::RCP<Epetra_Comm> comm);
   
   //! set blocks
   void SetBlocks(Teuchos::RCP<Epetra_CrsMatrix> a11, Teuchos::RCP<Epetra_CrsMatrix> a12,
				  Teuchos::RCP<Epetra_CrsMatrix> a21);
   
  
 public:
//-----//

friend class SppSimplePrec;
  
  //! constructor
  SaddlepointMatrix(Teuchos::RCP<Epetra_CrsMatrix> a11, Teuchos::RCP<Epetra_CrsMatrix> a12, 
            Teuchos::RCP<Epetra_CrsMatrix> a21, Teuchos::RCP<Epetra_Comm> comm);
      
  //! destructor
  virtual ~SaddlepointMatrix();

//!\name  Atribute set methods
//!@{

  //! If set true, transpose of this operator will be applied 
  //! not yet implemented, returns -1.
  int 	SetUseTranspose (bool UseTranspose)
    {return -1;}

//!@}

//!\name Mathematical functions
//!@{

  //! apply operator Y=Op*X
  int 	Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

  //! apply inverse operator (n/a)
  virtual int ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

  //! Returns the infinity norm of the global matrix.
  double NormInf () const {return normInf;}

//!@}

  
//!\name Atribute access functions
//!@{

  //!Returns a character string describing the operator.
  const char * 	Label () const {return label_.c_str();}
  
  //! Returns the current UseTranspose setting.
  bool 	UseTranspose () const {return false;}

  //! Returns true if the this object can provide an approximate Inf-norm, false otherwise.
  bool 	HasNormInf () const {return true;}
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & 	Comm () const {return *comm_;}
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & 	OperatorDomainMap () const 
    {
    return *domainMap;
    }
    
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & 	OperatorRangeMap () const 
    {
    return *rangeMap;
    }

  //! Returns the Epetra_Map object associated with the domain of this operator.
  Teuchos::RCP<Epetra_Map> GetDomainMap () 
    {
    return domainMap;
    }
    
  //! Returns the Epetra_Map object associated with the range of this operator.
  Teuchos::RCP<Epetra_Map> GetRangeMap ()
    {
    return rangeMap;
    }
  //!@}

  //! \name subblock accessor functions
  //!@{

  inline const Epetra_CrsMatrix& A11() const { return *A11_; }     
  inline Epetra_CrsMatrix& A11() { return *A11_; }     
  inline Teuchos::RCP<Epetra_CrsMatrix> GetA11() { return A11_; }     
  inline const Epetra_CrsMatrix& A12() const { return *A12_; }     
  inline Epetra_CrsMatrix& A12() { return *A12_; }     
  inline Teuchos::RCP<Epetra_CrsMatrix> GetA12() { return A12_; }     
  inline const Epetra_CrsMatrix& A21() const { return *A21_; }     
  inline Epetra_CrsMatrix& A21() { return *A21_; }     
  inline Teuchos::RCP<Epetra_CrsMatrix> GetA21() { return A21_; }     
  //!@}
  
 protected:
//--------//

  //! infinity norm is stored
  double normInf;

  //! label identifying this class
  std::string label_;

  //! domain map (map for vectors x in Ax=y)
  Teuchos::RCP<Epetra_Map> domainMap;
  
  //! range map (map for vectors y in Ax=y)
  Teuchos::RCP<Epetra_Map> rangeMap;

  //! subblocks
  Teuchos::RCP<Epetra_CrsMatrix> A11_,A12_,A21_;

  
  //! communicator
  Teuchos::RCP<Epetra_Comm> comm_;
  
protected:

  //! recompute infty-norm
  void recompute_normInf();
  
  //! apply operator to pre-split vector
  int Apply(const Epetra_Vector& x1, const Epetra_Vector& x2,
                   Epetra_Vector& y1,       Epetra_Vector& y2) const;
  
  };

//! operator class for the depth-averaged saddlepoint problem

//! This is a matrix class that represents our 2x2
//! block matrix (the depth-av'd  Saddlepoint problem,
//! SPP). It only supports multiplication by a vector
//! and construction from precomputed blocks.
//! \verbatim
//!                 ___                         
//!        | Auv    Guv |                       
//! Spp =  | ___        |                       
//!        | Duv    0   |                       
//!                                             
//! \endverbatim                                
//!                                             
//! \note: depth-averaging the velocity is now  
//!       no longer supported                   
//! \note: ___              ___                 
//!        Guv = Guv*Mzp1', Duv = Mzp2*Duv      
class SppDAMatrix : public SaddlepointMatrix
  {
  
 public:
//-----//

  
  //! constructor
  SppDAMatrix(Epetra_CrsMatrix& Mzp1, Epetra_CrsMatrix& Mzp2, 
              Epetra_CrsMatrix& Auv, 
              Epetra_CrsMatrix& Guv, Epetra_CrsMatrix& Duv, Teuchos::RCP<Epetra_Comm> comm_);
              
      
  //! destructor
  virtual ~SppDAMatrix();

  //! set new matrix Auv only (Duv and Guv are constant)
  void Update(Epetra_CrsMatrix& Auv);

  //! \name submatrix access
  //!@{
  
  //!
  inline const Epetra_CrsMatrix& mAuv() const {return A11();}
  //!
  inline const Epetra_CrsMatrix& mDuv() const {return A21();}
  //!
  inline const Epetra_CrsMatrix& mGuv() const {return A12();}
  //@}
  
  }; // end of class SppDAMatrix
  

//////////////////////////////////////////////////////////////////////


//! SIMPLE-type preconditioner for the saddlepoint problem

/*! implements the modified Simple, SL, and Simpler       
    preconditioners. They are modified in the sense that  
    we use a 2x2-blockdiagonal to approximate the inverse 
    of the 'A' matrix in order to catch the Coriolis force
                                                          
\begin{verbatim}                                          
                                                          
  to precondition                                         
                                                          
    A11 A12                                               
    A21  0                                                
                                                          
we use                                                    
                                                          
|A11 A11*D\A12|                                           
|A21    0     |                                           
                                                          
=                                                         
                                                          
|A11   0|  |I   D\A12    |                                
|A21   I|  |0 -A21*D\A12 |                                
                                                          
where D is the 2x2 blockdiagonal of A11.                  
                                                          
\end{verbatim}                                            
                                                          
*/
class SppSimplePrec : public Epetra_Operator
  {

  public:

  //! constructor                                              
  //!                                                          
  //! Spp: the saddlepoint problem to be preconditioned        
  //! scheme should be                                         
  //! - "SI" (Simple)                                          
  //! - "SL" (transposed Simple)                               
  //! - "SR" (Simpler or symmetrized Simple)                   
  //!                                                          
  //! zero_init: use 0 initial guess for inner solves          
  //!                                                          
  //! A11Solver|PrecList are used to construct a solver/       
  //!     preconditioner for A11 _unless_ A11Solver!=null      
  //!                                                          
  //! ChatXYZList are used likewise for the Schur-complement   
  //! SpaIList: this list contains the settings for the sparse 
  //! approximate inverse of A11 which should be used when     
  //! constructing the approximate Schur complement Chat.      
  //! it is a standard preconditioner list but only the choices
  //! "Block Jacobi" and "ParaSails" make sense and will work. 
  //! "Block Jacobi" selects the 2x2 block diagonal of A11 as  
  //! approximation whereas ParaSails can be used to build     
  //! better overall preconditioners.                          
	  SppSimplePrec(Teuchos::RCP<SaddlepointMatrix> Spp, Teuchos::ParameterList& params_,
					Teuchos::RCP<Epetra_Comm> comm_, Teuchos::RCP<AztecOO> A11Solver_,
					Teuchos::RCP<Epetra_Operator> A11Precond_,
					bool zero_init_=true);      
                
  //! Destructor
  virtual ~SppSimplePrec();
      
  //! Set transpose (n/a)
  
  //! The user should use the 'scheme' parameter passed to 
  //! the constructor to determine how the preconditioner is to be
  //! applied. 
  int SetUseTranspose(bool UseTranspose){return -1;}

  //! Apply preconditioner operator (n/a)
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
      
  //! Apply preconditioner operator inverse
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
      
  //! Computing infinity norm 
  double NormInf() const;
      
  //! Label
  const char* Label() const {return label_.c_str();}
      
  //! Transpose
  bool UseTranspose() const {return false;}
      
  //! Have norm-inf
  bool HasNormInf() const {return false;}
      
      //! get communicator
      
      /*! 
       Returns a pointer to the Epetra_Comm communicator associated 
       with this operator.
       */
      const Epetra_Comm& Comm() const {return *comm;}
      
      /*! 
       * \brief Returns the Epetra_Map object associated with the domain of 
       * this operator.
       */
      const Epetra_Map& OperatorDomainMap() const
        { 
        return *domainMap;
        }
      
      /*! 
       * \brief Returns the Epetra_Map object associated with the range of 
       * this operator.
       */
      virtual const Epetra_Map& OperatorRangeMap() const
        { 
        return *rangeMap;
        }
      
      
    protected:

      //! Label for this operator.
      std::string label_;

      //! communicator
      Teuchos::RCP<Epetra_Comm> comm;
      
      //! use zero starting guess for linear system solves?
      bool zero_init;

      //! Simple-scheme used: "SI","SL" or "SR"
      std::string scheme;
      
      //! domain- and rangemaps
      Teuchos::RCP<Epetra_Map> domainMap, rangeMap;

      //! The depth-averaged saddlepoint operator
      Teuchos::RCP<SaddlepointMatrix> Spp;
      
      //! 2x2 block diagonal of MAuv
      Teuchos::RCP<Epetra_CrsMatrix> BlockDiagA11;
      
      //! the Schur-complement matrix
      Teuchos::RCP<Epetra_CrsMatrix> Chat;

	  //! Options to use a scaling in the Chat solve
	  //! and remove zeros on its diagonals
	  bool scaleChat, fixSingularChat, printSingularChat, fixSingularA11;
	  double fixChatTol;
	  
	  //! scaling vector for Chat
      Teuchos::RCP<Epetra_Vector> scalingChat;
	  
      //! some temporary vectors
      Teuchos::RCP<Epetra_Vector> x1,x2,b1,b2,xtmp1,xtmp2,btmp1,btmp2;
      
#ifdef HAVE_ZOLTAN      
      //! object for repartitioning Chat
      Teuchos::RCP<Repart> RepartChat;
      Teuchos::RCP<Epetra_Vector> RepartChatVec1,RepartChatVec2,RepartChatVec3,RepartChatVec4;
#endif      
      //! Krylov solvers for the diagonal blocks (may be null,
      //! in which case only the preconditioner is applied once)
      Teuchos::RCP<AztecOO> A11Solver, ChatSolver;
      
      //! num iter for solves
      int nitA11,nitChat;

      //! num iter for solves
      double tolA11,tolChat;
      
      //! preconditioner for MAuv
      Teuchos::RCP<Epetra_Operator> A11Precond;
      
      //! preconditioner for the Schur complement Chat
      Teuchos::RCP<Epetra_Operator> ChatPrecond;
      
      //! fix pressure in these two local points (-1 if n/a)
      int fixp1, fixp2;
      
      //! fix pressure in fixp1/2 to this value
      double valp;
      
    private:


      //! after constructing the (approximate) Schur-Complement C^,
      //! we do something magical to cope with its singularity...
      void AdjustChat(Teuchos::RCP<Epetra_CrsMatrix> Chat);
            
      //! extract and invert 2x2 block diagonal from a CRS matrix
      
      //! bdiag should be allocated before and will be Filled() after the call.
      void ExtractInverseBlockDiagonal(const Epetra_CrsMatrix& A, Epetra_CrsMatrix& bdiag);
      
      //! apply SI or SL preconditioner inverse to a pre-split vector
      int ApplyInverse(Epetra_Vector& b1, Epetra_Vector& b2,
                        Epetra_Vector& x1, Epetra_Vector& x2, 
                        bool trans) const;
        
  };    //end of class SppSimplePrec

////////////////////////////////////////////////////////////////////////

}// namespace TRIOS

#endif //SADDLEPOINT_H
