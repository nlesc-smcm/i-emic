/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef TRIOS_BLOCKPRECONDITIONER_H
#define TRIOS_BLOCKPRECONDITIONER_H

#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Ifpack_Preconditioner.h"

// typedef'd Teuchos pointers

class Epetra_MultiVector;
class AztecOO;

namespace TRIOS {

//! identifiers for submatrices (used for accessing 'matrix arrays')
    enum MatrixID {_Auv,_Dw,_ATS,_BTSuv,_BTSw,_BwTS,_Gw,_Guv,_Duv};
// number of matrices in the arrays (see below)
#define _NUMSUBM 9


// forward declaration
    class Domain;
    class ApMatrix;
    class SaddlepointMatrix;
    class SppSimplePrec;
    class Repart;


    //! Block-ILU preconditioner for Trilinos-THCM

    /*! This class implements the Block-ILU preconditioner
      known as 'BGS' in THCM.

      For details about this preconditioner, please refer to
      "A tailored solver for bifurcation analysis of ocean-climate models"
      by Arie de Niet and Fred Wubs, J. Comp. Physics 2006

      \verbatim
      Notes on the notation used in this class
      ========================================
      \endverbatim

      The naming of variables follows the fortran implementation in THCM 6.0,
      at least to some extent. The factors L and U of the preconditioner are
      given by (eqn. 26/27 in the paper with the 'saddlepoint' problem replaced
      by the depth-averaged one):

      \verbatim


      | Ap    0     0    0   0      |      | I 0 0 0      Ap\BwTS            |
      | Guv  MAuv  MGuv  0   0      |      | 0 I 0 0 -Auv\Guv Ap\BwTS        |
 L =  |  0   MDuv   0    0   0      |, U = | 0 0 I 0         0               |
      |  0   Duv1   0    Aw  0      |      | 0 0 0 I Aw\Duv1\Auv\Guv Ap\BwTS |
      |  0   BTSvu  0  BTSw SchurTS |      | 0 0 0 0         I               |


      with SchurTS = ATS+(Buv-BwAw\Duv1) Auv\Guv Ap\BwTS.

      \endverbatim

      Generally the preconditioner attempts to keep the same indexing which is used
      in the global matrix, i.e. a vector x_uv will be based on a map containing
      entries of the form [0 1 6 7 12 13 ...] etc, and likewise for matrices.
      On the other hand, 'depth-averaged' quantities are numbered contiguously,
      i.e. after depth averaging, the same vector will be based on a map
      [0 1 2 ... n*m*2]. Depth-averaged quantities are denoted by 'bar' (mapPbar)
      or a 'z', i.e. xzp (which was the notation used in the ortran code).

      A term encountered frequently in this class is 'dummy' (W or P cell). It re-
      fers to grid cells were the vertical velocity and/or the pressure have no
      physical meaning because they are land cells (P and W) or the top layer of
      ocean cells (only dummy W). Consequently there are always n*m dummy W's and
      the number of dummy W's is always n*m + (number of dummy P's). Dummy points
      are 'ignored' by the preconditioner, which is implemented by removing them
      from all Epetra_Map objects. Maps without dummies generally have a '1'
      attached to their name. Some pressure 'P' maps are reduced so that the dummy
      W cells are removed as well. These are denoted by a '^' (i.e. mapPhat).

      Algorithm
      =========

      the full block-ILU factorization is not feasible because of the Schur-comple-
      ment operator SchurTS. We therefore use some permutation of the Jacobian which
      has few entries in the upper triangular block and then use Block-Gauss-Seidel.
      Three different permutations are implemented (see Arie's thesis, p. 104).
      The permutation used is selected by the option "Permutation" (1, 2, or 3):

           | Ap   0   0  BTS |       | Kuvp 0   0 Guv |       | Aw  0   0  Duv  |
      M1 = | Guv Kuvp 0   0  |, M2 = | Duv  Aw  0  0  |, M3 = | Bw ATS  0  Buv  |
           |  0  Duv  Aw  0  |       | Buv  Bw ATS 0  |       | 0  BTS Ap   0   |
           | Buv  0   Bw ATS |       |  0   0  BTS Ap |       | 0   0  Guv Kuvp |

    */
    class BlockPreconditioner: public Ifpack_Preconditioner
    {
    public:

        /*!
          \brief Constructor

          Performs the most general setup steps which do not need the Jacobian to be
          computed yet (not even the pattern). The domain object must be present, and
          its 'Decomp' function must have been called. The parameters are unfortunately
          not yet documented as they are subject to frequent change.
        */
        BlockPreconditioner(Teuchos::RCP<Epetra_CrsMatrix> jacobian,
                            Teuchos::RCP<Domain> domain_,
                            Teuchos::ParameterList& List);

        //! Destructor
        virtual ~BlockPreconditioner();

        //! Set transpose (not implemented => returns -1).
        int SetUseTranspose(bool UseTranspose) {return -1;}

        //! Apply block-ILU preconditioning operator (not implemented)
        int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

        //! Apply preconditioner operator inverse

        /*! The input and output vectors should be based on the standard
          'Solve' map which can be obtained from the domain object (or from
          the Jacobian, which should be based on the same map).
          Multiple RHS are NOT supported, so X and Y should actually be
          Epetra_Vector objects.
        */
        int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

        //! Computing infinity norm (not implemented)
        double NormInf() const;

        //! Label
        const char* Label() const {return label_.c_str();}

        //! Transposed? - returns false
        bool UseTranspose() const {return false;}

        //! Have norm-inf? returns false
        bool HasNormInf() const {return false;}

        /*!
         * \brief Returns a reference to the Epetra_Comm communicator associated
         * with this operator.
         */
        const Epetra_Comm& Comm() const {return *comm;}

        /*!
         * \brief Returns the Epetra_Map object associated with the domain of
         * this operator.
         */
        const Epetra_Map& OperatorDomainMap() const {return *domainMap;}

        /*!
         * \brief Returns the Epetra_Map object associated with the range of
         * this operator.
         */
        const Epetra_Map& OperatorRangeMap() const {return *rangeMap;}

        //! Compute preconditioner \f$M\f$.
        bool computePreconditioner(const Epetra_Vector& x,
                                   Epetra_Operator& Prec,
                                   Teuchos::ParameterList* p = NULL);


        //! convert Teuchos::ParameterList to aztec options/params (static helper function)
        static void ExtractAztecOptions(Teuchos::ParameterList& list, int* options, double* params);

        void setOutputStreams(Teuchos::RCP<std::ostream> innerOut, Teuchos::RCP<std::ostream> outerOut,
                              Teuchos::RCP<std::ostream> innerErr=Teuchos::null,
                              Teuchos::RCP<std::ostream> outerErr=Teuchos::null)
            {
                InnerStream=innerOut;
                OuterStream=outerOut;
                if (innerErr==Teuchos::null)
                    InnerErrorStream = InnerStream;
                else
                    InnerErrorStream = innerErr;
                if (outerErr==Teuchos::null)
                    OuterErrorStream = OuterStream;
                else
                    OuterErrorStream = outerErr;
            }

        //! \name Ifpack_Preconditinoer interface
        //!@{

        //! this constructor is compulsory for Ifpack_AdditiveSchwarz.
        //! It is not recommended, though: it gets its domain from the
        //! THCM instance and sets only default parameters.
        BlockPreconditioner(Epetra_RowMatrix* MatrixPtr);
        int SetParameters(Teuchos::ParameterList &List);
        int Initialize();
        bool IsInitialized() const;
        int Compute();
        bool IsComputed() const;
        double Condest() const;
        double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                       const int MaxIters = 1550,
                       const double Tol = 1e-9,
                       Epetra_RowMatrix* Matrix = 0);

        const Epetra_RowMatrix& Matrix() const;
        int NumInitialize() const;
        int NumCompute() const;
        int NumApplyInverse() const;
        double InitializeTime() const;
        double ComputeTime() const;
        double ApplyInverseTime() const;
        double InitializeFlops() const;
        double ComputeFlops() const;
        double ApplyInverseFlops() const;
        std::ostream& Print(std::ostream& os) const;
#ifdef TESTING
        //=============================================================================
        // Get some temporary access for debugging
        //=============================================================================
        void Test();
#endif
        //!@}

    protected:

        //! verbosity (0-10)
        int verbose;

        //! can be "ILU" or "Gauss-Seidel", obtained from
        //! parameter entry "Ocean Preconditioner"->"Block Precond"->"Scheme"
        std::string scheme;

        //! Permutation used (see class description)
        int permutation;

        //! relaxation parameter for GS-type preconditioner
        double DampingFactor;

        //! subtract checkerboard mode from pressure
        //! ("Pressure Correction" option, defaults to true)
        bool DoPresCorr;

        //! Reference to parameter list for controlling linear system solve.
        Teuchos::ParameterList lsParams;

        //! Label for this operator.
        std::string label_;

        //! the matrix of the linear system we want to precondition
        Teuchos::RCP<Epetra_CrsMatrix> jacobian;

        //! domain decomposition (constructed by OceanProblemInterface)
        Teuchos::RCP<Domain> domain;

        //! communicator
        Teuchos::RCP<Epetra_Comm> comm;

        //! domain and range maps
        Teuchos::RCP<Epetra_Map> domainMap, rangeMap;

        //! number of 'dummy' w points, i.e. grid cells
        //! where no w-equation is solved (certain land cells)
        int dumw;

        //! number of 'dummy' p points, i.e. grid cells
        //! where no p-equation is solved (certain land cells)
        int dump;

        //! boolean arrays that indicate which cells are dummy P and/or W:
        bool *is_dummyP, *is_dummyW;

        //! parts of the 'nun 6' map: uv, w, p, TS
        Teuchos::RCP<Epetra_Map> mapUV, mapW, mapP, mapTS;

        //! replicated versions for column maps (these are used when splitting the matrix)
        Teuchos::RCP<Epetra_Map> colmapUV, colmapW, colmapP, colmapTS;

        //! W and P maps without 'dummies'. Dummy P's are
        //! located in land cells, Dummy W's are land cells
        //! and the top layer of ocean cells.
        Teuchos::RCP<Epetra_Map> mapW1,mapP1,colmapW1,colmapP1;

        //! since there are per definition more dummy W points than
        //! dummy P's, this map contains only P points where there
        //! is a non dummy W as well (i.e. the top ocean layer is
        //! removed). This map is used to make P/W matrices square.
        Teuchos::RCP<Epetra_Map> mapPhat;

        //! maps for the depth-averaged pressure and horizontal velocities.
        //! These maps use contiguous indexing.
        Teuchos::RCP<Epetra_Map> mapPbar, mapUVbar;

        //! importers for the submaps
        Teuchos::RCP<Epetra_Import> importUV,importW,importP,importTS,importP1,importW1,importPhat;


        //! \name submatrices extracted from the Jacobian

        //! name    description/properties                 Fortran name/enum id
        //! ----    ----------------------                 ------------     ----
        //! A_{uv}   convection-diffusion/coriolis               Auv
        //! G_w      p-grad in w-eqn, upper tri                  Gw
        //! A_w      lower triangular                            Aw
        //! A_{ST}   conv-diff eq. for T and S                   ATS
        //! B_{uv}   block for the TS equation                   BTSuv
        //! B_w      block for the TS equation                   BTSw
        //! ~B~_ST                                               BwTS
        //! G_{uv}   discrete gradient                           Guv
        //! D_{uv}   discrete divergence                         Duv
        //!

        //!@{


        //! arary of matrices
        Teuchos::RCP<Epetra_CrsMatrix> SubMatrix[_NUMSUBM];

        //! for each matrix a pointer to the row map it needs
        Teuchos::RCP<Epetra_Map> SubMatrixRowMap[_NUMSUBM];

        //! for each matrix a pointer to the column map it needs
        //! (may be Teuchos::null)
        Teuchos::RCP<Epetra_Map> SubMatrixColMap[_NUMSUBM];

        //! for each matrix a pointer to the range map it needs
        Teuchos::RCP<Epetra_Map> SubMatrixRangeMap[_NUMSUBM];

        //! for each matrix a pointer to the domain map it needs
        Teuchos::RCP<Epetra_Map> SubMatrixDomainMap[_NUMSUBM];

        //! importers for the row map <-> total row map
        Teuchos::RCP<Epetra_Import> Importer[_NUMSUBM];

        //! labels to identify the matrices
        std::string SubMatrixLabel[_NUMSUBM];

        //@}

        /*! \name Diagonal Blocks for linear system solves

          whereas the matrices in the 'SubMatrix' array have
          a full colmap to allow convenient importing, the
          operators actually used in linear system solves must
          have rowmap=colmap
        */
        //@{


        //! Duv1 is the 'square part' of Duv, it maps from UV to W1 instead of P1
        Teuchos::RCP<Epetra_CrsMatrix> Duv1;

        //! Similar to Duv1, Aw is the 'square part' of Dw, so Aw:W1->W1
        Teuchos::RCP<Epetra_CrsMatrix> Aw;

        //! extra copies of diagonal blocks Auv and ATS with replaced colmaps.
        Teuchos::RCP<Epetra_CrsMatrix> Auv,ATS;

        //! orthogonal transform to make ATS easier to solve and transformed system
        Teuchos::RCP<Epetra_CrsMatrix> QTS,Arhomu;

        //! linear map for Arhomu (to allow using direct solvers via Amesos)
        Teuchos::RCP<Epetra_Map> Arhomu_linearmap;

        //! extra scaling vectors for ATS/Arhomu
        Teuchos::RCP<Epetra_Vector> ATS_leftscaling, ATS_rightscaling;

        //! \name operator representations for block matrices
        //!@{

        //! diagonal block matrix Ap, the square part of Gw (W1->Phat)
        //! it has it's own operator class so that we can wrap the
        //! triangular solve.
        Teuchos::RCP<ApMatrix> Ap;

        //! saddlepoint matrix [Auv,Guv; Duv,0].
        //! The pressure is always depth-averaged,
        //! depth-averaging the velocity field is optional
        Teuchos::RCP<SaddlepointMatrix> Spp;

        //!@}


        //@}

        //! depth-averaging operators for pressure and hor. velocity:
        Teuchos::RCP<Epetra_CrsMatrix> Mzp1, Mzp2;

        //! if true, all systems are solved with a 0 initial guess
        bool zero_init;

        //! singular vectors of the pressure.

        /*! these vectors represent 'checkerboard modes' of the pressure field,
          they are constructed such that Gw*svp = Guv*svp=0 and are used to project
          the spurious modes out of the computed pressure.
        */
        Teuchos::RCP<Epetra_Vector> svp1, svp2;

        //!\name Solvers and prexconditioners for subsystems
        //!@{

        /*! \note since we currently use standard GMRES in the outer iteration
          the number of steps done by inner solvers must be fixed. To avoid this
          problem, one should use flexible (F-) GMRES instead.
        */

        //! Krylov solver for the Saddlepoint Problem
        Teuchos::RCP<AztecOO> SppSolver;

        //! max number of iterations for Spp
        int nitSpp;

        //! required tolerance for Spp system
        double tolSpp;

        //! Krylov solver for Auv
        Teuchos::RCP<AztecOO> AuvSolver;

        //! max number of iterations for Auv
        int nitAuv;

        //! required tolerance for Spp system
        double tolAuv;

        //! Krylov solver for SchurTS (or ATS)
        Teuchos::RCP<AztecOO> ATSSolver;

        //! repartitioner for ATS
        Teuchos::RCP<Repart> RepartATS;

        //! max number of iterations for Spp
        int nitATS;

        //! required tolerance for Spp system
        double tolATS;

        //! preconditioner for Auv
        Teuchos::RCP<Epetra_Operator> AuvPrecond;

        //! preconditioner for SchurTS
        Teuchos::RCP<Epetra_Operator> ATSPrecond;

        //! preconditioner for Spp
        Teuchos::RCP<Epetra_Operator> SppPrecond;

        //@}

        //! file stream for outer iteration: Aztec output
        //! will be redirected to this stream _after_ the
        //! preconditioner has been applied
        Teuchos::RCP<std::ostream> OuterStream;

        //! file stream for outer iteration: Aztec error messages
        //! will be redirected to this stream _after_ the
        //! preconditioner has been applied
        Teuchos::RCP<std::ostream> OuterErrorStream;

        //! file streams to be used for inner AztecOO solves
        Teuchos::RCP<std::ostream> InnerStream,InnerErrorStream;

    protected:

        //! number of unknowns per grid cell (THCM)
        static const int dof_=6;

        //! extract all submatrices from the Jacobian
        void extract_submatrices(const Epetra_CrsMatrix&);

        //! builds solvers and blockmatrices
        void build_preconditioner(void);

        //! find dummy rows in a subset of rows of the matrix A.

        /*! Given a matrix A and a (sub-)map M, this function
          sets the flag array is_dummy to true if the corresponding
          local node is a dummy row, i.e. if it has only one entry
          and that entry is 1.0 on the diagonal of the matrix.
          Memory for the is_dummy array is allocated in this function.

          returns the number of dummies found.
        */
        int detect_dummies(const Epetra_CrsMatrix& A, const Epetra_Map& M, bool *is_dummy) const;

        //! builds the depth-averaging operators Mp as the matrix of singular vectors
        //! of a given gradient operator like Gw or Dw'.
        Teuchos::RCP<Epetra_CrsMatrix> build_singular_matrix(Teuchos::RCP<Epetra_CrsMatrix> Gw);

        //! construct singular vectors of P (two 'checkerboard modes')
        void build_svp();

        //! construct QTS*QTS=I and Arhomu=QTS*ATS*QTS so that Arhomu is easier
        //! to solve than ATS if convective adjustment is switched on.
        void setup_rhomu();

        //! in an extracted global matrix row, find the entry with
        //! indices[pos]=col.
        //! returns true if the entry was found, false otherwise.
        //! entry's position in the indices array is placed in pos.
        bool find_entry(int col, int *indices, int numentries, int& pos);

#ifdef TESTING
        //! check if pressure singular vectors are correct
        bool test_svp();
#endif
        //! if true, setup() has to be called. This is needed because
        //! the Jacobian may be unavailable when the constructor is called.
        bool needs_setup;

        //! for ifpack interface, not sure this is implemented correctly
        bool IsComputed_;

        //! this is called inside the constructor
        void Setup1();

        //! constructs the submaps/matrices (needs to be called only once, but not
        //! before a 'real' Jacobian is available)
        void Setup2();


        //! lower triangular solve with the factor L of the approximate Jacobian
        //! (Solve Lx=b for x). We have three versions of this function for the
        //! three permutations (see class description).
        void SolveLower1(const Epetra_Vector& buv, const Epetra_Vector& bw,
                         const Epetra_Vector& bp,  const Epetra_Vector& bTS,
                         Epetra_Vector& xuv, Epetra_Vector& xw,
                         Epetra_Vector& xp, Epetra_Vector& xTS) const;

        //! lower triangular solve with the factor L of the approximate Jacobian
        //! (Solve Lx=b for x). We have three versions of this function for the
        //! three permutations (see class description).
        void SolveLower2(const Epetra_Vector& buv, const Epetra_Vector& bw,
                         const Epetra_Vector& bp,  const Epetra_Vector& bTS,
                         Epetra_Vector& xuv, Epetra_Vector& xw,
                         Epetra_Vector& xp, Epetra_Vector& xTS) const;

        //! lower triangular solve with the factor L of the approximate Jacobian
        //! (Solve Lx=b for x). We have three versions of this function for the
        //! three permutations (see class description).
        void SolveLower3(const Epetra_Vector& buv, const Epetra_Vector& bw,
                         const Epetra_Vector& bp,  const Epetra_Vector& bTS,
                         Epetra_Vector& xuv, Epetra_Vector& xw,
                         Epetra_Vector& xp, Epetra_Vector& xTS) const;

        //! upper triangular solve with the factor U of the approximate Jacoibian
        //! (Solve Ux=b for x)
        void SolveUpper(const Epetra_Vector& buv, const Epetra_Vector& bw,
                        const Epetra_Vector& bp,  const Epetra_Vector& bTS,
                        Epetra_Vector& xuv, Epetra_Vector& xw,
                        Epetra_Vector& xp,  Epetra_Vector& xTS) const;

        //! solve linear system with ATS, satisfying integral condition
        //! for S if SRES==0.
        void SolveATS(Epetra_Vector& rhs, Epetra_Vector& sol,
                      double tol, int maxit) const;

        //! store Jacobian, rhs, start guess and all the preconditioner 'hardware'
        //! (i.e. depth-averaging operators etc) in an HDF5 file
        void dumpLinSys(const Epetra_Vector& x, const Epetra_Vector& b) const;



    }; // end of class BlockPreconditioner



//////////////////////////////////////////////////////////////////////////////////////////
//! \name Auxiliary operator classes for the Block-ILU preconditioner.

/*! Auxiliary operator classes (block systems, preconditioners etc.)                     //
  note: these are special purpose classes, they are not intended to be used outside    //
  class BlockPreconditioner.                                                                     //
*/
//!@{
//////////////////////////////////////////////////////////////////////////////////////////


//!
//! This is a matrix class that represents the matrix Ap. Ap can be
//! implemented either as a block system (as in Arie's original paper)
//! or simply as the 'square part' of Gw (as in the final thesis).
//! This class only supports application of the inverse of Ap,
//! which can be computed cheaply from a block-factorization or by
//! exploiting the triangular structure of Gw. We use the latter (and
//! simpler) implementation.
//!
//! *range and domain*
//!
//! Ap = Gw~:  P1 -> W1
//! (i.e. pressure and vert. vel. maps without dummies)
//!
//! \note internally the maps are adjusted to allow calls to 'ApplyInverse'
//!
    class ApMatrix
    {

    public:
//-----//


        //! constructor

        /*! apart from the matrix this constructor also craves some maps.
          '1' indicates that dummies should be removed from the standard maps;
          'Phat' are those p points which lie in cells where there is a non-dummy W.
        */
        ApMatrix(const Epetra_CrsMatrix &Gw,
                 const Epetra_CrsMatrix &Mp,
                 Teuchos::RCP<Epetra_Map> mapW1,
                 Teuchos::RCP<Epetra_Map> mapP1,
                 Teuchos::RCP<Epetra_Map> mapPhat,
                 Teuchos::RCP<Epetra_Comm> comm_,
                 char ApType_);

        //! destructor
        virtual ~ApMatrix();

        //! get range map (which is W1)
        const Epetra_Map& RangeMap() const {return *rangeMap;}

        //! get domain map (which is P1)
        const Epetra_Map& DomainMap() const {return *domainMap;}

        //! apply inverse operator x=Ap\b

        /*! Here b should be based on the 'W1' map,
          and X on the 'P1' map
        */
        int ApplyInverse (const Epetra_Vector &b, Epetra_Vector &x) const;


    protected:
//--------//

        //! range and domain map
        Teuchos::RCP<Epetra_Map> rangeMap, domainMap;

        //! modified Gw-matrix (maps are replaced)
        Teuchos::RCP<Epetra_CrsMatrix> Gw1, Gw2;

        //! modified nullspace of Gw
        Teuchos::RCP<Epetra_CrsMatrix> Mp1, Mp2;

        //! we need some importers
        Teuchos::RCP<Epetra_Import> importPbar, importPhat;

        //! to create temporary vectors etc we need some maps:
        Teuchos::RCP<Epetra_Map> mapW1, mapP1;

        char ApType;

    }; // end of class ApMatrix

// BILU is no longer supported, we keep the class below just in case
#if 0

//////////////////////////////////////////////////////////////////////

//! operator representing \f$ S_{TS}\f$

/*! This is a matrix class that represents our Schur-
  complement problem (SchurTS).
  It only supports multiplication by a vector
  and construction from precomputed blocks.
*/
    class SchurTSMatrix : public Epetra_Operator
    {

    public:
//-----//


        //! constructor
        SchurTSMatrix(Teuchos::RCP<Epetra_CrsMatrix> ATS_, Teuchos::RCP<Epetra_CrsMatrix> BTSuv_,
                      Teuchos::RCP<Epetra_CrsMatrix> BwTS_, Teuchos::RCP<Epetra_CrsMatrix> Aw_, Teuchos::RCP<Epetra_CrsMatrix> Duv1_,
                      Teuchos::RCP<Epetra_CrsMatrix> Guv_, Teuchos::RCP<Epetra_CrsMatrix> BTSw_,
                      Teuchos::RCP<AztecOO> SolverAuv, Teuchos::RCP<ApMatrix> Ap_,
                      Teuchos::RCP<Epetra_Map> mapTS_, Teuchos::RCP<Epetra_Comm> comm_);

        //! destructor
        virtual ~SchurTSMatrix();

//!\name  Atribute set methods
//!@{

        //! If set true, transpose of this operator will be applied (n/a)
        int   SetUseTranspose (bool UseTranspose){return -1;}

//!@}

//!\name Mathematical functions
//!@{

        //! apply operator Y=Op*X
        int   Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

        //! apply inverse operator (n/a)
        int   ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
            {
                (*info) << "ERROR: SchurTSMatrix::ApplyInverse is not available!\n";
                return -1;
            }

        //! Returns the infinity norm of the global matrix (n/a)
        double NormInf () const{return -1.0;}

//!@}

//!\name Atribute access functions


        //!Returns a character string describing the operator.
        const char *  Label () const {return label_.c_str();}

        //! Returns the current UseTranspose setting.
        bool  UseTranspose () const {return false;}

        //! Returns true if the this object can provide an approximate Inf-norm, false otherwise.
        bool  HasNormInf () const {return false;}

        //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
        const Epetra_Comm &   Comm () const {return *comm;}

        //! Returns the Epetra_Map object associated with the domain of this operator.
        const Epetra_Map &    OperatorDomainMap () const
            {
                return *domainMap;
            }

        //! Returns the Epetra_Map object associated with the range of this operator.
        const Epetra_Map &    OperatorRangeMap () const
            {
                return *rangeMap;
            }

        //\name set some properties of the inner Krylov solve with Auv
        //@{

        //! be verbose (or not)
        void SetPrintInnerSolve(bool what) {PrintInnerSolve=what;}

        //! number of iterations
        void SetAuvNumIter(int howmany) {nitAuv=howmany;}

        //! tolerance
        void SetAuvTol(double howmuch) {tolAuv=howmuch;}
        //@}

    protected:
//--------//


        //! label identifying this class
        std::string label_;

        //! communicator
        Teuchos::RCP<Epetra_Comm> comm;

        //! domain map (map for vectors x in Ax=y)
        Teuchos::RCP<Epetra_Map> domainMap;

        //! range map (map for vectors y in Ax=y)
        Teuchos::RCP<Epetra_Map> rangeMap;

        //! factors for this operator
        Teuchos::RCP<Epetra_CrsMatrix> ATS,BTSuv,BTSw,Duv1,Aw,Guv, BwTS;
        Teuchos::RCP<AztecOO> SolverAuv;

        //! factored matrix Ap
        Teuchos::RCP<ApMatrix> Ap;

        //! number of iterations for Auv system
        int nitAuv;

        //! required accuracy for Auv system
        double tolAuv;

        //! if true, the inner Krylov solve with Auv is printed
        bool PrintInnerSolve;


    };

#endif

}//namespace

#endif // OCEAN_PREC_H
