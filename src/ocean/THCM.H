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

//----------------------------------------------------------------------
// THCM is a singleton, there can be only one instance at a time.
// As base class Singleton is templated we must instantiate it
//----------------------------------------------------------------------

template <typename THCM>
Teuchos::RCP<THCM> Singleton<THCM>::instance;

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
//!  is A=df/du.
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
    //! Indices to denote different flux fields
    enum  FluxID
    {
        _Sal = 0,  // Total salinity flux
        _QSOA,     // ocean-atmos salinity flux
        _QSOS,     // ocean-seaice salinity flux
        _Temp,     // Total temperature flux
        _QSW,      // shortwave radiative flux
        _QSH,      // sensible heat flux
        _QLH,      // latent heat flux
        _QTOS,     // ocean-seaice flux
        _MSI       // seaice mask
    };

    //! Struct to hold pointers to derivative fields
    struct Derivatives
    {
        //! Derivative of the temperature equation with respect to a
        //! sea ice mask.
        Teuchos::RCP<Epetra_Vector> dFTdM;

        //! Derivative of the salinity equation with respect to sea
        //!ice heat flux.
        Teuchos::RCP<Epetra_Vector> dFSdQ;

        //! Derivative of the salinity equation with respect to a
        //! sea ice mask.
        Teuchos::RCP<Epetra_Vector> dFSdM;

        //! Derivative of the salinity equation with respect to the
        //! sea ice integral correction gamma.
        Teuchos::RCP<Epetra_Vector> dFSdG;
    };

    const Teuchos::ParameterList& getParameters();
    void setParameters(Teuchos::ParameterList&);

    static Teuchos::ParameterList getDefaultInitParameters();
    static Teuchos::ParameterList getDefaultParameters();

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
      we keep sigma set to 0. Use DiagB() to access the B matrix.

      If maskTest is true we compute the Jacobian just for testing the
      landmask. This means that we temporarily switch off restoring
      conditions and the integral condition.
    */
    bool evaluate (const Epetra_Vector& solnVector,
                   Teuchos::RCP<Epetra_Vector> rhsVector,
                   bool computeJac = false,
                   bool maskTest = false);

    //! only recompute the diagonal matrix B

    /*! the matrix B is used by THCM to 'switch off' some equations.
      It is diagonal and the diagonal entries can be accessed by the
      function DiagB(). We use the opposite sign as compared to THCM.
      The entries of B are (see assemble.f::fillcolB and
      usrc.F90::matrix): Ro for u,v 0 for w,p 1 for T,S
    */
    void evaluateB();

    //! Compute the forcing matrix used for the stochasic forcing
    bool computeForcing();

    //! ------------------- I-EMIC couplings --------------------
    //!Flags enabling the coupling with the I-EMIC, separated into a
    //!heat and a salinity flux contribution.
    //!
    //! 0: standalone thcm
    //! 1: accepting external forcing
    int coupledT_, coupledS_, coupledM_;

    int getCoupledT() { return coupledT_; }
    int getCoupledS() { return coupledS_; }

    //! Set atmosphere temperature in the ocean model
    void setAtmosphereT(Teuchos::RCP<Epetra_Vector> const &atmosT);

    //! Set atmosphere humidity field in the ocean model
    void setAtmosphereQ(Teuchos::RCP<Epetra_Vector> const &atmosQ);

    //! Set atmosphere humidity field in the ocean model
    void setAtmosphereA(Teuchos::RCP<Epetra_Vector> const &atmosA);

    //! Set atmosphere precipitation field in the ocean model
    void setAtmosphereP(Teuchos::RCP<Epetra_Vector> const &atmosP);

    //! Set sea ice mask
    void setSeaIceQ(Teuchos::RCP<Epetra_Vector> const &seaiceQ);

    //! Set sea ice surface temperature
    void setSeaIceM(Teuchos::RCP<Epetra_Vector> const &seaiceM);

    //! Set sea ice integral correction
    void setSeaIceG(Teuchos::RCP<Epetra_Vector> const &seaiceG);

    //! Set emip in the ocean model
    void setEmip(Teuchos::RCP<Epetra_Vector> const &emip, char mode = 'D');

    //! Set tatm in the ocean model
    void setTatm(Teuchos::RCP<Epetra_Vector> const &emip);

    //! Obtain non-overlapping emip
    Teuchos::RCP<Epetra_Vector> getEmip(char mode = 'D');

    //! Obtain shortwave radiation influence field
    Teuchos::RCP<Epetra_Vector> getSunO();

    //! Get a number of fluxes
    std::vector<Teuchos::RCP<Epetra_Vector> > getFluxes();

    Derivatives getDerivatives();

    //! Obtain local atmosphere temperature vector
    Teuchos::RCP<Epetra_Vector> getLocalAtmosT();

    //! Obtain local atmosphere humidity vector
    Teuchos::RCP<Epetra_Vector> getLocalAtmosQ();

    //! Obtain atmosphere humidity vector
    Teuchos::RCP<Epetra_Vector> getAtmosQ();

    //! Obtain local atmosphere precipitation vector
    Teuchos::RCP<Epetra_Vector> getLocalAtmosP();

    //! Obtain local ocean evaporation vector
    Teuchos::RCP<Epetra_Vector> getLocalOceanE();

    //! Obtain ocean evaporation vector
    Teuchos::RCP<Epetra_Vector> getOceanE();

    //! Return a shared pointer to the ocean's gathered landmask
    //! for use in the atmosphere
    std::shared_ptr<std::vector<int> > getLandMask();

    //! Return distributed mask read from maskName
    //! A fix vector can be supplied containing additional grid points that
    //! should be set to land.
    Teuchos::RCP<Epetra_IntVector>
    getLandMask(std::string const &maskName,
                Teuchos::RCP<Epetra_Vector> fix = Teuchos::null);

    //! Set local (distributed) landmask in THCM
    void setLandMask(Teuchos::RCP<Epetra_IntVector> landmask, bool init = true);

    //! Set global landmask in THCM
    void setLandMask(std::shared_ptr<std::vector<int> > landmask);

    //! ------------------------------------------------------------------
    //! Returns initial guess (Global/Solve form)
    Teuchos::RCP<Epetra_Vector> getSolution();

    //! returns the Jacobian matrix (Global/Solve form)
    Teuchos::RCP<Epetra_CrsMatrix> getJacobian();

    //! returns the Forcing matrix (Global/Solve form)
    Teuchos::RCP<Epetra_CrsMatrix> getForcing();

    //! get null space
    Teuchos::RCP<const Epetra_MultiVector> getNullSpace();

    //! get row scaling vector
    Teuchos::RCP<Epetra_Vector> getRowScaling() {return rowScaling_;}

    //! get column scaling vector
    Teuchos::RCP<Epetra_Vector> getColScaling() {return colScaling_;}

    //! Set a bifurcation parameter in the application physics
    bool setParameter(std::string label, double value);

    //! Get a bifurcation parameter from the application physics
    bool getParameter(std::string label, double& value);

    //! Write all THCM params in fort.7
    bool writeParams();

    //! returns the domain decomposition object
    Teuchos::RCP<TRIOS::Domain> GetDomain() {return domain_;}

    //! get the diagonal matrix B
    Teuchos::RCP<Epetra_Vector> DiagB(){return diagB_;}

    //! return the communicator object
    Teuchos::RCP<Epetra_Comm> GetComm(){return comm_;}

    //! get the SRES parameter (non-restoring salt condition)
    bool getSRES() const {return sres_;}

    //! get the TRES parameter (non-restoring temp condition)
    bool getTRES() const {return tres_;}

    //! get the its parameter (non-restoring salt condition)
    bool getITS() const {return its_;}

    //! get the ite parameter (non-restoring temp condition)
    bool getITE() const {return ite_;}

    //! get the salinity flux correction
    double getSCorr();

    //! ite_/its_=0 T/S forcing from data (Levitus)
    int ite_,its_;

    //! use internal temperature and salinity forcing
    bool internal_forcing_;

    //! get the global index of the row with the integral condition (if sres_==0)
    int getRowIntCon() const {return rowintcon_;}

    Teuchos::RCP<Epetra_Vector> getIntCondCoeff();

    //! Set the THCM flag 'vmix_fix' to 0 or 1. set the vmix_fix flag
    //! (required for controlling mixing and convective adjustment
    //! continuation/time-stepping). Note: vmix_fix doesn't have to be
    //! set if you have vmix_flag=1 in mix_imp.f (recommended).
    void fixMixing(int value);

    //! convert parameter name to integer (i.e. "Combined Forcing" => 19)
    int par2int(std::string const &label);

    // convert parameter index to std::string (i.e. 19 => "Combined Forcing")
    static std::string const int2par(int ind);

    //! Under non-restoring conditions: get an integral from vec and
    //! use that to correct the integral condition
    void setIntCondCorrection(Teuchos::RCP<Epetra_Vector> vec);

    //! Let THCM perform integral checks
    void integralChecks(Teuchos::RCP<Epetra_Vector> state,
                        double &salt_advection,
                        double &salt_diffusion);

private:

    //! object for domain-decomposition:
    Teuchos::RCP<TRIOS::Domain> domain_;

    //!\name maps and import/exprt objects for distributed data structures
    //! (obtained from the domain-object)
    //!@{
    //! overlapping map for 'THCM objects':
    Teuchos::RCP<Epetra_Map> assemblyMap_;
    Teuchos::RCP<Epetra_Map> assemblySurfaceMap_;
    Teuchos::RCP<Epetra_Map> assemblyVolumeMap_;

    //! Surface assembly to standardmap importer
    Teuchos::RCP<Epetra_Import> as2std_surf_;
    Teuchos::RCP<Epetra_Import> as2std_vol_;

    //! non-overlapping map for Trilinos objects:
    Teuchos::RCP<Epetra_Map> standardMap_;
    Teuchos::RCP<Epetra_Map> standardSurfaceMap_;
    Teuchos::RCP<Epetra_Map> standardVolumeMap_;

    //! non-overlapping load-balanced map for Trilinos objects
    //! (may contain non-rectangular subdomains)
    Teuchos::RCP<Epetra_Map> solveMap_;
    //!@}

    //! used only to define vector format (i.e. map), I think
    Teuchos::RCP<Epetra_Vector> initialSolution_;

    //! used to import current approximation to THCM:
    Teuchos::RCP<Epetra_Vector> localSol_;

    //! used to import atmosphere temperature into THCM
    Teuchos::RCP<Epetra_Vector> localAtmosT_;

    //! used to import atmosphere humidity field into THCM
    Teuchos::RCP<Epetra_Vector> localAtmosQ_;

    //! used to import albedo field into THCM
    Teuchos::RCP<Epetra_Vector> localAtmosA_;

    //! used to import atmosphere precipitation field into THCM
    Teuchos::RCP<Epetra_Vector> localAtmosP_;

    //! used to import sea ice heat flux into THCM
    Teuchos::RCP<Epetra_Vector> localSeaiceQ_;

    //! used to import sea ice mask into THCM
    Teuchos::RCP<Epetra_Vector> localSeaiceM_;

    //! used to import sea ice integral correction into THCM
    Teuchos::RCP<Epetra_Vector> localSeaiceG_;

    //! used to extract evaporation field from THCM
    Teuchos::RCP<Epetra_Vector> localOceanE_;

    //! used to import emip field into THCM
    Teuchos::RCP<Epetra_Vector> localEmip_;

    //! a surface overlapping vector without meaning
    Teuchos::RCP<Epetra_Vector> localSurfTmp_;

    //! used to import tatm field into THCM
    Teuchos::RCP<Epetra_Vector> localTatm_;

    //! used to export computed RHS vector from THCM:
    Teuchos::RCP<Epetra_Vector> localRhs_;

    //! the diagonal matrix B stored as a vector
    Teuchos::RCP<Epetra_Vector> localDiagB_, diagB_;

    //! Jacobian in globally assembled and load-balanced
    Teuchos::RCP<Epetra_CrsMatrix> jac_;

    //! Jacobian based on standard subdomains
    Teuchos::RCP<Epetra_CrsMatrix> localJac_, testJac_;

    //! Forcing in globally assembled and load-balanced
    Teuchos::RCP<Epetra_CrsMatrix> frc_;

    //! Forcing based on standard subdomains
    Teuchos::RCP<Epetra_CrsMatrix> localFrc_;

    //! type of scaling to be done
    std::string scalingType_;

    //! vectors storing the THCM scaling
    Teuchos::RCP<Epetra_Vector> rowScaling_, colScaling_,localRowScaling_,localColScaling_;

    //! (MPI) communicator
    Teuchos::RCP<Epetra_Comm> comm_;

    //! nullspace (p-vectors)
    Teuchos::RCP<Epetra_MultiVector> nullSpace_;

    //! global shared parameter list
    Teuchos::ParameterList paramList_;

    //! \name the matrix in THCM (CSR-format)
    /*! Memory for the Jacobian on the subdomain is allocated by the
      C++ code. In THCM, pointers are set to these locations so that
      the Fortran code fills these arrays directly. Afterwards the
      matrix is copied to the 'jac_' member, ignoring 'ghost' rows.
    */
    //!@{
    //! row pointer
    int* begA_;
    //! column indices
    int* jcoA_;
    //! values
    double* coA_;
    //! values of the diagonal matrix B
    double* coB_;
    //! row pointer
    int* begF_;
    //! column indices
    int* jcoF_;
    //! values
    double* coF_;
    //!@}

    //! global grid dimensions
    int n_,m_,l_;

    //! periodic domain in x-direction?
    bool periodic_;

    //! compute salinity integral
    bool compSalInt_;

    //! mixing?
    int vmix_;

    //! sres_=0: non-restoring salt forcing => integral condition in A and f
    int sres_;

    //! While sres_ = 0 (flux forcing in the local model), do not
    //! create integral condition.
    bool localSres_;

    //! tres_=0: non-restoring temperature forcing => integral condition in A and f
    int tres_;

    //! integral condition can have a negative or a positive sign
    int intSign_;

    //! which row is replaced by integral condition (global index of last row)
    int rowintcon_;

    //! correction for the integral condition based on the salinity
    //! integral of the initial state
    double intCorrection_;

    //! correction for the salinity flux
    double scorr_;

    //! vector with coefficients for integral condition (if sres_=0)
    Teuchos::RCP<Epetra_Vector> intcondCoeff_;

    //! sum of integration coefficients (total volume)
    double totalVolume_;

    //! pressure points where equation is replaced by P=0 (-1 means none)
    int rowPfix1_, rowPfix2_;

    //! asks THCM to recompute scaling vectors
    void RecomputeScaling(void);

    //! distribute land array after global initialization
    Teuchos::RCP<Epetra_IntVector> distributeLandMask(Teuchos::RCP<Epetra_IntVector> landm_glob);

    //! implement integral condition for S in Jacobian and B-matrix
    void intcond_S(Epetra_CrsMatrix& A, Epetra_Vector& B);

    //! flag to switch Dirichlet values P=0 on/off
    bool fixPressurePoints_;

    //! implement Dirichlet values P=0 in cells rowPfix1_/2 (if >=0)
    void fixPressurePoints(Epetra_CrsMatrix& A, Epetra_Vector& B);

    //! this subroutine defines the maximal matrix graph (pattern of nonzeros in the jacobian
    //! if convective adjustment occurs in all cells).
    Teuchos::RCP<Epetra_CrsGraph> CreateMaximalGraph(bool useSRES = true);


    //! private function used by CreateMaximalGraph()
    void insert_graph_entry(int* indices, int& pos,
                            int i, int j, int k, int var,
                            int N, int M, int L) const;
};

#endif
