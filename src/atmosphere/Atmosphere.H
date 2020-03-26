//!---------------------------------------------------------------------
/*

 This class takes care of the parallelization of the local atmosphere
 model AtmosLocal. For the domain decomposition we make use of the
 TRIOS_Domain stuff. This also ensures that, in a coupled case, we can
 move parallel objects around with maps that make sense.

 The physics are dealt with in AtmosLocal. Mos importantly in
 AtmosLocal::forcing() and AtmosLocal::computeRHS().
 
*/
//!--------------------------------------------------------------------

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_CrsGraph.h>
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>

#include "Model.H"
#include "AtmosLocal.H"
#include "TRIOS_Domain.H"
#include "GlobalDefinitions.H"

class Ocean;
class SeaIce;

class Atmosphere : public Model
{
public:    

    using CommPars      = AtmosLocal::CommPars;

protected:
    //! serial atmosphere class, working on local subdomains
    std::shared_ptr<AtmosLocal> atmos_;

    //! parameterlist
    ParameterList params_;

    //! grid size
    int n_,m_,l_, dim_;

    //! periodic bdc
    bool periodic_;

    //! auxiliary unknowns resulting in a bordered system
    int aux_;

    //! degrees of freedom
    int dof_;

    //! grid dimensions
    double xmin_,xmax_,ymin_,ymax_;

    //! TRIOS domain object for parallelization strategy
    Teuchos::RCP<TRIOS::Domain> domain_;

    //! standard map, without ghost nodes (non-overlapping).
    Teuchos::RCP<Epetra_Map> standardMap_;

    //! assembly map, with ghost nodes (overlapping).
    Teuchos::RCP<Epetra_Map> assemblyMap_;

    //! depth-averaged, single unknown for ocean surface temperature
    //! standard map, without ghost nodes (non-overlapping).
    Teuchos::RCP<Epetra_Map> standardSurfaceMap_;

    //! depth-averaged, single unknown for ocean surface temperature
    //! assembly map, with ghost nodes (overlapping).
    Teuchos::RCP<Epetra_Map> assemblySurfaceMap_;

    //! Surface assembly to standardmap importer
    Teuchos::RCP<Epetra_Import> as2std_surf_;

    //! parallel atmosphere state (overlapping)
    Teuchos::RCP<Epetra_Vector> localState_;

    //! parallel atmosphere rhs (non-overlapping)
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! parallel atmosphere rhs (overlapping)
    Teuchos::RCP<Epetra_Vector> localRHS_;

    // //! parallel diagonal B (mass matrix, non-overlapping)
    Teuchos::RCP<Epetra_Vector> diagB_;

    // //! parallel diagonal B (mass matrix, overlapping)
    Teuchos::RCP<Epetra_Vector> localDiagB_;

    //! parallel land surface temperature lst (non-overlapping)
    Teuchos::RCP<Epetra_Vector> lst_;

    //! parallel lst (overlapping)
    Teuchos::RCP<Epetra_Vector> localLST_;

    //! parallel sea surface temperature sst (non-overlapping)
    Teuchos::RCP<Epetra_Vector> sst_;

    //! parallel sst (overlapping)
    Teuchos::RCP<Epetra_Vector> localSST_;

    //! parallel atmosphere sol (non-overlapping)
    Teuchos::RCP<Epetra_Vector> sol_;

    //! parallel atmosphere sol (overlapping)
    Teuchos::RCP<Epetra_Vector> localSol_;

    //! parallel sea ice surface temperature sit (non-overlapping)
    Teuchos::RCP<Epetra_Vector> sit_;

    //! parallel sit (overlapping)
    Teuchos::RCP<Epetra_Vector> localSIT_;

    //! parallel sea ice mask (non-overlapping)
    Teuchos::RCP<Epetra_Vector> Msi_;

    //! parallel sea ice mask (overlapping)
    Teuchos::RCP<Epetra_Vector> localMSI_;

    //! parallel evaporation field (non-overlapping)
    Teuchos::RCP<Epetra_Vector> E_;

    //! parallel evaporation field (overlapping)
    Teuchos::RCP<Epetra_Vector> localE_;

    //! parallel precipitation field (non-overlapping)
    Teuchos::RCP<Epetra_Vector> P_;

    //! parallel precipitation distribution (non-overlapping)
    Teuchos::RCP<Epetra_Vector> Pdist_;

    //! parallel precipitation field (overlapping)
    Teuchos::RCP<Epetra_Vector> localP_;

    //! State component maps
    std::map<int, Teuchos::RCP<Epetra_Map> > Maps_;
    
    //! State component importers
    std::map<int, Teuchos::RCP<Epetra_Import> > Imps_;

    //! flag to disable integral condition
    bool useIntCondQ_;
    
    //! row in rhs and matrix for integral condition q
    int rowIntCon_;

    //! coefficients for integral condition
    Teuchos::RCP<Epetra_Vector> intcondCoeff_;

    //! coefficients for integral condition gathered on every proc
    Teuchos::RCP<Epetra_MultiVector> intcondGlob_;

    //! coefficients for precipitation integral
    Teuchos::RCP<Epetra_Vector> pIntCoeff_;

    //! total area for E,P
    double totalArea_;
    
    //! use idealized precipitation
    bool useFixedPrecip_;

    //! preconditioning initialization flag
    bool precInitialized_;

    //! preconditioning computation flag
    bool recomputePrec_;

    // //! mass matrix computation flag
    bool recompMassMat_;

    //! ifpack preconditioner object
    Teuchos::RCP<Ifpack_Preconditioner> precPtr_;

    //! Jacobian matrix
    Teuchos::RCP<Epetra_CrsMatrix> jac_;

    //! matrix graph, necessary to build Jacobian matrix
    Teuchos::RCP<Epetra_CrsGraph> matrixGraph_;

    //! global surface mask
    std::shared_ptr<std::vector<int> > surfmask_;

public:
    //! Constructor based on communicator and parameterlist
    Atmosphere(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);

    ~Atmosphere() {INFO("Atmosphere destructor");}
    
    std::string name() const { return "atmos"; }
    virtual int modelIdent() const { return 1; }
    
    //! Setup integral coefficients
    void setupIntCoeff();

    void computeRHS();

    void computeJacobian();

    //! fill epetra vector with constant P
    void fillP();
    void applyMatrix(Epetra_MultiVector const &in,
                     Epetra_MultiVector &out);

    void applyMassMat(Epetra_MultiVector const &in,
                      Epetra_MultiVector &out);

    void applyPrecon(Epetra_MultiVector const &in,
                     Epetra_MultiVector &out);

    //! build mass matrix
    void computeMassMat();
    
    //! initialize pr    
    void buildPreconditioner();

    void solve(Teuchos::RCP<Epetra_MultiVector> const &b);

    void setState(Teuchos::RCP<Epetra_Vector> input) { state_ = input; }

    //! Receive ocean data
    void synchronize(std::shared_ptr<Ocean> ocean);

    //! Receive ocean data
    void synchronize(std::shared_ptr<SeaIce> seaice);

    //! Meaningless: dummy implementation
    void synchronize(std::shared_ptr<Atmosphere> atmos) {}

    //! Meaningless: dummy implementation
    void pressureProjection(Teuchos::RCP<Epetra_Vector>) {}

    void setOceanTemperature(Teuchos::RCP<Epetra_Vector> sst);
    
    void setSeaIceTemperature(Teuchos::RCP<Epetra_Vector> sst);
    
    void setSeaIceMask(Teuchos::RCP<Epetra_Vector> sst);

    Teuchos::RCP<Epetra_Vector> getLocalSST() {return localSST_; }
    Teuchos::RCP<Epetra_Vector> getLandTemperature();

    void setLandMask(Utils::MaskStruct const &mask);

    // Atmosphere does not have this functionality, return empty landmask
    Utils::MaskStruct getLandMask()
        {
            WARNING("Not implemented", __FILE__, __LINE__);
            return Utils::MaskStruct();
        }

    //! set parameter
    void setPar(std::string const &parName, double value) {atmos_->setPar(parName, value); }

    //! Get continuation parameter
    double getPar(std::string const &parName) { return atmos_->getPar(parName); }

    //! Get parameters that we want to communicate
    void getCommPars(CommPars &parStruct) const;

    void preProcess();
    void postProcess();

    //! Gather important continuation data to use in summary
    std::string writeData(bool describe = false) const; 

    void printState(std::string const &fname);

    //! Let the atmos object use idealized values,
    //! mainly for testing purposes
    void idealized(double defP = 1.0);

    //! Obtain problem size (mainly for testing)
    int dim() { return dim_; }

    //! Obtain degrees of freedom
    int dof() { return dof_; }

    //! Return number of continuation parameters
    int npar() { return atmos_->npar(); }

    //! Return parameter name, supply index
    std::string int2par(int ind) const { return atmos_->int2par(ind); }

    Teuchos::RCP<Epetra_Vector> getState(char mode = 'C')
        { return Utils::getVector(mode, state_); }
    Teuchos::RCP<Epetra_Vector> getSolution(char mode = 'C')
        { return Utils::getVector(mode, sol_); }
    Teuchos::RCP<Epetra_Vector> getRHS(char mode = 'C')
        { return Utils::getVector(mode, rhs_); }
    Teuchos::RCP<Epetra_Vector> getMassMat(char mode = 'C')
        { return Utils::getVector(mode, diagB_); }

    Teuchos::RCP<Epetra_CrsMatrix> getJacobian() { return jac_; }

    //! Return pointer to domain object
    Teuchos::RCP<TRIOS::Domain> getDomain() { return domain_; }

    //! Factorized interface routine
    Teuchos::RCP<Epetra_Vector> interface(int XX) const;
    
    //! Get temperature at interface 
    Teuchos::RCP<Epetra_Vector> interfaceT() const;

    //! Get humidity at interface 
    Teuchos::RCP<Epetra_Vector> interfaceQ() const;

    //! Get albedo at interface 
    Teuchos::RCP<Epetra_Vector> interfaceA() const;

    //! Get precipitation at interface 
    Teuchos::RCP<Epetra_Vector> interfaceP() const;

    //! Compute evaporation at interface
    Teuchos::RCP<Epetra_Vector> interfaceE();

    //! Get integral condition coefficients for standard map
    Teuchos::RCP<Epetra_Vector> getIntCoeff() { return intcondCoeff_; }

    //! Get precipitation integral coefficients for surface map (nun=1)
    Teuchos::RCP<Epetra_Vector> getPIntCoeff(char mode = 'C')
        { return Utils::getVector(mode, pIntCoeff_); }

    Teuchos::RCP<Epetra_BlockMap> getStandardSurfaceMap() { return standardSurfaceMap_; }
    Teuchos::RCP<Epetra_BlockMap> getAssemblySurfaceMap() { return assemblySurfaceMap_; }

    //! Get ocean dependency
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Ocean> ocean);

    //! Get sea ice dependency
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<SeaIce> seaice);

    //! Dummy implementation
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Atmosphere> atmos)
        { return std::shared_ptr<Utils::CRSMat>(); }

    //! Get row at the interface, XX is 1-based, i and j 0-based. The
    //! result is also 0-based.
    int interface_row(int i, int j, int XX);

    //! Get row containing the humidity integral condition
    int getRowIntCon() { return rowIntCon_; }

    //! Get spatial precipitation distribution
    Teuchos::RCP<Epetra_Vector> getPdist() { return Pdist_; }
    
private:

    //! set precipitation distribution
    void setPdist();

    //! Private member to compute P if it is not an unknown in the
    //! model.
    Teuchos::RCP<Epetra_Vector> getP(Teuchos::RCP<Epetra_Vector> P = Teuchos::null) const;

    //! Private member to compute P_.
    void setP();

    //! Private member to compute E
    Teuchos::RCP<Epetra_Vector> getE(char mode = 'C');

    //! Private member to compute E_.
    void setE();

    //! create dependency graph for the Jacobian matrix
    //! function borrowed from THCM,
    //! --> could be factorized within the project when bored
    void createMatrixGraph();

    //! insert a dependency into the graph
    //! function borrowed from THCM,
    //! --> could be factorized within the project when bored
    void insert_graph_entry(int* indices, int& pos,
                            int i, int j, int k, int xx,
                            int N, int M, int L) const;

    void additionalImports(EpetraExt::HDF5 &HDF5, std::string const &filename){}

    //! HDF5-based save function for other components than the state
    //! and parameters.
    void additionalExports(EpetraExt::HDF5 &HDF5, std::string const &filename);

    //! Assemble fluxes from local model
    std::vector<Teuchos::RCP<Epetra_Vector> > getFluxes();

    //! Distribute parallel state in serial model
    void distributeState();

};

#endif
