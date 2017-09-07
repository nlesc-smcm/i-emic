#ifndef ATMOSPHEREPAR_H
#define ATMOSPHEREPAR_H

//! This class takes care of the parallelization of the atmosphere model.
//! For the domain decomposition we make use of the TRIOS_Domain stuff.
//! This als ensures that, in a coupled case, we can move parallel objects
//! around with maps that make sense.

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_CrsGraph.h>
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>

#include "Atmosphere.H"
#include "TRIOS_Domain.H"
#include "GlobalDefinitions.H"

class Ocean;

class AtmospherePar
{
public:
    using ParameterList = Teuchos::RCP<Teuchos::ParameterList>;
    using VectorPtr = Teuchos::RCP<Epetra_Vector>;

private:
    //! serial atmosphere class, working on local subdomains
    std::shared_ptr<Atmosphere> atmos_;

    //! parameterlist
    ParameterList params_;

    //! communicator object
    Teuchos::RCP<Epetra_Comm> comm_;

    //! grid size
    int n_,m_,l_, dim_;

    //! periodic bdc
    bool periodic_;

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

    //! parallel atmosphere state (non-overlapping)
    Teuchos::RCP<Epetra_Vector> state_;

    //! Index map for atmosphere T values
    Teuchos::RCP<Epetra_BlockMap> tIndexMap_;

    //! Index map for atmosphere Q values
    Teuchos::RCP<Epetra_BlockMap> qIndexMap_;

    //! parallel atmosphere temperature vector (non-overlapping)
    Teuchos::RCP<Epetra_Vector> atmosT_;

    //! parallel atmosphere humidity vector (non-overlapping)
    Teuchos::RCP<Epetra_Vector> atmosQ_;

    //! parallel atmosphere state (overlapping)
    Teuchos::RCP<Epetra_Vector> localState_;

    //! parallel atmosphere rhs (non-overlapping)
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! parallel atmosphere rhs (overlapping)
    Teuchos::RCP<Epetra_Vector> localRHS_;

    //! parallel sst (non-overlapping)
    Teuchos::RCP<Epetra_Vector> sst_;

    //! parallel sst (overlapping)
    Teuchos::RCP<Epetra_Vector> localSST_;

    //! parallel atmosphere sol (non-overlapping)
    Teuchos::RCP<Epetra_Vector> sol_;

    //! parallel atmosphere sol (overlapping)
    Teuchos::RCP<Epetra_Vector> localSol_;

    //! parallel evaporation field (non-overlapping)
    Teuchos::RCP<Epetra_Vector> E_;

    //! parallel evaporation field (overlapping)
    Teuchos::RCP<Epetra_Vector> localE_;

    //! parallel precipitation field (non-overlapping)
    Teuchos::RCP<Epetra_Vector> P_;

    //! parallel precipitation field (overlapping)
    Teuchos::RCP<Epetra_Vector> localP_;

    //! Atmosphere temperature importer
    Teuchos::RCP<Epetra_Import> atmosTimporter_;

    //! Atmosphere humidity importer
    Teuchos::RCP<Epetra_Import> atmosQimporter_;

    //! State hash
    size_t stateHash_;
    
    //! input file (HDF5)
    std::string const inputFile_;

    //! output file (HDF5)
    std::string const outputFile_;

    //! input/output flags
    bool loadState_, saveState_, saveEP_, saveEveryStep_;

    //! flag allowing to disable integral condition
    bool useIntCondQ_;
    
    //! row in rhs and matrix for integral condition q
    int rowIntCon_;

    //! coefficients for integral condition
    Teuchos::RCP<Epetra_Vector> intcondCoeff_;

    //! coefficients for precipitation integral
    Teuchos::RCP<Epetra_Vector> precipIntCo_;

    //! total area for E,P
    double totalArea_;
    
    //! preconditioning initialization flag
    bool precInitialized_;

    //! preconditioning computation flag
    bool recomputePrec_;

    //! ifpack preconditioner object
    Teuchos::RCP<Ifpack_Preconditioner> precPtr_;

    //! Jacobian matrix
    Teuchos::RCP<Epetra_CrsMatrix> jac_;

    //! matrix graph, necessary to build Jacobian matrix
    Teuchos::RCP<Epetra_CrsGraph> matrixGraph_;

    //! global surface mask
    std::shared_ptr<std::vector<int> > surfmask_;

public:
    //! Constructor based on parameterlist
    AtmospherePar(Teuchos::RCP<Epetra_Comm> comm, ParameterList params);

    ~AtmospherePar() {INFO("AtmospherePar destructor");}

    //! Setup integral coefficients
    void setupIntCoeff();

    void computeRHS();

    void computeJacobian();

    //! Compute E and P fields
    //! If stateHash is not updated, nothing happens here
    void computeEP();

    void applyMatrix(Epetra_MultiVector const &in,
                     Epetra_MultiVector &out);

    void applyPrecon(Epetra_MultiVector &in,
                     Epetra_MultiVector &out);

    //! initialize preconditioner
    void initializePrec();

    void solve(Teuchos::RCP<Epetra_MultiVector> const &b);

    void setState(Teuchos::RCP<Epetra_Vector> input) { state_ = input; }

    //! Receive ocean data
    void synchronize(std::shared_ptr<Ocean> ocean);

    void setOceanTemperature(Teuchos::RCP<Epetra_Vector> sst);

    Teuchos::RCP<Epetra_Vector> getLocalOceanT() {return localSST_; }

    void setLandMask(Utils::MaskStruct const &mask);

    //! set continuation parameter
    void setPar(double value) {atmos_->setPar(value); }

    //! set parameter
    void setPar(std::string const &parName, double value) {atmos_->setPar(parName, value); }

    //! Set continuation parameter name
    void setParName(std::string const &parName) { atmos_->setParName(parName); }

    //! Get continuation parameter
    double getPar() { return atmos_->getPar(); }

    //! Get continuation parameter
    double getPar(std::string const &parName) { return atmos_->getPar(parName); }

    //! Get derivatives of E and P wrt To and q
    void getEPconstants(double &qdim, double &nuq,
                        double &eta, double &dqso);
    
    //! Get continuation parameter name
    std::string const getParName() { return atmos_->getParName(); }

    void preProcess();
    void postProcess();

    void printState(std::string const &fname);

    //! Let the atmos object use idealized values,
    //! mainly for testing purposes
    void idealized();

    //! Obtain problem size
    int dim() { return dim_; }

    //! Obtain degrees of freedom
    int dof() { return dof_; }

    Teuchos::RCP<Epetra_Vector> getVector(char mode, Teuchos::RCP<Epetra_Vector> vec);

    Teuchos::RCP<Epetra_Vector> getState(char mode = 'C') { return getVector(mode, state_); }
    Teuchos::RCP<Epetra_Vector> getSolution(char mode = 'C') { return getVector(mode, sol_); }
    
    Teuchos::RCP<Epetra_Vector> getRHS(char mode = 'C')   { return getVector(mode, rhs_); }

    Teuchos::RCP<Epetra_Vector> getE(char mode = 'C')   { return getVector(mode, E_); }
    Teuchos::RCP<Epetra_Vector> getP(char mode = 'C')   { return getVector(mode, P_); }

    Teuchos::RCP<Epetra_Vector> getPrecipIntCo(char mode = 'C')
        { return getVector(mode, precipIntCo_); }

    Teuchos::RCP<Epetra_CrsMatrix> getJacobian() { return jac_; }

    //! Return pointer to domain object
    Teuchos::RCP<TRIOS::Domain> getDomain() { return domain_; }

    //! Get temperature at interface 
    Teuchos::RCP<Epetra_Vector> interfaceT();

    //! Get humidity at interface 
    Teuchos::RCP<Epetra_Vector> interfaceQ();

    //! Get precipitation at interface 
    Teuchos::RCP<Epetra_Vector> interfaceP();
    
    //! Get integral condition coefficients
    Teuchos::RCP<Epetra_Vector> getIntCoeff() { return intcondCoeff_; } 

    Teuchos::RCP<Epetra_BlockMap> getStandardSurfaceMap() { return standardSurfaceMap_; }
    Teuchos::RCP<Epetra_BlockMap> getAssemblySurfaceMap() { return assemblySurfaceMap_; }

    //! Get ocean dependency
    std::shared_ptr<Utils::CRSMat> getBlock(std::shared_ptr<Ocean> ocean);

    //! Get row at the interface, XX is 1-based, i and j 0-based
    int interface_row(int i, int j, int XX)
        { return FIND_ROW_ATMOS0(ATMOS_NUN_, n_, m_, l_, i, j, l_-1, ATMOS_TT_); }

private:

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

    //! HDF5-based load function for the state and parameters
    int loadStateFromFile(std::string const &name);
    
    //! HDF5-based save function for the state and parameters
    int saveStateToFile(std::string const &name);

    //! Copy HDF5
    void copyFiles();

    //! Distribute parallel state in serial model
    void distributeState();

};

#endif
