/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 *                                                                    *
 **********************************************************************/
/* -Messed up by Erik                                                 *
 **********************************************************************/

// for I-EMIC couplings
#include <math.h>
#include <iostream>
#include <sstream>
#include <memory>
#include <vector>
#include <algorithm>

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

#include "EpetraExt_MultiComm.h"
#include "EpetraExt_BlockVector.h"

// define global macros such as UU, _NUN_, INFO etc.
#include "THCMdefs.H"

// from TRIOS
#include "TRIOS_Domain.H"

// from trilinos_thcm
#include "THCM.H"

#ifdef DEBUGGING
#include "OceanGrid.H"
#endif

//============================================================================
extern "C" {

    void set_global_x(int*, double*);
    void set_global_xu(int*, double*);
    void set_global_y(int*, double*);
    void set_global_yv(int*, double*);
    void set_global_z(int*, double*);
    void set_global_zw(int*, double*);

    // usrc.F90
    _SUBROUTINE_(setparcs)(int* param, double* value);
    _SUBROUTINE_(getparcs)(int* param, double* value);
    _SUBROUTINE_(writeparams)();
    _SUBROUTINE_(rhs)(double* un, double* b);
    _SUBROUTINE_(setsres)(int* sres);
    _SUBROUTINE_(matrix)(double* un);
    _SUBROUTINE_(stochastic_forcing)();

    _SUBROUTINE_(init)(int* n, int* m, int* l, int* nmlglob,
                       double* xmin, double* xmax, double* ymin, double* ymax,
                       double* alphaT, double* alphaS,
                       int* ih, int* vmix, int* tap, int* rho_mixing,
                       int* coriolis_on,
                       int* periodic, int* landm,
                       double* taux, double* tauy, double* tatm, double* emip, double* spert);

    _SUBROUTINE_(set_landmask)(int* landm, int* periodic, int* reinit);

    _SUBROUTINE_(finalize)(void);

    // global.F90
    _MODULE_SUBROUTINE_(m_global,initialize)(int* N, int* M, int* L,
                                             double* Xmin, double* Xmax,
                                             double* Ymin, double* Ymax,
                                             double* hdim, double* qz,
                                             int* periodic, int* itopo, int* flat, int* rd_mask,
                                             int* TRES, int* SRES, int* iza, int* ite ,int* its, int* rd_spertm,
                                             int* coupled_T, int* coupled_S,
                                             int* forcing_type,
                                             const char *maskfile, const char *spertmaskfile,
                                             const char *windfile, const char *sstfile, const char *sssfile);

    _MODULE_SUBROUTINE_(m_global,finalize)(void);
    _MODULE_SUBROUTINE_(m_global,set_maskfile)(const char *maskfile);
    _MODULE_SUBROUTINE_(m_global,get_landm)(int* landm);
    _MODULE_SUBROUTINE_(m_global,get_current_landm)(int* landm);
    _MODULE_SUBROUTINE_(m_global,set_landm)(int* landm);
    _MODULE_SUBROUTINE_(m_global,get_windfield)(double* taux, double* tauy);
    _MODULE_SUBROUTINE_(m_global,get_temforcing)(double* tatm);
    _MODULE_SUBROUTINE_(m_global,get_salforcing)(double* emip);
    _MODULE_SUBROUTINE_(m_global,get_internal_temforcing)(double* temp);
    _MODULE_SUBROUTINE_(m_global,get_internal_salforcing)(double* salt);
    _MODULE_SUBROUTINE_(m_global,get_spert)(double* spert);

    _MODULE_SUBROUTINE_(m_usr,set_internal_forcing)(double*temp, double* salt);
    _MODULE_SUBROUTINE_(m_thcm_utils,get_landm)(int*);
    _MODULE_SUBROUTINE_(m_scaling,average_block)(double *db);
    _MODULE_SUBROUTINE_(m_scaling,compute)(double *db, double *rowscales, double* colscales);

    _SUBROUTINE_(fillcolb)(void);

    // CRS matrix allocation (module m_mat)
    _MODULE_SUBROUTINE_(m_mat,get_array_sizes)(int* nrows, int* nnz);
    _MODULE_SUBROUTINE_(m_mat,set_pointers)(int* nrows, int* nnz,
                                            int* begA_,int* jcoA_,double* coA_,
                                            double* coB_,
                                            int* begF_,int* jcoF_,double* coF_);

    // compute scaling factors for S-integral condition. Values is an n*m*l array
    _MODULE_SUBROUTINE_(m_thcm_utils,intcond_scaling)(double* values,int* indices,int* len);

    // compute weights for load balancing
    _MODULE_SUBROUTINE_(m_thcm_utils,loadbal_weights)(double* weights,double*,double*,double*);

    // sets the vmix_fix flag
    _MODULE_SUBROUTINE_(m_mix,set_vmix_fix)(int* vmix_fix);

    //---------------------- I-EMIC couplings--------------------------------------
    // Extensions created for communication within the I-EMIC
    //
    _MODULE_SUBROUTINE_(m_inserts, insert_atmosphere_t)(double *atmosT);
    _MODULE_SUBROUTINE_(m_inserts, insert_atmosphere_q)(double *atmosQ);
    _MODULE_SUBROUTINE_(m_inserts, insert_atmosphere_a)(double *atmosQ);
    _MODULE_SUBROUTINE_(m_inserts, insert_atmosphere_p)(double *atmosP);
    _MODULE_SUBROUTINE_(m_inserts, insert_seaice_q)(double *seaiceT);
    _MODULE_SUBROUTINE_(m_inserts, insert_seaice_m)(double *seaiceM);
    _MODULE_SUBROUTINE_(m_inserts, insert_seaice_g)(double *seaiceM);
    _MODULE_SUBROUTINE_(m_inserts, insert_emip)(double *emip);
    _MODULE_SUBROUTINE_(m_inserts, insert_adapted_emip)(double *emip);
    _MODULE_SUBROUTINE_(m_inserts, insert_emip_pert)(double *emip);
    _MODULE_SUBROUTINE_(m_inserts, insert_tatm)(double *emip);

    _MODULE_SUBROUTINE_(m_integrals, salt_advection)(double *un, double *check);
    _MODULE_SUBROUTINE_(m_integrals, salt_diffusion)(double *un, double *check);

    _MODULE_SUBROUTINE_(m_probe,  get_emip)(double *emip);
    _MODULE_SUBROUTINE_(m_probe,  get_suno)(double *suno);
    _MODULE_SUBROUTINE_(m_probe,  get_derivatives)(double *sol, double *dftdm,
                                                   double *dfsdq, double *dfsdm,
                                                   double *dfsdg);
    _MODULE_SUBROUTINE_(m_probe,  get_adapted_emip)(double *emip);
    _MODULE_SUBROUTINE_(m_probe,  get_emip_pert)(double *emip);
    _MODULE_SUBROUTINE_(m_probe,  get_salflux)(double *sol, double *salflux,
                                               double *scorr,
                                               double *qsoaflux, double *qsosflux);

    _MODULE_SUBROUTINE_(m_probe,  get_temflux)(double *sol, double *totflux,
                                               double *swflux, double *shflux,
                                               double *lhflux, double *siflux,
                                               double *simask);

    _MODULE_SUBROUTINE_(m_probe,  get_atmosphere_t)(double *atmosT);
    _MODULE_SUBROUTINE_(m_probe,  get_atmosphere_q)(double *atmosQ);
    _MODULE_SUBROUTINE_(m_probe,  get_atmosphere_p)(double *atmosP);
    _MODULE_SUBROUTINE_(m_probe,  compute_evap)(double *oceanE, double *x);

    _MODULE_SUBROUTINE_(m_global, get_land_temp)(double *land);
    //-----------------------------------------------------------------------------

    _SUBROUTINE_(write_levitus)(const char*);

}//extern

//=============================================================================
// constructor
THCM::THCM(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Comm> comm) :
    Singleton<THCM>(Teuchos::rcp(this, false)),
    comm_(comm),
    nullSpace_(Teuchos::null),
    paramList_("THCM Parameter List")
{
    DEBUG("### enter THCM::THCM ###");

    params.validateParametersAndSetDefaults(getDefaultInitParameters());
    paramList_.setParameters(params);

    std::string probdesc = paramList_.get<std::string>("Problem Description");

    n_ = paramList_.get<int>("Global Grid-Size n");
    m_ = paramList_.get<int>("Global Grid-Size m");
    l_ = paramList_.get<int>("Global Grid-Size l");

    std::stringstream ss;
    ss << "THCM (" << probdesc <<", "
       << n_ << "x" << m_ << "x" << l_ <<")";
    INFO(ss.str());
    this->SetLabel(ss.str().c_str());

    // default: north atlantic
    double xmin = paramList_.get<double>("Global Bound xmin") * PI_ / 180.0;
    double xmax = paramList_.get<double>("Global Bound xmax") * PI_ / 180.0;
    double ymin = paramList_.get<double>("Global Bound ymin") * PI_ / 180.0;
    double ymax = paramList_.get<double>("Global Bound ymax") * PI_ / 180.0;
    periodic_   = paramList_.get<bool>("Periodic");

    // sanity check
    double xdist = pow(cos(xmax)-cos(xmin), 2) + pow(sin(xmax)-sin(xmin), 2);
    if (xdist < 1e-2 && !periodic_)
    {
        WARNING("Periodic bdc disabled while \n"
                << " horizontal boundaries coincide. Distance: "
                << xdist, __FILE__, __LINE__);
    }

    double hdim  = paramList_.get<double>("Depth hdim");
    double qz    = paramList_.get<double>("Grid Stretching qz");
    int  itopo   = paramList_.get<int>("Topography");
    bool flat    = paramList_.get<bool>("Flat Bottom");
    compSalInt_  = paramList_.get<bool>("Compute salinity integral");

    bool rd_mask          = paramList_.get<bool>("Read Land Mask"); //== false in experiment0
    std::string mask_file = paramList_.get<std::string>("Land Mask");

    int ih             = paramList_.get<int>("Inhomogeneous Mixing");
    vmix_              = paramList_.get<int>("Mixing");
    bool rho_mixing    = paramList_.get<bool>("Rho Mixing");
    int tap            = paramList_.get<int>("Taper");
    double alphaT      = paramList_.get<double>("Linear EOS: alpha T");
    double alphaS      = paramList_.get<double>("Linear EOS: alpha S");
    tres_              = paramList_.get<int>("Restoring Temperature Profile");
    sres_              = paramList_.get<int>("Restoring Salinity Profile");
    localSres_         = paramList_.get<bool>("Local SRES Only");
    intSign_           = paramList_.get<int>("Salinity Integral Sign");
    ite_               = paramList_.get<int>("Levitus T");
    its_               = paramList_.get<int>("Levitus S");
    internal_forcing_  = paramList_.get<bool>("Levitus Internal T/S");
    coupledT_          = paramList_.get<int>("Coupled Temperature");
    coupledS_          = paramList_.get<int>("Coupled Salinity");
    coupledM_          = paramList_.get<int>("Coupled Sea Ice Mask");
    fixPressurePoints_ = paramList_.get<bool>("Fix Pressure Points");
    int coriolis_on    = paramList_.get<int>("Coriolis Force");
    int forcing_type   = paramList_.get<int>("Forcing Type");

    //------------------------------------------------------------------
    if ((coupledS_ == 1) && (sres_ == 1))
    {
        WARNING("Incompatible parameters: coupledS_ = " << coupledS_
                << " SRES = "
                << sres_ << " setting SRES = 0", __FILE__, __LINE__);
        sres_ = 0;
    }

    bool rd_spertm          = paramList_.get<bool>("Read Salinity Perturbation Mask");
    std::string spertm_file = paramList_.get<std::string>("Salinity Perturbation Mask");

    //------------------------------------------------------------------
    if (std::abs(intSign_) != 1)
    {
        ERROR("Invalid integral sign!", __FILE__, __LINE__);
    }

    // wind forcing method (0: data (trenberth), 1: zonally averaged, 2: idealized)
    int iza                = paramList_.get<int>("Wind Forcing Type");
    std::string windf_file = paramList_.get<std::string>("Wind Forcing Data");
    std::string temf_file  = paramList_.get<std::string>("Temperature Forcing Data");
    std::string salf_file  = paramList_.get<std::string>("Salinity Forcing Data");

    //------------------------------------------------------------------
    int dof = _NUN_; // number of unknowns, defined in THCMdefs.H

    // construct an object to decompose the domain:
    domain_ = Teuchos::rcp(new TRIOS::Domain(n_, m_, l_, dof, xmin, xmax, ymin, ymax,
                                            periodic_, hdim, qz, comm_));

    // perform a 2D decomposition of domain into rectangular boxes
    domain_->Decomp2D();

    // get a map object representing the subdomain (with ghost-nodes/overlap).
    // This map defines the nodes local to the THCM subdomain
    assemblyMap_ = domain_->GetAssemblyMap();

    // get a map object representing the subdomain (without ghost-nodes/overlap).
    // this is an intermediate representation between the assembly and the solve
    // phases.
    standardMap_ = domain_->GetStandardMap();

    // initialize THCM (allocate memory etc.)
    // for a subdomain including ghost-nodes:

    // the domain object knows the geometry of the subdomain:
    double xminloc = domain_->XminLoc();
    double xmaxloc = domain_->XmaxLoc();
    double yminloc = domain_->YminLoc();
    double ymaxloc = domain_->YmaxLoc();

    // and the number of grid points contained in it:
    int nloc = domain_->LocalN();
    int mloc = domain_->LocalM();
    int lloc = domain_->LocalL();

    // global settings for THCM
    int nglob_ = n_;
    int mglob_ = m_;
    int lglob_ = l_;

    // fortran routine that puts global info in module m_global
    // and allocates some arrays there if dimensions are non-zero
    // this has to be called _before_ init because init already uses
    // some of it (by calling 'forcing')

    // for portability reasons we rather pass integers to fortran than bools
    // (this was an issue on Huygens)
    int irho_mixing = (rho_mixing) ? 1 : 0;
    int iperiodic   = (periodic_ ) ? 1 : 0;
    int iflat       = (flat      ) ? 1 : 0;
    int ird_mask    = (rd_mask   ) ? 1 : 0;
    int ird_spertm  = (rd_spertm ) ? 1 : 0;

    INFO("THCM init: m_global::initialize...");
    INFO("    Mixing: vmix = " << vmix_);

    // In fortran object code this corresponds to the function
    //  __m_global_MOD_initialize
    F90NAME(m_global, initialize)(&nglob_, &mglob_, &lglob_,
                                  &xmin, &xmax, &ymin, &ymax, &hdim, &qz,
                                  &iperiodic, &itopo, &iflat, &ird_mask,
                                  &tres_, &sres_, &iza, &ite_, &its_, &ird_spertm,
                                  &coupledT_, &coupledS_,
                                  &forcing_type,
                                  mask_file.c_str(), spertm_file.c_str(),
                                  windf_file.c_str(), temf_file.c_str(), salf_file.c_str());

    INFO("THCM init: m_global::initialize... done");

    {
        int size;
        const TRIOS::Grid& grid = domain_->GetGlobalGrid();

        size = grid.x_.size();
        set_global_x(&size, grid.x_.get());
        size = grid.y_.size();
        set_global_y(&size, grid.y_.get());
        size = grid.z_.size();
        set_global_z(&size, grid.z_.get());
        size = grid.xu_.size();
        set_global_xu(&size, grid.xu_.get());
        size = grid.yv_.size();
        set_global_yv(&size, grid.yv_.get());
        size = grid.zw_.size();
        set_global_zw(&size, grid.zw_.get());
    }

    if (localSres_) // from here on we ignore the integral condition
        sres_ = 1;

    // read topography data and convert it to a global land mask
    DEBUG("Initialize land mask...");

    int I0 = 0; int I1 = n_+1;
    int J0 = 0; int J1 = m_+1;
    int K0 = 0; int K1 = l_+1;

    int i0=0, i1=-1, j0=0, j1=-1,k0=0,k1=-1;
    // global (gathered) map, all inds are on root proc, the ranges are
    if (comm->MyPID()==0)
    {
        i1 = I1; j1 = J1; k1=K1;
    }

// the next part of this file is devoted to taking arrays from m_global that
// have been read from files and putting them into m_usr, where they will be
// distributed and have two layers of overlap. This is done for many different
// arrays in exactly the same way, so it would call for some abstraction layer,
// but we currently do it separately for all arrays (TODO: make this general).
//
// The arrays are:
//
// landm: land mask, integer 3D
// taux,tauy: wind field, double 2D
// tatm: atmosphere temperature, double 2D
// emip: surface salinity forcing, double 2D
// temp: internal temperature forcing, double 3D
// salt: internal salinity forcing, double 3D

    DEBUG("Create gathered land map");
    Teuchos::RCP<Epetra_Map> landmap_glb =
        Utils::CreateMap(i0,i1,j0,j1,k0,k1,I0,I1,J0,J1,K0,K1,*comm);

    // sequential landm array on proc 0
    Teuchos::RCP<Epetra_IntVector> landm_glb =
        Teuchos::rcp(new Epetra_IntVector(*landmap_glb));

    int *landm;
    if (comm->MyPID()==0)
    {
        CHECK_ZERO(landm_glb->ExtractView(&landm));
        // make THCM fill the global landm array and put it into our C pointer location
        DEBUG("call m_global::get_landm");
        F90NAME(m_global,get_landm)(landm);
    }

    Teuchos::RCP<Epetra_IntVector> landm_loc = distributeLandMask(landm_glb);

    // import local landm-part to THCM
    CHECK_ZERO(landm_loc->ExtractView(&landm));

    // in the main part of THCM (except m_global) we set periodic
    // boundary conditions to .false. _unless_ we are running a
    // periodic problem on a single CPU in the x-direction:
    Teuchos::RCP<Epetra_Comm> xComm = domain_->GetProcRow(0);

    int perio = (periodic_ && xComm->NumProc() == 1);

    int nmlglob = n_*m_*l_;

    //--------------------------------------------------------------------------
    // read wind, temperature and salinity forcing and distribute it among
    // processors
    Teuchos::RCP<Epetra_Map> wind_map_loc   = domain_->CreateAssemblyMap(1,true);

    // Single-unknown surface maps
    standardSurfaceMap_ = domain_->CreateStandardMap(1,true);
    assemblySurfaceMap_ = domain_->CreateAssemblyMap(1,true);

    // Surface assembly/standard import strategy
    as2std_surf_ =
        Teuchos::rcp(new Epetra_Import(*assemblySurfaceMap_, *standardSurfaceMap_));

    // Single-unknown volume maps
    standardVolumeMap_ = domain_->CreateStandardMap(1,false);
    assemblyVolumeMap_ = domain_->CreateAssemblyMap(1,false);

    // Volume assembly/standard import strategy
    as2std_vol_ =
        Teuchos::rcp(new Epetra_Import(*assemblyVolumeMap_, *standardVolumeMap_));


    Teuchos::RCP<Epetra_Map> lev_map_loc    = wind_map_loc;
    Teuchos::RCP<Epetra_Map> intlev_map_loc = domain_->CreateAssemblyMap(1,false);
    Teuchos::RCP<Epetra_Vector> taux_loc    = Teuchos::rcp(new Epetra_Vector(*wind_map_loc));
    Teuchos::RCP<Epetra_Vector> tauy_loc    = Teuchos::rcp(new Epetra_Vector(*wind_map_loc));
    Teuchos::RCP<Epetra_Vector> tatm_loc    = Teuchos::rcp(new Epetra_Vector(*lev_map_loc));
    Teuchos::RCP<Epetra_Vector> emip_loc    = Teuchos::rcp(new Epetra_Vector(*lev_map_loc));
    Teuchos::RCP<Epetra_Vector> temp_loc    = Teuchos::rcp(new Epetra_Vector(*intlev_map_loc));
    Teuchos::RCP<Epetra_Vector> salt_loc    = Teuchos::rcp(new Epetra_Vector(*intlev_map_loc));
    Teuchos::RCP<Epetra_Vector> spert_loc   = Teuchos::rcp(new Epetra_Vector(*lev_map_loc));

    double* taux;
    double* tauy;
    double* tatm;
    double* emip;
    double* temp;
    double* salt;
    double* spert;

    CHECK_ZERO(taux_loc->ExtractView(&taux));
    CHECK_ZERO(tauy_loc->ExtractView(&tauy));
    CHECK_ZERO(tatm_loc->ExtractView(&tatm));
    CHECK_ZERO(emip_loc->ExtractView(&emip));
    CHECK_ZERO(temp_loc->ExtractView(&temp));
    CHECK_ZERO(salt_loc->ExtractView(&salt));
    CHECK_ZERO(spert_loc->ExtractView(&spert));

    DEBUG("Initialize Wind field...");

    Teuchos::RCP<Epetra_Map> wind_map_dist = domain_->CreateStandardMap(1,true);
    Teuchos::RCP<Epetra_Map> wind_map_root = Utils::Gather(*wind_map_dist,0);

    Teuchos::RCP<Epetra_Vector> taux_glob =
        Teuchos::rcp(new Epetra_Vector(*wind_map_root));
    Teuchos::RCP<Epetra_Vector> tauy_glob =
        Teuchos::rcp(new Epetra_Vector(*wind_map_root));

    double *taux_g, *tauy_g;
    CHECK_ZERO(taux_glob->ExtractView(&taux_g));
    CHECK_ZERO(tauy_glob->ExtractView(&tauy_g));

    if (comm->MyPID()==0)
    {
        std::cout << " obtaining windfield" << std::endl;
        F90NAME(m_global,get_windfield)(taux_g,tauy_g);
    }

    // distribute wind fields
    Teuchos::RCP<Epetra_MultiVector> taux_dist =
        Utils::Scatter(*taux_glob,*wind_map_dist);
    Teuchos::RCP<Epetra_MultiVector> tauy_dist =
        Utils::Scatter(*tauy_glob,*wind_map_dist);

    // import overlap
    Teuchos::RCP<Epetra_Import> wind_loc2dist =
        Teuchos::rcp(new Epetra_Import(*wind_map_loc,*wind_map_dist));
    CHECK_ZERO(taux_loc->Import(*taux_dist,*wind_loc2dist,Insert));
    CHECK_ZERO(tauy_loc->Import(*tauy_dist,*wind_loc2dist,Insert));

    double tauxmax, tauymax;
    double tauxmin, tauymin;
    CHECK_ZERO(taux_dist->MaxValue(&tauxmax));
    CHECK_ZERO(tauy_dist->MaxValue(&tauymax));
    CHECK_ZERO(taux_dist->MinValue(&tauxmin));
    CHECK_ZERO(tauy_dist->MinValue(&tauymin));

    INFO("Zonal wind forcing from data ranges between: ["
         << tauxmin << ".." << tauxmax << "]");
    INFO("Meridional wind forcing from data ranges between: ["
         << tauymin << ".." << tauymax << "]");

    DEBUG("Initialize Temperature and Salinity forcing...");

    Teuchos::RCP<Epetra_Map> lev_map_dist  = domain_->CreateStandardMap(1,true);

    Teuchos::RCP<Epetra_Map> lev_map_root  = Utils::Gather(*lev_map_dist,0);

    Teuchos::RCP<Epetra_Map> intlev_map_dist  = domain_->CreateStandardMap(1,false);
    Teuchos::RCP<Epetra_Map> intlev_map_root  = Utils::Gather(*intlev_map_dist,0);

    Teuchos::RCP<Epetra_Vector> tatm_glob     =
        Teuchos::rcp(new Epetra_Vector(*lev_map_root));
    Teuchos::RCP<Epetra_Vector> emip_glob     =
        Teuchos::rcp(new Epetra_Vector(*lev_map_root));
    Teuchos::RCP<Epetra_Vector> temp_glob     =
        Teuchos::rcp(new Epetra_Vector(*intlev_map_root));
    Teuchos::RCP<Epetra_Vector> salt_glob     =
        Teuchos::rcp(new Epetra_Vector(*intlev_map_root));
    Teuchos::RCP<Epetra_Vector> spert_glob    =
        Teuchos::rcp(new Epetra_Vector(*lev_map_root));

    double *tatm_g, *emip_g, *spert_g, *temp_g, *salt_g;
    CHECK_ZERO(tatm_glob->ExtractView(&tatm_g));
    CHECK_ZERO(emip_glob->ExtractView(&emip_g));
    CHECK_ZERO(temp_glob->ExtractView(&temp_g));
    CHECK_ZERO(salt_glob->ExtractView(&salt_g));
    CHECK_ZERO(spert_glob->ExtractView(&spert_g));

    if (comm->MyPID() == 0)
    {
        F90NAME(m_global, get_temforcing)(tatm_g);
        F90NAME(m_global, get_salforcing)(emip_g);
        if (internal_forcing_)
        {
            F90NAME(m_global,get_internal_temforcing)(temp_g);
            F90NAME(m_global,get_internal_salforcing)(salt_g);
        }
        else
        {
            temp_glob->PutScalar(0.0);
            salt_glob->PutScalar(0.0);
        }
        F90NAME(m_global,get_spert)(spert_g);
    }

    // distribute levitus fields
    Teuchos::RCP<Epetra_MultiVector> tatm_dist     =
        Utils::Scatter(*tatm_glob, *lev_map_dist);
    Teuchos::RCP<Epetra_MultiVector> emip_dist     =
        Utils::Scatter(*emip_glob, *lev_map_dist);
    Teuchos::RCP<Epetra_MultiVector> temp_dist     =
        Utils::Scatter(*temp_glob, *intlev_map_dist);
    Teuchos::RCP<Epetra_MultiVector> salt_dist     =
        Utils::Scatter(*salt_glob, *intlev_map_dist);
    Teuchos::RCP<Epetra_MultiVector> spert_dist    =
        Utils::Scatter(*spert_glob, *lev_map_dist);

    // import overlap
    Teuchos::RCP<Epetra_Import> lev_loc2dist =
        Teuchos::rcp(new Epetra_Import(*lev_map_loc, *lev_map_dist));
    Teuchos::RCP<Epetra_Import> intlev_loc2dist =
        Teuchos::rcp(new Epetra_Import(*intlev_map_loc, *intlev_map_dist));

    CHECK_ZERO(tatm_loc->Import(*tatm_dist, *lev_loc2dist, Insert));
    CHECK_ZERO(emip_loc->Import(*emip_dist, *lev_loc2dist, Insert));
    CHECK_ZERO(temp_loc->Import(*temp_dist, *intlev_loc2dist, Insert));
    CHECK_ZERO(salt_loc->Import(*salt_dist, *intlev_loc2dist, Insert));
    CHECK_ZERO(spert_loc->Import(*spert_dist, *lev_loc2dist, Insert));
    double tatmmax, emipmax;
    double tatmmin, emipmin;
    CHECK_ZERO(tatm_dist->MaxValue(&tatmmax));
    CHECK_ZERO(emip_dist->MaxValue(&emipmax));
    CHECK_ZERO(tatm_dist->MinValue(&tatmmin));
    CHECK_ZERO(emip_dist->MinValue(&emipmin));

    INFO("Temperature forcing from data ranges between: ["
         << tatmmin <<".." << tatmmax << "]");
    INFO("Salinity forcing from data ranges between: ["
         << emipmin <<".." <<emipmax <<"]");

////////////////////////////////////////////////////////////////////////////////

    // initialize THCM subdomain
    DEBUG("call init..."); // in usrc.F90
    FNAME(init)(&nloc, &mloc, &lloc, &nmlglob,
                &xminloc, &xmaxloc, &yminloc, &ymaxloc,
                &alphaT,&alphaS,
                &ih, &vmix_, &tap, &irho_mixing,
                &coriolis_on,
                &perio, landm,
                taux, tauy, tatm, emip, spert);

    INFO("   initialize THCM subdomain done");

    if (internal_forcing_)
    {
        F90NAME(m_usr,set_internal_forcing)(temp,salt);
    }

    // get a map object for constructing vectors without overlap
    // (load-balanced, used for solve phase)
    solveMap_ = domain_->GetSolveMap();

    // Create internal vectors
    initialSolution_ = Teuchos::rcp(new Epetra_Vector(*solveMap_));
    diagB_           = Teuchos::rcp(new Epetra_Vector(*solveMap_));
    localDiagB_      = Teuchos::rcp(new Epetra_Vector(*standardMap_));
    localRhs_        = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));
    localSol_        = Teuchos::rcp(new Epetra_Vector(*assemblyMap_));

    // 2D overlapping interface fields
    localAtmosT_     = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosQ_     = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosA_     = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localAtmosP_     = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localSeaiceQ_    = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localSeaiceM_    = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localSeaiceG_    = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localOceanE_     = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localEmip_       = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localSurfTmp_    = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));
    localTatm_       = Teuchos::rcp(new Epetra_Vector(*assemblySurfaceMap_));

    // allocate mem for the CSR matrix in THCM.
    // first ask how big it should be:
    int nrows, nnz;
    DEBUG("call get_array_sizes...");
    F90NAME(m_mat,get_array_sizes)(&nrows,&nnz);
    INFO("Allocating Fortran CSR arrays, nrows=" << nrows << ", nnz=" << nnz);

    // allocate the memory
    begA_ = new int[nrows+1];
    coA_  = new double[nnz];
    jcoA_ = new int[nnz];

    coB_  = new double[nrows];

    begF_ = new int[nrows+1];
    coF_  = new double[nrows];
    jcoF_ = new int[nrows];

    // give THCM the opportunity to set its pointers to
    // the new memory block
    DEBVAR("call set_pointers...");
    F90NAME(m_mat,set_pointers)(&nrows,&nnz,begA_,jcoA_,coA_,coB_,begF_,jcoF_,coF_);

    // Initialize integral condition row, correction and coefficients
    rowintcon_     = -1;
    intCorrection_ = 0.0;

    // Initialize salinity flux correction
    scorr_ = 0.0;

    // Obtain integral condition coefficients

    int N = domain_->GlobalN();
    int M = domain_->GlobalM();
    int L = domain_->GlobalL();

    int Nic = paramList_.get<int>("Integral row coordinate i");
    if (Nic == -1) {
        Nic = N-1;
    }
    int Mic = paramList_.get<int>("Integral row coordinate j");
    if (Mic == -1) {
        Mic = M-1;
    }
    int midx;     // mask index
    int mval = 1; // mask value
    int tmp  = 0;
    if (sres_ == 0)
    {
        if (comm->MyPID() == 0)
        {
            midx = FIND_ROW2(1, N+2, M+2, L+2, Nic + 1, Mic + 1, L, 1);
            tmp = (*landm_glb)[midx];
        }
        comm->SumAll(&tmp, &mval, 1);

        if (mval != 0) // throw ERROR if integral row on land
        {
            if (comm->MyPID() == 0)
            {
                for (int j = M+1; j >= 0; --j)
                {
                    for (int i = 0; i != N+2; ++i)
                    {
                        midx = FIND_ROW2(1, N+2, M+2, L+2, i, j, L, 1);
                        std::cout << (*landm_glb)[midx];
                    }
                    std::cout << std::endl;
                }
            }

            ERROR("Integral row coordinates ("
                  << Nic << "," << Mic << ") give a land point! \n"
                  << "  Please give better coordinates in xml.",
                  __FILE__, __LINE__);
        }

        // location of integral condition
        rowintcon_ = FIND_ROW2(_NUN_,N,M,L,Nic,Mic,L-1,SS);

        INFO("THCM: integral condition for S is in global row " << rowintcon_);

    }

    // Initialize integral coefficients
    intcondCoeff_ = Teuchos::rcp(new Epetra_Vector(*solveMap_));

    // Obtain integral coefficients
    getIntCondCoeff();

    // create a graph describing the maximal matrix pattern.
    // Note that in LOCA we can't change the pattern of the matrix
    // during the continuation process as we pass pointers to LOCA
    // which can't be changed, whereas changing the pattern would
    // require building a whole new matrix. As we can't predict
    // where convective adjustment will happen, we assume it hap-
    // pens everywhere.
    Teuchos::RCP<Epetra_CrsGraph> localMatrixGraph = CreateMaximalGraph();
    Teuchos::RCP<Epetra_CrsGraph> testMatrixGraph  = CreateMaximalGraph(false);
    Teuchos::RCP<Epetra_CrsGraph> matrixGraph      = localMatrixGraph;

    if (solveMap_!=standardMap_)
    {
        matrixGraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*solveMap_,20));
        Teuchos::RCP<Epetra_Import> import = Teuchos::rcp(new Epetra_Import(*standardMap_,*solveMap_));
        DEBUG("Migrate graph to solve map...");
        CHECK_ZERO(matrixGraph->Export(*localMatrixGraph,*import,Insert));
        CHECK_ZERO(matrixGraph->FillComplete());
    }

    localJac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *localMatrixGraph));
    localJac_->SetLabel("Local Jacobian");

    testJac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *testMatrixGraph));
    testJac_->SetLabel("Testing Jacobian");

    jac_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *matrixGraph));
    jac_->SetLabel("Jacobian");

    localFrc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *standardMap_, 1));
    localFrc_->SetLabel("Local Forcing");
    frc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *solveMap_, 1));
    frc_->SetLabel("Forcing");

// we can select two points where the continuity equation will be replaced by
// P(i,j,k) = 0. This is experimental, we hope to fix the divergence problem in the 4D case
//  like this
//      int N = domain_->GlobalN();
//      int M = domain_->GlobalM();
//      int L = domain_->GlobalL();
//

    if (fixPressurePoints_)
    {
        rowPfix1_ = FIND_ROW2(_NUN_,N,M,L,N-1,M-1,L-1,PP);
        rowPfix2_ = FIND_ROW2(_NUN_,N,M,L,N-2,M-1,L-1,PP);
    }
    else
    {
        rowPfix1_ = -1;
        rowPfix2_ = -1;
    }

    // build vector with integral coefficients
    this->evaluateB();

    scalingType_ = paramList_.get<std::string>("Scaling");

    if (scalingType_ == "THCM")
    {
        // construct the scaling object. The scaling is computed by THCM (m_scaling)
        // and passed on to Trilinos:
        rowScaling_ = Teuchos::rcp(new Epetra_Vector(*solveMap_));
        rowScaling_->SetLabel("Row Scaling");
        colScaling_ = Teuchos::rcp(new Epetra_Vector(*solveMap_));
        colScaling_->SetLabel("Col Scaling");
        localRowScaling_ = Teuchos::rcp(new Epetra_Vector(*assemblyMap_) );
        localColScaling_ = Teuchos::rcp(new Epetra_Vector(*assemblyMap_) );

        rowScaling_->PutScalar(1.0);
        colScaling_->PutScalar(1.0);
    }

    Teuchos::ParameterList tmpList = paramList_.sublist("Starting Parameters");

    for (auto& pair : tmpList) {
        std::string key = pair.first;
        double val = pair.second.getValue(&val);
        if (!std::isnan(val)) setParameter(key, val);
        else {
            getParameter(key, val);
            paramList_.sublist("Starting Parameters").remove(key);
        }

        paramList_.sublist("Starting Parameters").get(key, val);
    }

    params = paramList_;
}

//=============================================================================
THCM::~THCM()
{
    INFO("THCM destructor");
    FNAME(finalize)();
    if (comm_->MyPID()==0)
    {
        F90NAME(m_global,finalize)();
    }

    delete [] jcoA_;
    delete [] coA_;
    delete [] begA_;

    delete [] coB_;

    delete [] jcoF_;
    delete [] coF_;
    delete [] begF_;
    // the rest is handled by Teuchos::rcp's
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getSolution()
{
    return initialSolution_;
}

//=============================================================================
Teuchos::RCP<Epetra_CrsMatrix> THCM::getJacobian()
{
    return jac_;
}

//=============================================================================
// Compute and get the forcing
bool THCM::computeForcing()
{
    if (localFrc_->Filled())
        localFrc_->PutScalar(0.0); // set all matrix entries to zero

    const int maxlen = 1; // only 1 element per row
    int indices[maxlen];
    double values[maxlen];

    int index, numentries;

    int NumMyElements = assemblyMap_->NumMyElements();
    int imax = NumMyElements;
    std::vector<int> MyElements;

    FNAME(stochastic_forcing)();

    for (int i = 0; i < imax; i++)
    {
#ifndef NO_INTCOND
        if (sres_ == 0 && assemblyMap_->GID(i) == rowintcon_)
            continue;
#endif
        if (domain_->IsGhost(i, _NUN_))
            continue;

        index = begF_[i]; // note that these arrays use 1-based indexing
        numentries = begF_[i+1] - index;
        for (int j = 0; j <  numentries ; j++)
        {
            indices[j] = assemblyMap_->GID(jcoF_[index-1+j] - 1);
            values[j]  = coF_[index - 1 + j];
            MyElements.push_back(indices[j]);
        }

        if (localFrc_->Filled())
        {
            int ierr = localFrc_->ReplaceGlobalValues(assemblyMap_->GID(i), numentries,
                                                     values, indices);

            // ierr == 3 probably means not all row entries are replaced,
            // does not matter because we zeroed them.
            if (((ierr!=0) && (ierr!=3)))
            {
                std::cout << "\n ERROR " << ierr;
                std::cout << "\n myPID " << comm_->MyPID();
                std::cout <<"\n while inserting/replacing values in local Jacobian"
                          << std::endl;

                INFO(" ERROR while inserting/replacing values in local Jacobian");

                int GRID = assemblyMap_->GID(i);
                std::cout << " GRID: " << GRID << std::endl;
                std::cout << " max GRID: " << assemblyMap_->GID(imax-1) << std::endl;
                std::cout << " number of entries: " << numentries << std::endl;

                std::cout << " entries: ";
                for (int j = 0; j < numentries; j++)
                    std::cout << "(" << indices[j] << " " << values[j] << ") ";
                std::cout << std::endl;

                std::cout << " NumMyElements:        " << NumMyElements << std::endl;
                std::cout << " i:                    " << i << std::endl;
                std::cout << " imax:                 " << imax << std::endl;
                std::cout << " maxlen:               " << maxlen << std::endl;

                std::cout << " row:                  " << GRID << std::endl;
                int LRID = localFrc_->LRID(GRID);
                std::cout << " LRID:                 " << LRID << std::endl;
                std::cout << " graph inds in LRID:   "
                          << localFrc_->Graph().NumMyIndices(LRID) << std::endl;

                int ierr2 = localFrc_->ExtractGlobalRowCopy(
                    assemblyMap_->GID(i), maxlen, numentries, values, indices);

                std::cout << "\noriginal row: " << std::endl;
                std::cout << "number of entries: " << numentries << std::endl;
                std::cout << "entries: ";

                for (int j=0; j < numentries; j++)
                    std::cout << "(" << indices[j] << " " << values[j] << ") ";
                std::cout << std::endl;

                CHECK_ZERO(ierr2);
            }
        }
        else
        {
            CHECK_ZERO(localFrc_->InsertGlobalValues(assemblyMap_->GID(i), numentries,
                                                    values, indices));
        }
    } //i-loop over rows

    auto last = std::unique(MyElements.begin(), MyElements.end());
    Epetra_Map colMap(-1, (int)std::distance(MyElements.begin(), last),
                      &MyElements[0], 0, *comm_);
    CHECK_ZERO(localFrc_->FillComplete(colMap, *standardMap_));

    // redistribute according to solveMap_ (may be load-balanced)
    // standard and solve maps are equal
    domain_->Standard2Solve(*localFrc_, *frc_);     // no effect
    CHECK_ZERO(frc_->FillComplete(colMap, *solveMap_));

    return true;
}

Teuchos::RCP<Epetra_CrsMatrix> THCM::getForcing()
{
    return frc_;
}

//=============================================================================
// Compute Jacobian and/or RHS.
bool THCM::evaluate(const Epetra_Vector& soln,
                    Teuchos::RCP<Epetra_Vector> tmp_rhs,
                    bool computeJac,
                    bool maskTest)
{
    if (compSalInt_)
    {
        double intcond;
        CHECK_ZERO(intcondCoeff_->Dot(soln,&intcond));
        intcond -= intCorrection_; // apply correction

        // INFO("Salinity integral condition: " << intcond);

    }

    if (!(soln.Map().SameAs(*solveMap_)))
    {
        ERROR("Map of solution vector not same as solve-map ",__FILE__,__LINE__);
    }


    // convert to standard distribution and
    // import values from ghost-nodes on neighbouring subdomains:
    domain_->Solve2Assembly(soln,*localSol_);


    int NumMyElements = assemblyMap_->NumMyElements();

//  DEBUG("=== evaluate: input vector");
//  DEBUG( (domain_->Gather(*soln,0)) )

    // extract an array to pass on to THCM:
    // We use 'Copy' mode, which is clean but possibly slow.
    // Probably 'View' would be allowable as well.
    double* solution;
    localSol_->ExtractView(&solution);
    if(tmp_rhs!=Teuchos::null)
    {
        // INFO("Compute RHS...");
        // build rhs simultaneously on each process
        double* RHS;
        CHECK_ZERO(localRhs_->ExtractView(&RHS));
        TIMER_START("Ocean: compute rhs: fortran part");
        // compute right-hand-side on whole subdomain (by THCM)
        FNAME(rhs)(solution, RHS);
        TIMER_STOP("Ocean: compute rhs: fortran part");

        // export overlapping rhs to unique-id global rhs vector,
        // and load-balance for solve phase:
        domain_->Assembly2Solve(*localRhs_,*tmp_rhs);

        // Negating the RHS... Instead of scaling the RHS, scaling the
        // Jacobian corresponds better to how the equations would be
        // written.
        CHECK_ZERO(tmp_rhs->Scale(-1.0));

#ifndef NO_INTCOND
        if ((sres_ == 0) && !maskTest)
        {
            double intcond;
            //TODO: check which is better:
            CHECK_ZERO(intcondCoeff_->Dot(soln,&intcond));
            //std::cout << " dot product: "<<intcond << std::endl;
            //intcond = 0.0;
            if (tmp_rhs->Map().MyGID(rowintcon_))
            {
                (*tmp_rhs)[tmp_rhs->Map().LID(rowintcon_)] =
                    intSign_ * (intcond - intCorrection_);
            }
        }
#endif
        if (rowPfix1_>=0)
        {
            if (tmp_rhs->Map().MyGID(rowPfix1_))
            {
                (*tmp_rhs)[tmp_rhs->Map().LID(rowPfix1_)]=0.0;
            }
        }
        if (rowPfix2_>=0)
        {
            if (tmp_rhs->Map().MyGID(rowPfix2_))
            {
                (*tmp_rhs)[tmp_rhs->Map().LID(rowPfix2_)]=0.0;
            }
        }
    }
    if(computeJac)
    {
        // INFO("Compute Jacobian...");
        Teuchos::RCP<Epetra_CrsMatrix> tmpJac;
        if (maskTest) // Use Jacobian based on testing graph
            tmpJac = testJac_;
        else // Use Jacobian based on standard graph
            tmpJac = localJac_;

        tmpJac->PutScalar(0.0); // set all matrix entries to zero
        localDiagB_->PutScalar(0.0);

        //Call the fortran routine, providing the solution vector,
        //and get back the three vectors of the sparse Jacobian (CSR form)
        TIMER_START("Ocean: compute jacobian: fortran part");

        // If we test the mask we need non-restoring conditions in the matrix
        if (maskTest)
        {
            int tmp_sres = 0;
            FNAME(setsres)(&tmp_sres);
        }

        FNAME(matrix)(solution);

        // Restore from the testing config
        if (maskTest)
            FNAME(setsres)(&sres_);

        TIMER_STOP("Ocean: compute jacobian: fortran part");

        const int maxlen = _NUN_*_NP_+1;    //nun*np+1 is max nonzeros per row
        int indices[maxlen];
        double values[maxlen];

        int index, numentries;

        int imax = NumMyElements;

        for (int i = 0; i < imax; i++)
        {
            if (!domain_->IsGhost(i, _NUN_) &&
                ( ( assemblyMap_->GID(i) != rowintcon_ ) || maskTest ) )
            {
                index = begA_[i]; // note that these arrays use 1-based indexing
                numentries = begA_[i+1] - index;
                for (int j = 0; j <  numentries ; j++)
                {
                    indices[j] = assemblyMap_->GID(jcoA_[index-1+j] - 1);
                    values[j]  = coA_[index - 1 + j];
                }

                int ierr = tmpJac->ReplaceGlobalValues(assemblyMap_->GID(i), numentries,
                                                         values, indices);

                // ierr == 3 probably means not all row entries are replaced,
                // does not matter because we zeroed them.
                if (((ierr!=0) && (ierr!=3)))
                {
                    std::stringstream ss;
                    ss << "graph_pid" << comm_->MyPID();
                    std::ofstream file(ss.str());
                    file << tmpJac->Graph();

                    std::cout << "\n ERROR " << ierr;
                    std::cout << ((ierr == 2) ? ": value excluded" : "") << std::endl;
                    std::cout << "\n myPID " << comm_->MyPID();
                    std::cout <<"\n while inserting/replacing values in local Jacobian"
                              << std::endl;

                    INFO(" ERROR while inserting/replacing values in local Jacobian");

                    int GRID = assemblyMap_->GID(i);
                    std::cout << " GRID: " << GRID << std::endl;
                    std::cout << " max GRID: " << assemblyMap_->GID(imax-1) << std::endl;
                    std::cout << " number of entries: " << numentries << std::endl;

                    std::cout << " entries: ";
                    for (int j = 0; j < numentries; j++)
                        std::cout << "(" << indices[j] << " " << values[j] << ") ";
                    std::cout << std::endl;

                    std::cout << " NumMyElements:        " << NumMyElements << std::endl;
                    std::cout << " i:                    " << i << std::endl;
                    std::cout << " imax:                 " << imax << std::endl;
                    std::cout << " maxlen:               " << maxlen << std::endl;

                    std::cout << " row:                  " << GRID << std::endl;
                    std::cout << " have rowintcon:       " << tmpJac->MyGRID(rowintcon_)
                              << std::endl;
                    std::cout << " rowintcon:            " << rowintcon_ << std::endl;
                    std::cout << " assembly rowintcon:   " << assemblyMap_->LID(rowintcon_)
                              << std::endl;
                    std::cout << " standard rowintcon:   " << standardMap_->LID(rowintcon_)
                              << std::endl;
                    int LRID = tmpJac->LRID(GRID);
                    std::cout << " LRID:                 " << LRID << std::endl;
                    std::cout << " graph inds in LRID:   "
                              << tmpJac->Graph().NumMyIndices(LRID) << std::endl;

                    int ierr2 = tmpJac->ExtractGlobalRowCopy
                        (assemblyMap_->GID(i), maxlen, numentries, values, indices);

                    std::cout << "\noriginal row: " << std::endl;
                    std::cout << "number of entries: " << numentries << std::endl;
                    std::cout << "entries: ";

                    for (int j=0; j < numentries; j++)
                        std::cout << "(" << indices[j] << " " << values[j] << ") ";
                    std::cout << std::endl;

                    CHECK_ZERO(ierr2);
                }

                // reconstruct the diagonal matrix B
                int lid = standardMap_->LID(assemblyMap_->GID(i));
                double mass_param = 1.0;
                (*localDiagB_)[lid] = coB_[i] * mass_param;
            } //not a ghost?
        } //i-loop over rows

#ifndef NO_INTCOND
        if ((sres_ == 0) && !maskTest)
        {
            this->intcond_S(*tmpJac,*localDiagB_);
        }
#endif

        if (fixPressurePoints_)
            this->fixPressurePoints(*tmpJac,*localDiagB_);

        CHECK_ZERO(tmpJac->FillComplete());

        // redistribute according to solveMap_ (may be load-balanced)
        // standard and solve maps are equal
        domain_->Standard2Solve(*localDiagB_, *diagB_); // no effect
        domain_->Standard2Solve(*tmpJac, *jac_);     // no effect
        CHECK_ZERO(jac_->FillComplete());

        if (scalingType_ == "THCM")
        {
            DEBUG(" THCM:  RecomputeScaling()");
            this->RecomputeScaling();

#if 0
            std::ofstream ofs1("rowScaling_.txt");
            ofs1 << *rowScaling_;
            ofs1.close();
            std::ofstream ofs2("colScaling_.txt");
            ofs2 << *colScaling_;
            ofs2.close();
#endif
        }

    } // matrix
    return true;
}

// just reconstruct the diagonal matrix B from THCM
void THCM::evaluateB(void)
{
    int NumMyElements = assemblyMap_->NumMyElements();

    DEBUG("Construct matrix B...");

    localDiagB_->PutScalar(0.0);
    FNAME(fillcolb)();
    for (int i = 0; i < NumMyElements; i++)
    {
        if (!domain_->IsGhost(i, _NUN_))
        {
            // reconstruct the diagonal matrix B
            int lid = standardMap_->LID(assemblyMap_->GID(i));
            (*localDiagB_)[lid] = coB_[i];
        } // not a ghost?
    } // i-loop over rows

    if (fixPressurePoints_)
    {
        for (int i=1;i<=2;i++) {
            int row = (i==1)? rowPfix1_: rowPfix2_;
            if (localDiagB_->Map().MyGID(row))
            {
                int lid = localDiagB_->Map().LID(row);
                (*localDiagB_)[lid] = 0.0; // no more time-dependence for this P-point
            }
        }
    }

#ifndef NO_INTCOND
    if (sres_ == 0)
    {
        if (localDiagB_->Map().MyGID(rowintcon_))
        {
            (*localDiagB_)[localDiagB_->Map().LID(rowintcon_)]=0.0;
        }
    }
#endif
    domain_->Standard2Solve(*localDiagB_,*diagB_);
}

//==================================================================
// Get current global landmask including borders
std::shared_ptr<std::vector<int> > THCM::getLandMask()
{
    // length of landmask array
    int dim = (n_+2)*(m_+2)*(l_+2);

    // Create landmask array
    std::shared_ptr<std::vector<int> > landm =
        std::make_shared<std::vector<int> >(dim, 0);

    // Let THCM fill the landmask array on proc = 0
    if (comm_->MyPID() == 0)
        F90NAME(m_global, get_current_landm)(&(*landm)[0]);

    // Get the MpiComm from Epetra
    Epetra_MpiComm const MpiComm =
        dynamic_cast<Epetra_MpiComm const &>(*comm_);

    // Broadcast the landmask
    MPI_Bcast(&(*landm)[0], dim, MPI_INTEGER, 0, MpiComm.GetMpiComm());

    return landm;
}

//=============================================================================
// Get distributed land mask based on maskName
Teuchos::RCP<Epetra_IntVector> THCM::getLandMask(std::string const &maskName,
                                                 Teuchos::RCP<Epetra_Vector> fix)
{
    // Create gathered map for land mask
    // All indices are on root process
    int I0 = 0; int I1 = n_+1;
    int J0 = 0; int J1 = m_+1;
    int K0 = 0; int K1 = l_+1;

    int i0 = 0, i1 = -1, j0 = 0, j1 = -1,k0 = 0,k1 = -1;
    if (comm_->MyPID() == 0)
    {
        i1 = I1; j1 = J1; k1=K1;
    }

    Teuchos::RCP<Epetra_Map> landmap_glb =
        Utils::CreateMap(i0,i1,j0,j1,k0,k1,I0,I1,J0,J1,K0,K1,*comm_);

    // Create sequential landmask array on proc 0
    Teuchos::RCP<Epetra_IntVector> landm_glb =
        Teuchos::rcp(new Epetra_IntVector(*landmap_glb));

    // Get global landmask from fortran
    int *landm;
    if (comm_->MyPID()==0)
    {
        F90NAME(m_global,set_maskfile)(maskName.c_str());

        CHECK_ZERO(landm_glb->ExtractView(&landm));

        // Let THCM fill the global landm array and put it into our C pointer location
        if (maskName == "current")
            F90NAME(m_global,get_current_landm)(landm);
        else
            F90NAME(m_global,get_landm)(landm);
    }

    // Fixing landmask
    if (fix != Teuchos::null)
    {
        // Gather fix on proc 0
        Teuchos::RCP<Epetra_MultiVector> fix0 =
            Utils::Gather(*fix, 0);

        int len = fix0->Map().NumMyElements();

        std::vector<double> fix1(len);
        (*fix0)(0)->ExtractCopy(&fix1[0]);

        int i,j,k;
        if (comm_->MyPID() == 0 && len > 0)
        {
            int pos = 0;
            int idx = 0;
            for (k = K0+1; k < K1; ++k)
            {
                for (j = J0+1; j < J1; ++j)
                {
                    for (i = I0+1; i < I1; ++i)
                    {
                        idx = k*(m_+2)*(n_+2) + j*(n_+2) + i;

                        // magic number 2 indicates too little
                        // contributions in the corresponding matrix
                        // row
                        if (fix1[pos] == 2)
                            landm[idx] = 1; // adjust landmask

                        pos++;
                    }
                }
            }

            // Setting fixed global landmask in THCM
            F90NAME(m_global, set_landm)(landm);
        }
    }

    // Return distributed landmask
    Teuchos::RCP<Epetra_IntVector> landm_loc = distributeLandMask(landm_glb);

    return landm_loc;
}

//=============================================================================
// Set distributed landmask in THCM
// set_landmask takes care of a few reinitializations if requested
void THCM::setLandMask(Teuchos::RCP<Epetra_IntVector> landmask, bool init)
{
    // in the main part of THCM (except m_global) we set periodic
    // boundary conditions to .false. _unless_ we are running a
    // periodic problem on a single CPU in the x-direction:
    Teuchos::RCP<Epetra_Comm> xComm = domain_->GetProcRow(0);
    int perio   = (periodic_ && xComm->NumProc() == 1);

    int *landm;
    CHECK_ZERO(landmask->ExtractView(&landm));

    int reinit = (init) ? 1 : 0;
    FNAME(set_landmask)(landm, &perio, &reinit);
}

//=============================================================================
// Set global landmask in THCM
void THCM::setLandMask(std::shared_ptr<std::vector<int> > landmask)
{
    if (comm_->MyPID() == 0)
        F90NAME(m_global, set_landm)(&(*landmask)[0]);
}

//=============================================================================
void THCM::setAtmosphereT(Teuchos::RCP<Epetra_Vector> const &atmosT)
{
    CHECK_MAP(atmosT, standardSurfaceMap_);
    // Standard2Assembly
    // Import atmosT into local atmosT
    CHECK_ZERO(localAtmosT_->Import(*atmosT, *as2std_surf_, Insert));

    double *locAtmosT;

    localAtmosT_->ExtractView(&locAtmosT);
    F90NAME(m_inserts, insert_atmosphere_t)( locAtmosT );
}

//=============================================================================
void THCM::setAtmosphereQ(Teuchos::RCP<Epetra_Vector> const &atmosQ)
{
    CHECK_MAP(atmosQ, standardSurfaceMap_);

    // Standard2Assembly
    // Import atmosQ into local atmosQ
    CHECK_ZERO( localAtmosQ_->Import(*atmosQ, *as2std_surf_, Insert) );

    double *tmpAtmosQ;
    localAtmosQ_->ExtractView(&tmpAtmosQ);
    F90NAME(m_inserts, insert_atmosphere_q)( tmpAtmosQ );
}

//=============================================================================
void THCM::setAtmosphereA(Teuchos::RCP<Epetra_Vector> const &atmosA)
{
    CHECK_MAP(atmosA, standardSurfaceMap_);

    // Standard2Assembly
    // Import atmosQ into local atmosQ
    CHECK_ZERO( localAtmosA_->Import(*atmosA, *as2std_surf_, Insert) );

    double *tmpAtmosA;
    localAtmosA_->ExtractView(&tmpAtmosA);
    F90NAME(m_inserts, insert_atmosphere_a)( tmpAtmosA );
}

//=============================================================================
void THCM::setAtmosphereP(Teuchos::RCP<Epetra_Vector> const &atmosP)
{
    CHECK_MAP(atmosP, standardSurfaceMap_);

    // Import atmosP into local atmosP
    CHECK_ZERO(localAtmosP_->Import(*atmosP, *as2std_surf_, Insert));

    double *tmpAtmosP;
    localAtmosP_->ExtractView(&tmpAtmosP);

    F90NAME(m_inserts, insert_atmosphere_p)( tmpAtmosP );
}

//=============================================================================
void THCM::setSeaIceQ(Teuchos::RCP<Epetra_Vector> const &seaiceQ)
{
    CHECK_MAP(seaiceQ, standardSurfaceMap_);
    CHECK_ZERO(localSeaiceQ_->Import(*seaiceQ, *as2std_surf_ ,Insert));
    double *Q;
    localSeaiceQ_->ExtractView(&Q);
    F90NAME(m_inserts, insert_seaice_q)( Q );
}

//=============================================================================
void THCM::setSeaIceM(Teuchos::RCP<Epetra_Vector> const &seaiceM)
{
    CHECK_MAP(seaiceM, standardSurfaceMap_);
    CHECK_ZERO(localSeaiceM_->Import(*seaiceM, *as2std_surf_ ,Insert));
    double *M;

    if (!coupledM_)
        localSeaiceM_->PutScalar(0.0); // disable coupling with mask

    localSeaiceM_->ExtractView(&M);
    F90NAME(m_inserts, insert_seaice_m)( M );
}

//=============================================================================
void THCM::setSeaIceG(Teuchos::RCP<Epetra_Vector> const &seaiceG)
{
    CHECK_MAP(seaiceG, standardSurfaceMap_);
    CHECK_ZERO(localSeaiceG_->Import(*seaiceG, *as2std_surf_ ,Insert));
    double *G;
    localSeaiceG_->ExtractView(&G);
    F90NAME(m_inserts, insert_seaice_g)( G );
}

//=============================================================================
//FIXME: superfluous?? ->setAtmosphereT()
void THCM::setTatm(Teuchos::RCP<Epetra_Vector> const &tatm)
{

    if (!(tatm->Map().SameAs(*standardSurfaceMap_)))
    {
        // INFO("THCM::setAtmosphereEP: atmosP map -> standardSurfaceMap_");
        CHECK_ZERO(tatm->ReplaceMap(*standardSurfaceMap_));
    }

    // Standard2Assembly
    // Import atmosP into local atmosP
    CHECK_ZERO(localTatm_->Import(*tatm, *as2std_surf_, Insert));

    double *tmpTatm;
    localTatm_->ExtractView(&tmpTatm);

    F90NAME(m_inserts, insert_tatm)( tmpTatm );
}

//=============================================================================
void THCM::setEmip(Teuchos::RCP<Epetra_Vector> const &emip, char mode)
{

    if (!(emip->Map().SameAs(*standardSurfaceMap_)))
    {
        INFO("THCM::setEmip: emip map -> standardSurfaceMap_");
        CHECK_ZERO(emip->ReplaceMap(*standardSurfaceMap_));
    }

    // Standard2Assembly
    // Import atmosP into local atmosP
    CHECK_ZERO(localSurfTmp_->Import(*emip, *as2std_surf_, Insert));

    double *tmpEmip;
    localSurfTmp_->ExtractView(&tmpEmip);

    if (mode == 'A')
    {
        F90NAME(m_inserts, insert_adapted_emip )( tmpEmip );
    }
    else if (mode == 'P')
    {
        F90NAME(m_inserts, insert_emip_pert )( tmpEmip );
    }
    else
    {
        F90NAME(m_inserts, insert_emip )( tmpEmip );
    }
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getSunO()
{
    double *suno;
    Epetra_Vector localSunO(*assemblySurfaceMap_);
    localSunO.ExtractView(&suno);
    F90NAME(m_probe, get_suno) ( suno );
    Teuchos::RCP<Epetra_Vector> out =
        Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    CHECK_ZERO(out->Export(localSunO, *as2std_surf_, Zero));
    return out;
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getEmip(char mode)
{
    double* tmpEmip;
    localSurfTmp_->ExtractView(&tmpEmip);

    if (mode == 'A')
    {
        F90NAME(m_probe, get_adapted_emip )( tmpEmip );
    }
    else if (mode == 'P')
    {
        F90NAME(m_probe, get_emip_pert )( tmpEmip );
    }
    else
    {
        F90NAME(m_probe, get_emip )( tmpEmip );
    }


    Teuchos::RCP<Epetra_Vector> emip =
        Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

    // Export assembly map surface emip
    CHECK_ZERO(emip->Export(*localEmip_, *as2std_surf_, Zero));

    return emip;
}

//=============================================================================
std::vector<Teuchos::RCP<Epetra_Vector> > THCM::getFluxes()
{
    int numFluxes = _MSI+1;

    std::vector<Teuchos::RCP<Epetra_Vector> > fluxes;
    std::vector<Epetra_Vector> localFluxes;

    for (int i = 0; i != numFluxes; ++i)
    {
        fluxes.push_back(Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_)));
        localFluxes.push_back(Epetra_Vector(*assemblySurfaceMap_));
    }

    std::vector<double*> tmpPtrs(numFluxes);

    for (int i = 0; i != numFluxes; ++i)
        localFluxes[i].ExtractView(&tmpPtrs[i]);

    double* solution;
    localSol_->ExtractView(&solution);

    F90NAME(m_probe, get_salflux )( solution, tmpPtrs[_Sal],  &scorr_,
                                    tmpPtrs[_QSOA], tmpPtrs[_QSOS] );

    F90NAME(m_probe, get_temflux )( solution, tmpPtrs[_Temp], tmpPtrs[_QSW],
                                    tmpPtrs[_QSH], tmpPtrs[_QLH], tmpPtrs[_QTOS],
                                    tmpPtrs[_MSI] );

    for (int i = 0; i != numFluxes; ++i)
    {
        CHECK_ZERO(fluxes[i]->Export(localFluxes[i], *as2std_surf_, Zero));
    }

    return fluxes;
}

//=============================================================================
THCM::Derivatives THCM::getDerivatives()
{
    Derivatives d;

    double *solution;
    localSol_->ExtractView(&solution);

    Epetra_Vector local_dftdm(*assemblySurfaceMap_);
    Epetra_Vector local_dfsdq(*assemblySurfaceMap_);
    Epetra_Vector local_dfsdm(*assemblySurfaceMap_);
    Epetra_Vector local_dfsdg(*assemblySurfaceMap_);

    double *dftdm, *dfsdq, *dfsdm, *dfsdg;
    local_dftdm.ExtractView(&dftdm);
    local_dfsdq.ExtractView(&dfsdq);
    local_dfsdm.ExtractView(&dfsdm);
    local_dfsdg.ExtractView(&dfsdg);
    F90NAME(m_probe, get_derivatives)( solution, dftdm, dfsdq, dfsdm, dfsdg);

    d.dFTdM = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    d.dFSdQ = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    d.dFSdM = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));
    d.dFSdG = Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

    d.dFTdM->Export(local_dftdm, *as2std_surf_, Zero);
    d.dFSdQ->Export(local_dfsdq, *as2std_surf_, Zero);
    d.dFSdM->Export(local_dfsdm, *as2std_surf_, Zero);
    d.dFSdG->Export(local_dfsdg, *as2std_surf_, Zero);

    return d;
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getLocalAtmosT()
{
    double *tmpAtmosT;
    localAtmosT_->ExtractView(&tmpAtmosT);
    F90NAME(m_probe, get_atmosphere_t )( tmpAtmosT );
    return localAtmosT_;
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getLocalAtmosQ()
{
    double *tmpAtmosQ;
    localAtmosP_->ExtractView(&tmpAtmosQ);
    F90NAME(m_probe, get_atmosphere_q )( tmpAtmosQ );
    return localAtmosQ_;
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getAtmosQ()
{
    Teuchos::RCP<Epetra_Vector> atmosQ =
        Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

    // Export assembly map surface evaporation to standard surface map
    CHECK_ZERO(atmosQ->Export(*getLocalAtmosQ(), *as2std_surf_, Zero));
    return atmosQ;
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getLocalAtmosP()
{
    double *tmpAtmosP;
    localAtmosP_->ExtractView(&tmpAtmosP);
    F90NAME(m_probe, get_atmosphere_p )( tmpAtmosP );
    return localAtmosP_;
}

//============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getLocalOceanE()
{
    double *tmpOceanE;
    localOceanE_->ExtractView(&tmpOceanE);

    // localsol should contain something meaningful
    double* solution;
    localSol_->ExtractView(&solution);

    F90NAME( m_probe, compute_evap )( tmpOceanE, solution );

    return localOceanE_;
}

//============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getOceanE()
{
    Teuchos::RCP<Epetra_Vector> oceanE =
        Teuchos::rcp(new Epetra_Vector(*standardSurfaceMap_));

    // Export assembly map surface evaporation to standard surface map
    CHECK_ZERO(oceanE->Export(*getLocalOceanE(), *as2std_surf_, Zero));
    return oceanE;
}

//=============================================================================
// Recompute scaling for the linear system
void THCM::RecomputeScaling(void)
{

    DEBUG("Compute new scaling...");

    // (1) THCM block scaling

    // array to be filled by F90 routine
    // (average diagonal block, 6x6)
    double ldb[_NUN_*_NUN_];
    double gdb[_NUN_*_NUN_];

    int len = localRowScaling_->MyLength();
    double *rowscal = new double[len];
    double *colscal = new double[len];

    // compute local average diagonal block
    F90NAME(m_scaling,average_block)(ldb);

    comm_->SumAll(ldb,gdb,_NUN_*_NUN_);

    for (int i=0;i<_NUN_*_NUN_;i++)
    {
        gdb[i]/=comm_->NumProc();
    }
    // compute row- and column scaling
    F90NAME(m_scaling,compute)(gdb,rowscal,colscal);

    // put them back in assembly vectors:
    // note: the scaling matrices in Trilinos are
    // defined as the inverse of those in THCM
    for (int i=0;i<len;i++)
    {
        (*localRowScaling_)[i] = 1.0/rowscal[i];
        (*localColScaling_)[i] = 1.0/colscal[i];
    }

    // kick out the ghost nodes:
    domain_->Assembly2Solve(*localRowScaling_,*rowScaling_);
    domain_->Assembly2Solve(*localColScaling_,*colScaling_);

    // make sure T and S are scaled the same way in each cell
    // we need this because of our special block scaling for the
    // ATS matrix in the preconditioner.
    for (int i = TT-1; i < rowScaling_->MyLength(); i += _NUN_)
    {
        double mean = 0.5*((*rowScaling_)[i]+(*rowScaling_)[i+1]);
        (*rowScaling_)[i] = mean;
        (*rowScaling_)[i+1] = mean;
        mean = 0.5*((*colScaling_)[i]+(*colScaling_)[i+1]);
        (*colScaling_)[i] = mean;
        (*colScaling_)[i+1] = mean;
    }

    delete [] rowscal;
    delete [] colscal;
}

//=============================================================================
// convert parameter name to integer
int THCM::par2int(std::string const &label)
{
    // parameter numbering in fortran code
    int TIME   =  0;
    int AL_T   =  1; int RAYL   =  2; int EK_V   =  3; int EK_H   =  4;
    int ROSB   =  5; int MIXP   =  6; int RESC   =  7; int SPL1   =  8;
    int HMTP   =  9; int SUNP   = 10; int PE_H   = 11; int PE_V   = 12;
    int P_VC   = 13; int LAMB   = 14; int SALT   = 15; int WIND   = 16;
    int TEMP   = 17; int BIOT   = 18; int COMB   = 19; int ARCL   = 20;
    int NLES   = 21; int IFRICB = 22; int CONT   = 23; int ENER   = 24;
    int ALPC   = 25; int CMPR   = 26; int FPER   = 27; int SPER   = 28;
    int MKAP   = 29; int SPL2   = 30;

    if      (label == "Time")                            return TIME;
    else if (label == "AL_T")                            return AL_T;
    else if (label == "Rayleigh-Number")                 return RAYL;
    else if (label == "Rossby-Number")                   return ROSB;
    else if (label == "Vertical Ekman-Number")           return EK_V;
    else if (label == "Horizontal Ekman-Number")         return EK_H;
    else if (label == "MIXP")                            return MIXP;
    else if (label == "RESC")                            return RESC;
    else if (label == "SPL1")                            return SPL1;
    else if (label == "Salinity Homotopy")               return HMTP;
    else if (label == "Solar Forcing")                   return SUNP;
    else if (label == "Vertical Peclet-Number")          return PE_V;
    else if (label == "Horizontal Peclet-Number")        return PE_H;
    else if (label == "P_VC")                            return P_VC;
    else if (label == "LAMB")                            return LAMB;
    else if (label == "ARCL")                            return ARCL;
    else if (label == "Salinity Forcing")                return SALT;
    else if (label == "Wind Forcing")                    return WIND;
    else if (label == "Temperature Forcing")             return TEMP;
    else if (label == "Nonlinear Factor")                return BIOT;
    else if (label == "Combined Forcing")                return COMB;
    else if (label == "CONT")                            return CONT;
    else if (label == "IFRICB")                          return IFRICB;
    else if (label == "NLES")                            return NLES;
    else if (label == "CMPR")                            return CMPR;
    else if (label == "ALPC")                            return ALPC;
    else if (label == "Energy")                          return ENER;
    else if (label == "Salinity Perturbation")           return SPER;
    else if (label == "Flux Perturbation")               return FPER;
    else if (label == "MKAP")                            return MKAP;
    else if (label == "SPL2")                            return SPL2;

    return -1;
}

//=============================================================================
// convert parameter name to integer
std::string const THCM::int2par(int index)
{
    std::string label = "Invalid Parameter Index";

    // parameter numbering in fortran code
    int TIME   =  0;
    int AL_T   =  1; int RAYL   =  2; int EK_V   =  3; int EK_H   =  4;
    int ROSB   =  5; int MIXP   =  6; int RESC   =  7; int SPL1   =  8;
    int HMTP   =  9; int SUNP   = 10; int PE_H   = 11; int PE_V   = 12;
    int P_VC   = 13; int LAMB   = 14; int SALT   = 15; int WIND   = 16;
    int TEMP   = 17; int BIOT   = 18; int COMB   = 19; int ARCL   = 20;
    int NLES   = 21; int IFRICB = 22; int CONT   = 23; int ENER   = 24;
    int ALPC   = 25; int CMPR   = 26; int FPER   = 27; int SPER   = 28;
    int MKAP   = 29; int SPL2   = 30;

    if      (index==TIME)   label = "Time";
    else if (index==AL_T)   label = "AL_T";
    else if (index==RAYL)   label = "Rayleigh-Number";
    else if (index==EK_V)   label = "Vertical Ekman-Number";
    else if (index==EK_H)   label = "Horizontal Ekman-Number";
    else if (index==ROSB)   label = "Rossby-Number";
    else if (index==MIXP)   label = "MIXP";
    else if (index==RESC)   label = "RESC";
    else if (index==SPL1)   label = "SPL1";
    else if (index==HMTP)   label = "Salinity Homotopy";
    else if (index==SUNP)   label = "Solar Forcing";
    else if (index==PE_V)   label = "Vertical Peclet-Number";
    else if (index==PE_H)   label = "Horizontal Peclet-Number";
    else if (index==SALT)   label = "Salinity Forcing";
    else if (index==WIND)   label = "Wind Forcing";
    else if (index==TEMP)   label = "Temperature Forcing";
    else if (index==BIOT)   label = "Nonlinear Factor";
    else if (index==COMB)   label = "Combined Forcing";
    else if (index==NLES)   label = "NLES";
    else if (index==ARCL)   label = "ARCL";
    else if (index==IFRICB) label = "IFRICB";
    else if (index==CONT)   label = "CONT";
    else if (index==P_VC)   label = "P_VC";
    else if (index==LAMB)   label = "LAMB";
    else if (index==CMPR)   label = "CMPR";
    else if (index==ALPC)   label = "ALPC";
    else if (index==ENER)   label = "Energy";
    else if (index==MKAP)   label = "MKAP";
    else if (index==SPL2)   label = "SPL2";
    else if (index==FPER)   label = "Flux Perturbation";
    else if (index==SPER)   label = "Salinity Perturbation";
    else
    {
        ERROR("Parameter index is invalid!",__FILE__,__LINE__);
    }
    return label;
}

//=============================================================================
bool THCM::setParameter(std::string label, double value)
{
    int param = par2int(label);
    if (param > 0 && param <= _NPAR_) // time (0) and exp/seas (31/32) are not passed to THCM
    {
        FNAME(setparcs)(&param,&value);
    }
    else
    {
        ERROR("Invalid Parameter",__FILE__,__LINE__);
    }
    return true;
}

//=============================================================================
bool THCM::getParameter(std::string label, double& value)
{
    int param = par2int(label);
    if (param>0 && param<=_NPAR_) // time (0) and exp (_NPAR_+1) are not passed to THCM
    {
        FNAME(getparcs)(&param,&value);
    }
    // The rest is not implemented
    return true;
}

//=============================================================================
bool THCM::writeParams()
{
    FNAME(writeparams)();
    return true;
}

//=============================================================================
Teuchos::RCP<Epetra_IntVector> THCM::distributeLandMask(Teuchos::RCP<Epetra_IntVector> landm_glb)
{

    DEBUG("Create local (land-)maps...");

    // create a non-overlapping distributed map
    int i0 = domain_->FirstRealI()+1; // 'grid-style' indexing is 1-based
    int i1 = domain_->LastRealI()+1;
    int j0 = domain_->FirstRealJ()+1;
    int j1 = domain_->LastRealJ()+1;
    int k0 = domain_->FirstRealK()+1;
    int k1 = domain_->LastRealK()+1;

    // add global boundary cells
    if (i0 == 1) i0-- ;
    if (i1 == n_) i1++;
    if (j0 == 1) j0-- ;
    if (j1 == m_) j1++;
    if (k0 == 1) k0-- ;
    if (k1 == l_) k1++;

    int I0 = 0; int I1 = n_+1;
    int J0 = 0; int J1 = m_+1;
    int K0 = 0; int K1 = l_+1;

    DEBUG("create landmap without overlap...");
    Teuchos::RCP<Epetra_Map> landmap_loc0 = Utils::CreateMap(i0,i1,j0,j1,k0,k1,
                                                             I0,I1,J0,J1,K0,K1,*comm_);

    // create an overlapping distributed map
    i0 = domain_->FirstI()+1; // 'grid-style' indexing is 1-based
    i1 = domain_->LastI()+1;
    j0 = domain_->FirstJ()+1;
    j1 = domain_->LastJ()+1;
    k0 = domain_->FirstK()+1;
    k1 = domain_->LastK()+1;

    //add the boundary cells i=0,n_+1 etc (this is independent of overlap)
    i0--; i1++; j0--; j1++; k0--; k1++;

    DEBUG("create landmap with overlap...");
    Teuchos::RCP<Epetra_Map> landmap_loc = Utils::CreateMap(i0,i1,j0,j1,k0,k1,
                                                            I0,I1,J0,J1,K0,K1,*comm_);

    DEBUG("Create local vectors...");

    // distributed non-overlapping version of landm
    Teuchos::RCP<Epetra_IntVector> landm_loc0 =
        Teuchos::rcp(new Epetra_IntVector(*landmap_loc0));

    // distributed overlapping version of landm
    Teuchos::RCP<Epetra_IntVector> landm_loc =
        Teuchos::rcp(new Epetra_IntVector(*landmap_loc));

    DEBUG("Create importers...");

    // scatter to non-overlapping vector
    Teuchos::RCP<Epetra_Import> scatter, exchange;
    const Epetra_BlockMap& landmap_glb = landm_glb->Map();
    DEBUG("Create Gather-Import");
    scatter = Teuchos::rcp(new Epetra_Import(landmap_glb,*landmap_loc0));
    // exchange overlap
    DEBUG("Create Overlap-Import");
    exchange = Teuchos::rcp(new Epetra_Import(*landmap_loc,*landmap_loc0));


    // this is a 'scatter' operation to an overlapping distribution
    DEBUG("Distribute landmask...");
    // this helps to identify errors
    landm_loc0->PutValue(-999);
    landm_loc->PutValue(42);

    // scatter
    CHECK_ZERO(landm_loc0->Export(*landm_glb, *scatter,Insert));

    // get boundaries correct
    CHECK_ZERO(landm_loc->Import(*landm_loc0, *exchange,Insert));

    return landm_loc;
}

//=============================================================================
void THCM::setIntCondCorrection(Teuchos::RCP<Epetra_Vector> vec)
{
    if (sres_ != 0)
    {
        WARNING("This should not be called when SRES!=0",
                __FILE__, __LINE__);
        return;
    }
    else
    {
        // compute salinity integral, put it in correction
        CHECK_ZERO(intcondCoeff_->Dot(*vec, &intCorrection_));

        if (std::abs(intCorrection_) > 1e-8)
        {
            INFO("THCM integral correction: " << intCorrection_);
        }
    }
}

//=============================================================================
//! Let THCM perform integral checks
void THCM::integralChecks(Teuchos::RCP<Epetra_Vector> state,
                          double &salt_advection,
                          double &salt_diffusion)
{
    if (!(state->Map().SameAs(*standardMap_)))
    {
        ERROR("Map of input vector not same as standard map ",__FILE__,__LINE__);
    }

    // Create vectors for integral coefficients
    Teuchos::RCP<Epetra_Vector> globalCoeff =
        Teuchos::rcp(new Epetra_Vector(*standardVolumeMap_));
    Teuchos::RCP<Epetra_Vector> localCoeff =
        Teuchos::rcp(new Epetra_Vector(*assemblyVolumeMap_));

    // Create pointer to view of local coefficients
    double *localCoeffView;
    localCoeff->ExtractView(&localCoeffView);

    // Create local state
    Teuchos::RCP<Epetra_Vector> localState =
        Teuchos::rcp(new Epetra_Vector(*assemblyMap_));

    // Import into local state
    domain_->Solve2Assembly(*state, *localState);

    // Create pointer to view of local state entries
    double *localStateView;
    localState->ExtractView(&localStateView);

    // --------------------------------------------
    // Compute salt advection volume integral
    F90NAME( m_integrals, salt_advection )( localStateView, localCoeffView );

    // Export local coefficients into global entries
    CHECK_ZERO(globalCoeff->Export(*localCoeff, *as2std_vol_, Zero));

    // Compute integral on global non-overlapping domain
    double localInt = 0.0;
    for (int i = 0; i != globalCoeff->Map().NumMyElements(); ++i)
        localInt += (*globalCoeff)[i];

    // Sum over subdomains, obtain resulting integral
    comm_->SumAll(&localInt, &salt_advection, 1);

    std::cout << "Salt advection integral, PID = " << comm_->MyPID() << " lSum = " << localInt
              << " gSum = " << salt_advection << std::endl;

    // Reset local coefficients and integral
    localCoeff->PutScalar(0.0);
    localInt = 0.0;

    // --------------------------------------------
    // Compute salt diffusion volume integral
    F90NAME( m_integrals, salt_diffusion )( localStateView, localCoeffView );

    // Export local coefficients into global non overlapping entries
    CHECK_ZERO(globalCoeff->Export(*localCoeff, *as2std_vol_, Zero));

    // std::ofstream file;
    // std::stringstream ss;
    // ss << "globalCoeff" << comm_->MyPID();
    // file.open(ss.str());
    // globalCoeff->Print(file);
    // file.close();

    // Compute integral on subdomains
    for (int i = 0; i != globalCoeff->Map().NumMyElements(); ++i)
        localInt += (*globalCoeff)[i];

    // Sum over subdomains
    comm_->SumAll(&localInt, &salt_diffusion, 1);

    std::cout << "Salt diffusion integral, PID = " << comm_->MyPID() << " lSum = " << localInt
              << " gSum = " << salt_diffusion << std::endl;
}

//=============================================================================
// implement integral condition for S in Jacobian and B-matrix
void THCM::intcond_S(Epetra_CrsMatrix& A, Epetra_Vector& B)
{
    int N=domain_->GlobalN();
    int M=domain_->GlobalM();
    int L=domain_->GlobalL();

    int root = comm_->NumProc()-1;

    Teuchos::RCP<Epetra_MultiVector> intcond_glob =
        Utils::Gather(*intcondCoeff_, root);

    if (A.MyGRID(rowintcon_))
    {
        if (comm_->MyPID()!=root)
        {
            ERROR("S-integral condition should be on last processor!",__FILE__,__LINE__);
        }
        int lid = B.Map().LID(rowintcon_);
        B[lid]  = 0.0;   // no more time-dependence for this S-point
        int len = N*M*L;

        double *values = new double[len];
        int *indices   = new int[len];

        int pos=0;
        for (int i=0;i<N;i++)
            for (int j=0;j<M;j++)
                for (int k=0;k<L;k++)
                {
                    int gid = FIND_ROW2(_NUN_,N,M,L,i,j,k,SS);
                    indices[pos] = gid;

                    values[pos] = intSign_ * (*intcond_glob)[0][gid];
                    pos++;
                }

        /*
          len=1;
          indices[0]=rowintcon_;
          values[0]=1.0;
        */
        int ierr;
        if (A.Filled())
        {
            ierr = A.ReplaceGlobalValues(rowintcon_,len,values,indices);
        }
        else
        {
            ierr = A.InsertGlobalValues(rowintcon_,len,values,indices);
        }
        if (ierr != 0)
        {
            std::stringstream ss;
            ss << "graph_pid" << comm_->MyPID();
            std::ofstream file(ss.str());
            file << A.Graph();

            std::cout << "Insertion ERROR! " << ierr << ", filled = "
                      << A.Filled() << std::endl;

            std::cout << "\n ERROR " << ierr;
            std::cout << ((ierr == 2) ? ": value excluded" : "") << std::endl;
            std::cout << "\n myPID " << comm_->MyPID();

            std::cout << " while inserting/replacing values in local Jacobian" << std::endl;
            std::cout << "  GRID: " << rowintcon_ << std::endl;
            ERROR("Error during insertion/replacing of values in local Jacobian",
                  __FILE__, __LINE__);
        }

        delete []  values;
        delete []  indices;
    }
    else if (comm_->MyPID()==root)
    {
        ERROR("S-integral condition should be on last processor!",__FILE__,__LINE__);
    }
}

//=============================================================================
void THCM::fixPressurePoints(Epetra_CrsMatrix& A, Epetra_Vector& B)
{
    for (int i=1;i<=2;i++) {
        int row = (i==1)? rowPfix1_: rowPfix2_;
        if (A.MyGRID(row))
        {
            int lidB = B.Map().LID(row);
            int lidA = A.RowMap().LID(row);
            B[lidB] = 0.0; // no more time-dependence for this P-point

            int numEntries = A.NumMyEntries(lidA);
            double *vals   = new double[numEntries];
            int *inds      = new int[numEntries];

            // Extract current row and zero out except diagonal
            CHECK_NONNEG(A.ExtractGlobalRowCopy(row, numEntries, numEntries, vals, inds));
            for (int i = 0; i != numEntries; ++i)
            {
                if (inds[i] == row)
                    vals[i] = 1.0;
                else
                    vals[i] = 0.0;
            }

            if (A.Filled())
            {
                CHECK_NONNEG(A.ReplaceGlobalValues(row,numEntries,vals,inds));
            }
            else
            {
                CHECK_NONNEG(A.InsertGlobalValues(row,numEntries,vals,inds));
            }

            delete [] vals;
            delete [] inds;
        }
    }//for two pressure dirichlet values
}

//=============================================================================
Teuchos::RCP<Epetra_CrsGraph> THCM::CreateMaximalGraph(bool useSRES)
{
    DEBUG("Constructing maximal matrix graph...");
    int n=domain_->LocalN();
    int m=domain_->LocalM();
    int l=domain_->LocalL();
    int ndim = standardMap_->NumMyElements();
    int *numEntriesPerRow = new int[ndim];
    //int *landm = new int[n*m*l];
    //F90NAME(m_thcm_utils,get_landm)(landm);
    //const int LAND=1;
    for (int k=1;k<=l;k++)
        for (int j=1;j<=m;j++)
            for (int i=1;i<=n;i++)
            {
                int lidU = FIND_ROW2(_NUN_,n,m,l,i-1,j-1,k-1,UU);
                int gidU = assemblyMap_->GID(lidU);
                if (standardMap_->MyGID(gidU)) // otherwise: ghost cell, not in Jacobian
                {
                    int lid0 = standardMap_->LID(gidU)-1;
                    numEntriesPerRow[lid0+UU] = 24;
                    numEntriesPerRow[lid0+VV] = 22;
                    numEntriesPerRow[lid0+WW] = 7;
                    numEntriesPerRow[lid0+PP] = 11;
                    numEntriesPerRow[lid0+TT] = 20;
                    numEntriesPerRow[lid0+SS] = 20;
                }
            }

    Teuchos::RCP<Epetra_CrsGraph> graph
        = Teuchos::rcp(new Epetra_CrsGraph(Copy,*standardMap_,numEntriesPerRow,false));

    DEBVAR(rowintcon_);
    DEBVAR(sres_);

    int indices[24];
    int N = domain_->GlobalN();
    int M = domain_->GlobalM();
    int L = domain_->GlobalL();

    int I0 = domain_->FirstRealI();
    int J0 = domain_->FirstRealJ();
    int K0 = domain_->FirstRealK();
    int I1 = domain_->LastRealI();
    int J1 = domain_->LastRealJ();
    int K1 = domain_->LastRealK();
    int pos; // counts nonzero's per row and keeps track of position
    for (int k=K0; k<=K1; k++)
        for (int j=J0; j<=J1; j++)
            for (int i=I0; i<=I1; i++)
            {
                int gidU = FIND_ROW2(_NUN_,N,M,L,i,j,k,UU);
                int gid0 = gidU-1;

                // U-equation
                pos=0;

                // u-u: 7-point stencil
                insert_graph_entry(indices,pos,i,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,UU,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,UU,N,M,L);

                // u-v: 5-point stencil TODO: why is it 5 in U and 3 in V-eqn?
                insert_graph_entry(indices,pos,i,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,VV,N,M,L);

                // u-w: cell-averages at level k and k-1
                insert_graph_entry(indices,pos,i,j,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i+1,j+1,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,WW,N,M,L);

                insert_graph_entry(indices,pos,i,j,k-1,WW,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k-1,WW,N,M,L);
                insert_graph_entry(indices,pos,i+1,j+1,k-1,WW,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k-1,WW,N,M,L);

                // u-p: gradient
                insert_graph_entry(indices,pos,i,j,k,PP,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,PP,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,PP,N,M,L);
                insert_graph_entry(indices,pos,i+1,j+1,k,PP,N,M,L);

//#ifdef DEBUGGING
#if 0
#define DEBUG_GRAPH_ROW(var)                                            \
                std::cout << "graph row "<<i<<" "<<j<<" "<<k<<" "<<var;     \
                std::cout << " (gid "<<(gid0+var)<<"):"<<std::endl;         \
                std::cout << "predicted length: "<<numEntriesPerRow[standardMap_->LID(gid0+var)]; \
                std::cout <<", actual length: "<<pos<<std::endl;        \
                for (int pp=0;pp<pos;pp++) std::cout << (indices)[pp]<<" ";     \
                std::cout << std::endl;
#else
#define DEBUG_GRAPH_ROW(var)
#endif
//"
                DEBUG_GRAPH_ROW(UU)
                    CHECK_ZERO(graph->InsertGlobalIndices(gid0+UU,pos,indices));

                // v-equation
                pos=0;

                // v-v: 7-point stencil
                insert_graph_entry(indices,pos,i,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,VV,N,M,L);

                // v-u: 3-point stencil (why not 5, as in u-v?)
                insert_graph_entry(indices,pos,i,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,UU,N,M,L);

                // v-w: cell-averages at level k and k-1
                insert_graph_entry(indices,pos,i,j,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i+1,j+1,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,WW,N,M,L);

                insert_graph_entry(indices,pos,i,j,k-1,WW,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k-1,WW,N,M,L);
                insert_graph_entry(indices,pos,i+1,j+1,k-1,WW,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k-1,WW,N,M,L);

                // v-p: gradient
                insert_graph_entry(indices,pos,i,j,k,PP,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,PP,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,PP,N,M,L);
                insert_graph_entry(indices,pos,i+1,j+1,k,PP,N,M,L);

                DEBUG_GRAPH_ROW(VV)
                    CHECK_ZERO(graph->InsertGlobalIndices(gid0+VV,pos,indices));

                // w-equation
                pos=0;

                insert_graph_entry(indices,pos,i,j,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i,j,k,PP,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,PP,N,M,L);
                insert_graph_entry(indices,pos,i,j,k,TT,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,TT,N,M,L);
                insert_graph_entry(indices,pos,i,j,k,SS,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,SS,N,M,L);

                DEBUG_GRAPH_ROW(WW)
                    CHECK_ZERO(graph->InsertGlobalIndices(gid0+WW,pos,indices));

                // continuity-equation
                pos=0;

                insert_graph_entry(indices,pos,i,j,k,PP,N,M,L);

                insert_graph_entry(indices,pos,i,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i-1,j-1,k,UU,N,M,L);

                insert_graph_entry(indices,pos,i,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i-1,j-1,k,VV,N,M,L);

                insert_graph_entry(indices,pos,i,j,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,WW,N,M,L);

                DEBUG_GRAPH_ROW(PP)
                    CHECK_ZERO(graph->InsertGlobalIndices(gid0+PP,pos,indices));

                // T-equation
                pos=0;

                // T-T: 7-point stencil
                insert_graph_entry(indices,pos,i,j,k,TT,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,TT,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,TT,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,TT,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,TT,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,TT,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,TT,N,M,L);

                // T-U/V/W:
                insert_graph_entry(indices,pos,i,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i-1,j-1,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,UU,N,M,L);

                insert_graph_entry(indices,pos,i,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i-1,j-1,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,VV,N,M,L);

                insert_graph_entry(indices,pos,i,j,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,WW,N,M,L);

                // T-S: these terms may enter due to convective adjustment
                insert_graph_entry(indices,pos,i,j,k,SS,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,SS,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,SS,N,M,L);

                DEBUG_GRAPH_ROW(TT)
                    CHECK_ZERO(graph->InsertGlobalIndices(gid0+TT,pos,indices));

                // S-equation
                pos=0;

                // S-S: 7-point stencil
                insert_graph_entry(indices,pos,i,j,k,SS,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,SS,N,M,L);
                insert_graph_entry(indices,pos,i+1,j,k,SS,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,SS,N,M,L);
                insert_graph_entry(indices,pos,i,j+1,k,SS,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,SS,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,SS,N,M,L);

                // S-U/V/W:
                insert_graph_entry(indices,pos,i,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i-1,j-1,k,UU,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,UU,N,M,L);

                insert_graph_entry(indices,pos,i,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i-1,j,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i-1,j-1,k,VV,N,M,L);
                insert_graph_entry(indices,pos,i,j-1,k,VV,N,M,L);

                insert_graph_entry(indices,pos,i,j,k,WW,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,WW,N,M,L);

                // S-T: these terms may enter due to convective adjustment
                insert_graph_entry(indices,pos,i,j,k,TT,N,M,L);
                insert_graph_entry(indices,pos,i,j,k-1,TT,N,M,L);
                insert_graph_entry(indices,pos,i,j,k+1,TT,N,M,L);
#ifndef NO_INTCOND
                if ( (sres_ == 0) &&
                     ( (gid0+SS) == rowintcon_ ) &&
                     useSRES)
                    continue;
#endif
                DEBUG_GRAPH_ROW(SS)
                    CHECK_ZERO(graph->InsertGlobalIndices(gid0+SS,pos,indices));
            }

#ifndef NO_INTCOND
    if ((sres_ == 0) && useSRES)
    {
        int grid = rowintcon_;
        if (standardMap_->MyGID(grid))
        {
            int len = N*M*L;
            int *inds = new int[len];
            pos=0;
            for (int i=0;i<N;i++)
                for (int j=0;j<M;j++)
                    for (int k=0;k<L;k++)
                    {
                        int gcid = FIND_ROW2(_NUN_,N,M,L,i,j,k,SS);
                        inds[pos] = gcid;
                        pos++;
                    }
            CHECK_NONNEG(graph->InsertGlobalIndices(grid,len,inds));
            delete [] inds;
        }
    }
#endif
    CHECK_ZERO(graph->FillComplete());

    //delete [] landm;
    delete [] numEntriesPerRow;

    return graph;
}

//=============================================================================
void THCM::insert_graph_entry(int* indices, int& pos,
                              int i, int j, int k, int xx,
                              int N, int M, int L) const
{
    int ii=i; // if x-boundary is periodic i may be out of bounds.
    // ii will be adjusted in that case:
    if (domain_->IsPeriodic())
    {
        ii = MOD((double)i,(double)N);
    }
    if ((ii>=0) && (j>=0) && (k>=0) &&
        (ii< N) && (j< M) && (k< L) )
    {
        indices[pos++] = FIND_ROW2(_NUN_,N,M,L,ii,j,k,xx);
    }
}

//=============================================================================
double THCM::getSCorr()
{
    getFluxes();
    return scorr_;
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getIntCondCoeff()
{
    intcondCoeff_->PutScalar(0.0);
    Teuchos::RCP<Epetra_Vector> intcond_tmp =
        Teuchos::rcp(new Epetra_Vector(*assemblyMap_));

    int nml = (domain_->LocalN())*(domain_->LocalM())*(domain_->LocalL());

    double *values = new double[nml];
    int *indices   = new int[nml];
    int len;
    F90NAME(m_thcm_utils, intcond_scaling)(values, indices, &len);

    for (int i=0;i<len;i++)
    {
        (*intcond_tmp)[indices[i]-1] = values[i];
    }

    delete [] values;
    delete [] indices;

    domain_->Assembly2Solve(*intcond_tmp,*intcondCoeff_);

    intcondCoeff_->Norm1(&totalVolume_);

    INFO("  total volume: " << totalVolume_);

    return intcondCoeff_;
}

//=============================================================================
// set vmix_fix
void THCM::fixMixing(int value)
{
    if (vmix_ == 2)
    {
        INFO(" ** fixing vmix_fix: " << value << " **");
        F90NAME(m_mix, set_vmix_fix)(&value);
    }
}

//=============================================================================
extern "C" {

// this is a cheat for the fortran routine fsint from forcing.F90
    void thcm_forcing_integral_(double* qfun2, double* y, int* landm, double* fsint)
    {
        Teuchos::RCP<Epetra_Comm> comm = THCM::Instance().GetComm();
        Teuchos::RCP<TRIOS::Domain> domain_ = THCM::Instance().GetDomain();

        int n = domain_->LocalN();
        int m = domain_->LocalM();
        int l = domain_->LocalL();

        double lsint = 0.0, lfsint=0.0, sint;

        int i0 = domain_->FirstRealI()-domain_->FirstI();
        int j0 = domain_->FirstRealJ()-domain_->FirstJ();

        int i1 = domain_->LastRealI()-domain_->FirstI();
        int j1 = domain_->LastRealJ()-domain_->FirstJ();

        for (int j=j0; j<=j1; j++)
        {
            for (int i=i0;i<=i1;i++)
            {
                //note: the fortran landm array is 0-based, so we add a 1 to i,j,k
                int pl = FIND_ROW2(1,n+2,m+2,l+2,i+1,j+1,l,1);
                int pq = FIND_ROW2(1,n,m,1,i,j,0,1);
                lfsint = qfun2[pq] * cos(y[j]) * (1 - landm[pl]) + lfsint;
                lsint  = cos(y[j]) * (1 - landm[pl]) + lsint;
            }
        }

        CHECK_ZERO( comm->SumAll(&lfsint,fsint,1) );
        CHECK_ZERO( comm->SumAll(&lsint,&sint,1) );

        *fsint = *fsint/sint;
    }

    //------------------------------------------------------------------
    // helper function to let the fortran code participate in error handling
    void thcm_throw_error_(char *msg)
    {
        ERROR(msg, "the fortran code", "somewhere");
    }
}

//=============================================================================
Teuchos::ParameterList
THCM::getDefaultInitParameters()
{
    Teuchos::ParameterList result = getDefaultParameters();
    result.setName("THCM Default Init Parameters");

    result.set("Problem Description", "Unnamed",
               "A descriptive name to identify the settings");

    result.set("Global Grid-Size n", 16,
               "Global number of grid point in the x-direction (west-east).");
    result.set("Global Grid-Size m", 16,
               "Global number of grid point in the y-direction (south-north).");
    result.set("Global Grid-Size l", 16,
               "Global number of grid point in the z-direction (depth).");

    Teuchos::RCP<Teuchos::EnhancedNumberValidator<double> > longitude_validator(
        new Teuchos::EnhancedNumberValidator<double>(-360.0, 360.0));

    Teuchos::RCP<Teuchos::EnhancedNumberValidator<double> > latitude_validator(
        new Teuchos::EnhancedNumberValidator<double>(-90.0, 90.0));

    // default: north atlantic
    result.set("Global Bound xmin", 286.0, "Western domain bound", longitude_validator);
    result.set("Global Bound xmax", 350.0, "Eastern domain bound", longitude_validator);
    result.set("Global Bound ymin", 10.0, "Southern domain bound", latitude_validator);
    result.set("Global Bound ymax", 74.0, "Northern domain bound", latitude_validator);
    result.set("Periodic", false, "Periodic boundary conditions in the x-direction");

    result.get("Depth hdim", 4000.0);
    result.get("Grid Stretching qz", 1.0);
    result.get("Topography", 1);
    result.get("Flat Bottom", false);
    result.get("Compute salinity integral", true);

    result.get("Read Land Mask", false); //== false in experiment0
    result.get("Land Mask","no_mask_specified");

    result.get("Inhomogeneous Mixing", 0);
    result.get("Mixing", 1);
    result.get("Rho Mixing", true);
    result.get("Taper", 1);
    result.get("Linear EOS: alpha T", 1.0e-4);
    result.get("Linear EOS: alpha S", 7.6e-4);
    result.get("Restoring Temperature Profile", 1);
    result.get("Restoring Salinity Profile", 1);
    result.get("Local SRES Only", false);
    result.get("Salinity Integral Sign", -1);
    result.get("Levitus T", 1);
    result.get("Levitus S", 1);
    result.get("Levitus Internal T/S", false);
    result.get("Coupled Temperature", 0);
    result.get("Coupled Salinity", 0);
    result.get("Coupled Sea Ice Mask", 1);
    result.get("Fix Pressure Points", false);
    result.get("Coriolis Force", 1);
    result.get("Forcing Type", 0);

    result.get("Read Salinity Perturbation Mask", false);
    result.get("Salinity Perturbation Mask", "no_mask_specified");

    result.get("Wind Forcing Type", 2);

    result.get("Wind Forcing Data", "wind/trtau.dat");
    result.get("Temperature Forcing Data", "levitus/new/t00an1");
    result.get("Salinity Forcing Data", "levitus/new/s00an1");

    result.get("Integral row coordinate i", -1);
    result.get("Integral row coordinate j", -1);

    result.get("Scaling","THCM");

    return result;
}

Teuchos::ParameterList
THCM::getDefaultParameters()
{
    Teuchos::ParameterList result("THCM Default Parameters");

    Teuchos::ParameterList& startParams = result.sublist("Starting Parameters");
    // Start from 1 since 0 (Time) isn't settable
    for (int i=1; i<= _NPAR_; i++)
    {
        std::string label = int2par(i);
        startParams.get(label, std::numeric_limits<double>::quiet_NaN());
    }

    return result;
}

const Teuchos::ParameterList& THCM::getParameters()
{ return paramList_; }

void THCM::setParameters(Teuchos::ParameterList& newParams)
{
    newParams.validateParameters(getDefaultParameters());
    paramList_.setParameters(newParams);

    for (auto& pair : paramList_.sublist("Starting Parameters")) {
        double val = pair.second.getValue(&val);
        if (!std::isnan(val)) setParameter(pair.first, val);
    }
}

Teuchos::RCP<const Epetra_MultiVector> THCM::getNullSpace()
{
    if (nullSpace_==Teuchos::null)
    {
        nullSpace_ = Teuchos::rcp(new Epetra_MultiVector
                                 (*standardMap_,2,true) );

        // the svp's are fairly easy to construct, they are
        // so-called 'checkerboard' modes' in the x-y planes.
        // we first construct them for the standard rectan-
        // gular subdomains and then export them to the load-
        // balanced 'solve' map TODO: typically the two are the
        // same (unless load balancing is active), so we currently
        // don't actually do the export. Load balancing is likely
        // to be discarded in the near future anyway.

        // loop over all non-ghost subdomain cells:
        int pos=PP-1;
        for (int k=domain_->FirstRealK();k<=domain_->LastRealK();k++)
            for (int j=domain_->FirstRealJ();j<=domain_->LastRealJ();j++)
                for (int i=domain_->FirstRealI();i<=domain_->LastRealI();i++)
                {
                    if ((i+j)%2)
                    {
                        (*(*nullSpace_)(0))[pos] = 1;
                    }
                    else
                    {
                        (*(*nullSpace_)(1))[pos] = 1;
                    }
                    pos+=_NUN_;
                }

        double nrm1,nrm2;
        CHECK_ZERO((*nullSpace_)(0)->Norm2(&nrm1));
        CHECK_ZERO((*nullSpace_)(1)->Norm2(&nrm2));
        CHECK_ZERO((*nullSpace_)(0)->Scale(1.0/nrm1));
        CHECK_ZERO((*nullSpace_)(1)->Scale(1.0/nrm2));
    }
    return nullSpace_;
}
