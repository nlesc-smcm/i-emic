/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

// for I-EMIC couplings
#include <math.h>
#include <iostream>
#include <memory>
#include <vector>

#include "Teuchos_StandardCatchMacros.hpp"

#include "Epetra_Comm.h"
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


//=============================================================================
extern "C" {

	// usrc.F90
	_SUBROUTINE_(setparcs)(int*,double*);
	_SUBROUTINE_(getparcs)(int*,double*);
	_SUBROUTINE_(writeparams)();
	_SUBROUTINE_(rhs)(double*,double*);
	_SUBROUTINE_(matrix)(double*,double*,double*);

	// input:   n,m,l,nmlglob
	//          xmin,xmax,ymin,ymax,
	//          periodic,landm,
	//          taux,tauy,tatm,emip,spert
	_SUBROUTINE_(init)(int*,int*,int*,int*,
                       double*,double*,double*,double*,
                       int*,int*,
                       double*,double*,double*,double*,double*);
	_SUBROUTINE_(finalize)(void);

	// global.F90
    // input:   N,M,L,
	//          Xmin,Xmax,Ymin,Ymax,hdim,qz,
	//          alphaT,alphaS,
	//          ih,vmix_GLB,tap,rho_mixing,
	//          periodic,itopo,flat,rd_mask,
	//          TRES,SRES,iza,ite,its,rd_spertm
	//          coupled_atm
	_MODULE_SUBROUTINE_(m_global,initialize)(int*,int*,int*,
                                             double*,double*,double*,double*,double*,double*,
                                             double*,double*,
                                             int*,int*,int*,int*,
                                             int*,int*,int*,int*,
                                             int*,int*,int*,int*,int*,int*,
		                                     int*);
	
	_MODULE_SUBROUTINE_(m_global,finalize)(void);
	_MODULE_SUBROUTINE_(m_global,get_landm)(int*);
	_MODULE_SUBROUTINE_(m_global,get_monthly_forcing)(double* tatm, double* emip,
													  double* taux, double* tauy, int* month);
	_MODULE_SUBROUTINE_(m_global,get_monthly_internal_forcing)(double* temp, double* salt,
															   int* month);
	_MODULE_SUBROUTINE_(m_global,get_windfield)(double* taux, double* tauy);
	_MODULE_SUBROUTINE_(m_global,get_temforcing)(double* tatm);
	_MODULE_SUBROUTINE_(m_global,get_salforcing)(double* emip);
	_MODULE_SUBROUTINE_(m_global,get_internal_temforcing)(double* temp);
	_MODULE_SUBROUTINE_(m_global,get_internal_salforcing)(double* salt);
	_MODULE_SUBROUTINE_(m_global,get_spert)(double* spert);
	_MODULE_SUBROUTINE_(m_monthly,set_forcing)(double*tatm, double* emip, double* taux,
											   double* tauy, int* month);
	_MODULE_SUBROUTINE_(m_monthly,set_internal_forcing)(double*temp, double* salt, int* month);
	_MODULE_SUBROUTINE_(m_usr,set_internal_forcing)(double*temp, double* salt);
	_MODULE_SUBROUTINE_(m_thcm_utils,get_landm)(int*);
	_MODULE_SUBROUTINE_(m_scaling,average_block)(double *db);
	_MODULE_SUBROUTINE_(m_scaling,compute)(double *db, double *rowscales, double* colscales);

	_SUBROUTINE_(fillcolb)(void);

	// CRS matrix allocation (module m_mat)
	_MODULE_SUBROUTINE_(m_mat,get_array_sizes)(int* nrows, int* nnz);
	_MODULE_SUBROUTINE_(m_mat,set_pointers)(int* nrows, int* nnz, int* beg,
											int* jco,double* co,double* coB);

	// compute scaling factors for S-integral condition. Values is an n*m*l array
	_MODULE_SUBROUTINE_(m_thcm_utils,intcond_scaling)(double* values,int* indices,int* len);

	// compute weights for load balancing
	_MODULE_SUBROUTINE_(m_thcm_utils,loadbal_weights)(double* weights,double*,double*,double*);

	// sets the vmix_fix flag
	_MODULE_SUBROUTINE_(m_mix,set_vmix_fix)(int* vmix_fix);

	// for time-dependent forcing (gamma* is a continuation parameter for wind, T and S):
	_MODULE_SUBROUTINE_(m_monthly,update_forcing)(double* t,
												  double* gammaw,double* gammat, double* gammas);
	_MODULE_SUBROUTINE_(m_monthly,update_internal_forcing)(double* t,
														   double* gammat, double* gammas);

	//---------------------- I-EMIC couplings--------------------------------------
	// Extensions created for communication within the I-EMIC
	//
	_MODULE_SUBROUTINE_(m_inserts, insert_atmosphere)(double *atmos);
	//-----------------------------------------------------------------------------
	
	_SUBROUTINE_(write_levitus)(const char*);

}//extern


//=============================================================================
THCM::THCM(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Comm> comm) :
    Singleton<THCM>(Teuchos::rcp(this, false)),
    Comm(comm),
	nullSpace(Teuchos::null),
	paramList(params)
	
{
	DEBUG("### enter THCM::THCM ###");

	std::string probdesc = paramList.get("Problem Description","Unnamed");

	n = paramList.get("Global Grid-Size n", 64);
	m = paramList.get("Global Grid-Size m", 64);
	l = paramList.get("Global Grid-Size l", 16);

	//=================================================
	// TODO: implement atmosphere layer
	la = 0;
	//=================================================

	std::stringstream ss;
	ss << "THCM (" << probdesc <<", "
	   << n << "x" << m << "x" << l <<")";
	INFO(ss.str());
	this->SetLabel(ss.str().c_str());

	// default: north atlantic
	xmin = paramList.get("Global Bound xmin", 286.0) * PI_ / 180.0;
	xmax = paramList.get("Global Bound xmax", 350.0) * PI_ / 180.0;
	ymin = paramList.get("Global Bound ymin", 10.0)  * PI_ / 180.0;
	ymax = paramList.get("Global Bound ymax", 74.0)  * PI_ / 180.0;
	periodic = paramList.get("Periodic"  , false);
	hdim     = paramList.get("Depth hdim", 4000.0);
	double qz      = paramList.get("Grid Stretching qz", 1.0);
	int    itopo   = paramList.get("Topography", 1);
	bool   flat    = paramList.get("Flat Bottom", false);
	bool   rd_mask = paramList.get("Read Land Mask", false); //== false in experiment0

	if (rd_mask)
    {
		std::string mask_file = paramList.get("Land Mask","no_mask_specified");
		// we put the name of the desired mask in a file so it can be
		// obtained from there by the fortran code:
		if (Comm->MyPID() == 0)
		{
			std::ofstream mfs("mask_name.txt", std::ios::trunc);
			mfs << mask_file;
		}
    }

	ih               = paramList.get("Inhomogeneous Mixing",0);
	vmix_GLB         = paramList.get("Mixing",1);
	rho_mixing       = paramList.get("Rho Mixing",true);
	tap              = paramList.get("Taper",1);
	alphaT           = paramList.get("Linear EOS: alpha T",1.0e-4);
	alphaS           = paramList.get("Linear EOS: alpha S",7.6e-4);
	tres             = paramList.get("Restoring Temperature Profile",1);
	sres             = paramList.get("Restoring Salinity Profile",1);
	ite              = paramList.get("Levitus T",2);
	its              = paramList.get("Levitus S",1);
	internal_forcing = paramList.get("Levitus Internal T/S",false);
	bool rd_spertm   = paramList.get("Read Salinity Perturbation Mask",false);
	coupled_atm      = paramList.get("Coupled Atmosphere", 0);

	if (rd_spertm)
    {
		std::string spertm_file = paramList.get("Salinity Perturbation Mask",
												"no_mask_specified");
		if (Comm->MyPID()==0)
		{
			std::ofstream mfs("spertm_name.txt",std::ios::trunc);
			mfs << spertm_file;
		}
    }

	iza  = paramList.get("Wind Forcing",1);

	int dof = _NUN_; // number of unknowns, defined in THCMdefs.H

	// construct an object to decompose the domain:
	domain = Teuchos::rcp(new TRIOS::Domain(n, m, l, dof, xmin, xmax, ymin, ymax,
											periodic, hdim, Comm));

	// perform a 2D decomposition of domain into rectangular boxes
	domain->Decomp2D();

	// get a map object representing the subdomain (with ghost-nodes/overlap).
	// This map defines the nodes local to the THCM subdomain
	AssemblyMap = domain->GetAssemblyMap();

	// get a map object representing the subdomain (without ghost-nodes/overlap).
	// this is an intermediate representation between the assembly and the solve
	// phases.
	StandardMap = domain->GetStandardMap();

	// initialize THCM (allocate memory etc.)
	// for a subdomain including ghost-nodes:

	// the domain object knows the geometry of the subdomain:
	double xminloc = domain->XminLoc();
	double xmaxloc = domain->XmaxLoc();
	double yminloc = domain->YminLoc();
	double ymaxloc = domain->YmaxLoc();
	// double zminloc = domain->ZminLoc();
	// double zmaxloc = domain->ZmaxLoc();

	DEBVAR(xmin);  //== Output variable into debug stream
	DEBVAR(xminloc);
	DEBVAR(xmax);
	DEBVAR(xmaxloc);
	DEBVAR(ymin);
	DEBVAR(yminloc);
	DEBVAR(ymax);
	DEBVAR(ymaxloc);
	DEBVAR(hdim);

	// and the number of grid points contained in it:
	int nloc = domain->LocalN();
	int mloc = domain->LocalM();
	int lloc = domain->LocalL();

	// global settings for THCM
	// memory for I/O is allocated only on the root process (pid 0)
	int nglob_ = 0;
	int mglob_ = 0;
	int lglob_ = 0;

	if (Comm->MyPID() == 0) // this one is responsible for I/O
    {
		nglob_ = n;
		mglob_ = m;
		lglob_ = l;
    }

	// fortran routine that puts global info in module m_global
	// and allocates some arrays there if dimensions are non-zero
	// this has to be called _before_ init because init already uses
	// some of it (by calling 'forcing')

	// for portability reasons we rather pass integers to fortran than bools
	// (this was an issue on Huygens)
	int irho_mixing = (rho_mixing) ? 1 : 0;
	int iperiodic   = (periodic  ) ? 1 : 0;
	int iflat       = (flat      ) ? 1 : 0;
	int ird_mask    = (rd_mask   ) ? 1 : 0;
	int ird_spertm  = (rd_spertm ) ? 1 : 0;

	DEBUG("call m_global::initialize...");
	//== In fortran object code this corresponds to the function __m_global_MOD_initialize
	F90NAME(m_global, initialize)(&nglob_, &mglob_, &lglob_,
								  &xmin, &xmax, &ymin, &ymax, &hdim, &qz,
								  &alphaT, &alphaS,
								  &ih, &vmix_GLB, &tap, &irho_mixing,
								  &iperiodic, &itopo, &iflat, &ird_mask,
								  &tres, &sres, &iza, &ite, &its, &ird_spertm,
								  &coupled_atm);

	// read topography data and convert it to a global land mask
	DEBUG("Initialize land mask...");

	int I0 = 0; int I1 = n+1;
	int J0 = 0; int J1 = m+1;
	int K0 = 0; int K1 = l+la+1;

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
	Teuchos::RCP<Epetra_Comm> xComm = domain->GetProcRow(0);

	int perio = (periodic && xComm->NumProc() == 1);

	int nmlglob = n*m*l;

	//--------------------------------------------------------------------------
	// read wind, temperature and salinity forcing and distribute it among
	// processors
	Teuchos::RCP<Epetra_Map> wind_map_loc   = domain->CreateAssemblyMap(1,true);
	Teuchos::RCP<Epetra_Map> lev_map_loc    = wind_map_loc;
	Teuchos::RCP<Epetra_Map> intlev_map_loc = domain->CreateAssemblyMap(1,false);
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

	Teuchos::RCP<Epetra_Map> wind_map_dist = domain->CreateStandardMap(1,true);
	Teuchos::RCP<Epetra_Map> wind_map_root = Utils::Gather(*wind_map_dist,0);

	Teuchos::RCP<Epetra_Vector> taux_glob     = Teuchos::rcp(new Epetra_Vector(*wind_map_root));
	Teuchos::RCP<Epetra_Vector> tauy_glob     = Teuchos::rcp(new Epetra_Vector(*wind_map_root));

	double *taux_g, *tauy_g;
	CHECK_ZERO(taux_glob->ExtractView(&taux_g));
	CHECK_ZERO(tauy_glob->ExtractView(&tauy_g));

	if (comm->MyPID()==0)
    {
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

	Teuchos::RCP<Epetra_Map> lev_map_dist  = domain->CreateStandardMap(1,true);
	Teuchos::RCP<Epetra_Map> lev_map_root  = Utils::Gather(*lev_map_dist,0);

	Teuchos::RCP<Epetra_Map> intlev_map_dist  = domain->CreateStandardMap(1,false);
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
		if (internal_forcing)
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
                &perio, landm,
                taux, tauy, tatm, emip, spert);

    if (internal_forcing)
	{
		F90NAME(m_usr,set_internal_forcing)(temp,salt);
	}

#ifdef DEBUGGING
	OceanGrid G(domain);
	G.ImportLandMask(*landm_loc);
	//DEBVAR(G);
#endif

#if 0
//#ifdef TESTING
	{
		std::stringstream ss;
		ss << "levitus_"<<Comm->MyPID()<<".txt";
		const char* fname = ss.str().c_str();
		INFO("store levitus fields in "<<fname);
		FNAME(write_levitus)(fname);
	}
#endif

	bool time_dep_forcing = paramList.get("Time Dependent Forcing",false);
	if (time_dep_forcing)
    {
		// read and distribute Levitus data
		this->SetupMonthlyForcing();
    }

	// get a map object for constructing vectors without overlap
	// (load-balanced, used for solve phase)
	SolveMap = domain->GetSolveMap();

	// Create internal vectors
	initialSolution = Teuchos::rcp(new Epetra_Vector(*SolveMap));
	diagB           = Teuchos::rcp(new Epetra_Vector(*SolveMap));
	localDiagB      = Teuchos::rcp(new Epetra_Vector(*StandardMap));
	localRhs        = Teuchos::rcp(new Epetra_Vector(*AssemblyMap));
	localSol        = Teuchos::rcp(new Epetra_Vector(*AssemblyMap));

	// allocate mem for the CSR matrix in THCM.

	int nrows, nnz;

	// first ask how big it should be:
	DEBUG("call get_array_sizes...");
	F90NAME(m_mat,get_array_sizes)(&nrows,&nnz);

	DEBUG("Allocating Fortran CSR arrays, nrows="<<nrows<<", nnz="<<nnz);

	// allocate the memory
	begA = new int[nrows+1];
	coA  = new double[nnz];
	jcoA = new int[nnz];
	coB  = new double[nrows];

	// give THCM the opportunity to set its pointers to
	// the new memory block
	DEBVAR("call set_pointers...");
	F90NAME(m_mat,set_pointers)(&nrows,&nnz,begA,jcoA,coA,coB);

	rowintcon_=-1;
#ifndef NO_INTCOND
	if (sres==0)
    {
		int N=domain->GlobalN();
		int M=domain->GlobalM();
		int L=domain->GlobalL();
		rowintcon_ = FIND_ROW2(_NUN_,N,M,L,N-1,M-1,L-1,SS);
		INFO("integral condition for S is in global row "<<rowintcon_);
		intcond_coeff = Teuchos::rcp(new Epetra_Vector(*SolveMap));
		intcond_coeff->PutScalar(0.0);
		Teuchos::RCP<Epetra_Vector> intcond_tmp = Teuchos::rcp(new Epetra_Vector(*AssemblyMap));
		int nml = (domain->LocalN())*(domain->LocalM())*(domain->LocalL());
		double *values = new double[nml];
		int *indices = new int[nml];
		int len;
		F90NAME(m_thcm_utils,intcond_scaling)(values,indices,&len);
		for (int i=0;i<len;i++)
		{
			(*intcond_tmp)[indices[i]-1] = values[i];
		}
		delete [] values;
		delete [] indices;
		domain->Assembly2Solve(*intcond_tmp,*intcond_coeff);
#ifdef STORE_MATRICES
		std::ofstream ofs("intcond.txt");
		ofs << *intcond_coeff;
		ofs.close();
#endif
    }
#endif
	// create a graph describing the maximal matrix pattern.
	// Note that in LOCA we can't change the pattern of the matrix
	// during the continuation process as we pass pointers to LOCA
	// which can't be changed, whereas changing the pattern would
	// require building a whole new matrix. As we can't predict
	// where convective adjustment will happen, we assume it hap-
	// pens everywhere.
	localMatrixGraph=this->CreateMaximalGraph();

	if (SolveMap!=StandardMap)
    {
		MatrixGraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*SolveMap,20));
		Teuchos::RCP<Epetra_Import> import = Teuchos::rcp(new Epetra_Import(*StandardMap,*SolveMap));
		DEBUG("Migrate graph to solve map...");
		CHECK_ZERO(MatrixGraph->Export(*localMatrixGraph,*import,Insert));
		CHECK_ZERO(MatrixGraph->FillComplete());
    }
	else
    {
		MatrixGraph = localMatrixGraph;
    }
	localJac = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *localMatrixGraph));
	localJac->SetLabel("Local Jacobian");
	Jac = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *MatrixGraph));
	Jac->SetLabel("Jacobian");

	// we do not allow these to be modified, use the OceanModel interface
	// (computeShiftedMatrix, setXdot)
	// for implementing Time integration or eigenvalue things.
	sigmaUVTS=0.0;
//  sigmaWP=1.0e-14;
	sigmaWP=0.0;

// we can select two points where the continuity equation will be replaced by
// P(i,j,k) = 0. This is experimental, we hope to fix the divergence problem in the 4D case
//  like this
// 	int N = domain->GlobalN();
// 	int M = domain->GlobalM();
// 	int L = domain->GlobalL();
//
// rowPfix1 = FIND_ROW2(_NUN_,N,M,L,N-1,M-1,L-1,PP);
// rowPfix2 = FIND_ROW2(_NUN_,N,M,L,N-2,M-1,L-1,PP);
	
	rowPfix1=-1;
	rowPfix2=-1;

	// build vectonr with integral coefficients
	this->evaluateB();

	scaling_type=paramList.get("Scaling","THCM");

	if (scaling_type=="THCM")
    {
		// construct the scaling object. The scaling is computed by THCM (m_scaling)
		// and passed on to Trilinos:
		row_scaling = Teuchos::rcp(new Epetra_Vector(*SolveMap));
		row_scaling->SetLabel("Row Scaling");
		col_scaling = Teuchos::rcp(new Epetra_Vector(*SolveMap));
		col_scaling->SetLabel("Col Scaling");
		local_row_scaling = Teuchos::rcp(new Epetra_Vector(*AssemblyMap) );
		local_col_scaling = Teuchos::rcp(new Epetra_Vector(*AssemblyMap) );

		row_scaling->PutScalar(1.0);
		col_scaling->PutScalar(1.0);
    }

}

//=============================================================================
THCM::~THCM()
{
	this->printTiming(std::cout);
	DEBUG("Destroy THCM...");

	FNAME(finalize)();
	if (Comm->MyPID()==0)
    {
		F90NAME(m_global,finalize)();
    }

	delete [] jcoA;
	delete [] coA;
	delete [] begA;
	delete [] coB;

	// the rest is handled by Teuchos::rcp's
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> THCM::getSolution()
{
	return initialSolution;
}

//=============================================================================
Teuchos::RCP<Epetra_CrsMatrix> THCM::getJacobian()
{
	return Jac;
}

//=============================================================================
// Compute Jacobian and/or RHS.
bool THCM::evaluate(const Epetra_Vector& soln,
				    Teuchos::RCP<Epetra_Vector> tmp_rhs,
				    bool computeJac)
{

#if defined(TESTING) && !defined(NO_INTCOND)
	if (sres==0)
	{
		double intcond;
		CHECK_ZERO(intcond_coeff->Dot(soln,&intcond));
		/*
		  double dx = (xmax-xmin)/domain->GlobalN();
		  double dy = (ymax-ymin)/domain->GlobalM();
		  double dz = 1.0/domain->GlobalL();
		  intcond = intcond*dx*dy*dz;
		*/
		INFO("Salinity integral condition (should be 0): " << intcond);
	}
#endif

	if (!(soln.Map().SameAs(*SolveMap)))
    {
		ERROR("Map of solution vector not same as solve-map ",__FILE__,__LINE__);
    }

	// convert to standard distribution and
	// import values from ghost-nodes on neighbouring subdomains:
	domain->Solve2Assembly(soln,*localSol);


	int NumMyElements = AssemblyMap->NumMyElements();

//  DEBUG("=== evaluate: input vector");
//  DEBUG( (domain->Gather(*soln,0)) )

	// extract an array to pass on to THCM:
	// We use 'Copy' mode, which is clean but possibly slow.
	// Probably 'View' would be allowable as well.
	double* solution;
	localSol->ExtractView(&solution);
	if(tmp_rhs!=Teuchos::null)
    {
		// INFO("Compute RHS...");
		// build rhs simultaneously on each process
		double* RHS;
		CHECK_ZERO(localRhs->ExtractView(&RHS));
		TIMER_START("Ocean: compute rhs: fortran part");
		// compute right-hand-side on whole subdomain (by THCM)
		FNAME(rhs)(solution, RHS); 
		TIMER_STOP("Ocean: compute rhs: fortran part");
		// export overlapping rhs to unique-id global rhs vector,
		// and load-balance for solve phase:
		domain->Assembly2Solve(*localRhs,*tmp_rhs);

		CHECK_ZERO(tmp_rhs->Scale(-1.0));

#ifndef NO_INTCOND
		if (sres==0)
		{
			int lastrow = rowintcon_;
			double intcond;
			//TODO: check which is better:
			CHECK_ZERO(intcond_coeff->Dot(soln,&intcond));
			//std::cout << " dot product: "<<intcond << std::endl;
			//intcond = 0.0;
			if (tmp_rhs->Map().MyGID(lastrow))
			{
				(*tmp_rhs)[tmp_rhs->Map().LID(lastrow)]=intcond;
			}
		}
#endif
		if (rowPfix1>=0)
		{
			if (tmp_rhs->Map().MyGID(rowPfix1))
			{
				(*tmp_rhs)[tmp_rhs->Map().LID(rowPfix1)]=0.0;
			}
		}
		if (rowPfix2>=0)
		{
			if (tmp_rhs->Map().MyGID(rowPfix2))
			{
				(*tmp_rhs)[tmp_rhs->Map().LID(rowPfix2)]=0.0;
			}
		}
    }
	if(computeJac)
    {
		// INFO("Compute Jacobian...");

		localJac->PutScalar(0.0); // set all matrix entries to zero
		localDiagB->PutScalar(0.0);

		if (sigmaUVTS||(sigmaWP>1.0e-14)) {
			ERROR("We do not allow THCM to shift the matrix anymore!",
				  __FILE__,__LINE__);
		}

		//Call the fortran routine, providing the solution vector,
		//and get back the three vectors of the sparse Jacobian (CSR form)
		TIMER_START("Ocean: compute jacobian: fortran part");
		FNAME(matrix)(solution,&sigmaUVTS,&sigmaWP);
		TIMER_STOP("Ocean: compute jacobian: fortran part");
		
		const int maxlen = _NUN_*_NP_+1;    //nun*np+1 is max nonzeros per row
		int indices[maxlen];
		double values[maxlen];

		int index, numentries;

		int imax = NumMyElements;
#ifndef NO_INTCOND
		if (sres==0)
		{
			if (Jac->MyGRID(rowintcon_))
			{
				imax--;
			}
		}
#endif
		for (int i = 0; i<imax; i++)
		{
			if (!domain->IsGhost(i))
			{
				index = begA[i]; // note that these arrays use 1-based indexing
				numentries = begA[i+1] - index;
				for (int j = 0; j <  numentries ; j++)
				{
					indices[j] = AssemblyMap->GID(jcoA[index-1+j] - 1);
					values[j] = coA[index - 1 + j];
				}

				int ierr = localJac->ReplaceGlobalValues(AssemblyMap->GID(i), numentries,
														 values, indices);

				if ((ierr!=0) && (ierr!=3))
				{
					(std::cout) << "\nERROR " << ierr;
					(std::cout) <<" while inserting/replacing values in local Jacobian" << std::endl;
					DEBUG(" ERROR while inserting/replacing values in local Jacobian");

					INFO("GRID: "<<AssemblyMap->GID(i))
						INFO("number of entries: "<<numentries);
					(std::cout) << "entries: ";
					for (int j=0;j<numentries;j++) (std::cout) << "("<<indices[j]<<" "<<values[j]<<") ";

					CHECK_ZERO(localJac->ExtractGlobalRowCopy(AssemblyMap->GID(i),maxlen,numentries,values,indices));
					INFO("\noriginal row: ");
					INFO("number of entries: "<<numentries);
					(std::cout) << "entries: ";
					for (int j=0;j<numentries;j++) (std::cout) << "("<<indices[j]<<" "<<values[j]<<") ";
					(std::cout) << std::endl;

					// ierr == 3 probably means not all row entries are replaced,
					// does not matter because we zeroed them.
				}

				// reconstruct the diagonal matrix B
				int lid = StandardMap->LID(AssemblyMap->GID(i));
				double mass_param = 1.0;
				this->getParameter("Mass", mass_param);
				(*localDiagB)[lid] = -coB[i] * mass_param;
			} //not a ghost?
		} //i-loop over rows

#ifndef NO_INTCOND
		if (sres == 0)
		{
			this->intcond_S(*localJac,*localDiagB);
		}
#endif
		//this->fixPressurePoints(*localJac,*localDiagB);
		CHECK_ZERO(localJac->FillComplete());

		// redistribute according to SolveMap (may be load-balanced)
		domain->Standard2Solve(*localDiagB, *diagB);
		domain->Standard2Solve(*localJac, *Jac);
		CHECK_ZERO(Jac->FillComplete());

		if (scaling_type == "THCM")
		{
			DEBUG(" THCM:  RecomputeScaling()");
			this->RecomputeScaling();

#if 0
			std::ofstream ofs1("row_scaling.txt");
			ofs1 << *row_scaling;
			ofs1.close();
			std::ofstream ofs2("col_scaling.txt");
			ofs2 << *col_scaling;
			ofs2.close();
#endif
		}


    } // matrix
	return true;
}


// just reconstruct the diagonal matrix B from THCM
void THCM::evaluateB(void)
{

	int NumMyElements = AssemblyMap->NumMyElements();

	DEBUG("Construct matrix B...");

	localDiagB->PutScalar(0.0);
	FNAME(fillcolb)();
	for (int i = 0; i<NumMyElements; i++)
    {
		if (!domain->IsGhost(i))
		{
			// reconstruct the diagonal matrix B
			int lid = StandardMap->LID(AssemblyMap->GID(i));
			(*localDiagB)[lid] = -coB[i];
		}//not a ghost?
    }//i-loop over rows
#ifndef NO_INTCOND
	if (sres==0)
    {
		int lastrow = (domain->GlobalN())*(domain->GlobalM())*(domain->GlobalL());
		if (localDiagB->Map().MyGID(lastrow))
		{
			(*localDiagB)[localDiagB->Map().LID(lastrow)]=0.0;
		}
    }
#endif
	domain->Standard2Solve(*localDiagB,*diagB);
}	
	
//=============================================================================
// I-EMIC stuff
//=============================================================================
std::shared_ptr<std::vector<int> > THCM::getLandMask()
{

	// length of landmask array
	int dim = (n+2)*(m+2)*(l+la+2);

	// Create landmask array
	std::shared_ptr<std::vector<int> > landm =
		std::make_shared<std::vector<int> >(dim, 0);
	
	// Let THCM fill the landmask array on proc = 0
	if (Comm->MyPID() == 0)
		F90NAME(m_global,get_landm)(&(*landm)[0]);
	
#ifdef HAVE_MPI
	// Get the MpiComm from Epetra
	Epetra_MpiComm const MpiComm =
		dynamic_cast<Epetra_MpiComm const &>(*Comm); 
	MPI_Bcast(&(*landm)[0], dim, MPI_INTEGER, 0, MpiComm.GetMpiComm());
#endif 
	return landm;
}
		
//=============================================================================
void THCM::setAtmosphere(std::vector<double> const &atmosvec)
{
	// Create a gather map
	Teuchos::RCP<Epetra_Map> atmos_map_dist = domain->CreateStandardMap(1, true);
	Teuchos::RCP<Epetra_Map> atmos_map_root = Utils::Gather(*atmos_map_dist, 0);

	// Copy atmosphere to non-const vector
	std::vector<double> atmosCpy(atmosvec);
	
	// Insert the atmosphere
	Teuchos::RCP<Epetra_Vector> atmos_glob =
		Teuchos::rcp(new Epetra_Vector(Copy, *atmos_map_root, &atmosCpy[0]));

	// Distribute the atmosphere
	Teuchos::RCP<Epetra_MultiVector> atmos_dist =
		Utils::Scatter(*atmos_glob, *atmos_map_dist);

	// Create assembly vector
	Teuchos::RCP<Epetra_Map> atmos_map_loc  = domain->CreateAssemblyMap(1, true);
	Teuchos::RCP<Epetra_Vector> atmos_loc =
		Teuchos::rcp(new Epetra_Vector(*atmos_map_loc));

	// Import overlap
	Teuchos::RCP<Epetra_Import> atmos_loc2dist =
		Teuchos::rcp(new Epetra_Import(*atmos_map_loc, *atmos_map_dist));

	// Insert the assembly into THCM
	atmos_loc->Import(*atmos_dist, *atmos_loc2dist, Insert);
	double *atmos;
	atmos_loc->ExtractView(&atmos);
	F90NAME(m_inserts, insert_atmosphere)(atmos);	
}	

//=============================================================================
void THCM::setAtmosphereTest()
{
	// This is a test.
	// We build an idealized atmosphere on the assembly map to give to THCM
	Teuchos::RCP<Epetra_Map> atmos_map_loc = domain->CreateAssemblyMap(1, true);
	Teuchos::RCP<Epetra_Vector> atmos_loc  =
		Teuchos::rcp(new Epetra_Vector(*atmos_map_loc));

	int nloc       =  domain->LocalN();
	int mloc       =  domain->LocalM();
	double yminLoc =  domain->YminLoc();
	double ymaxLoc =  domain->YmaxLoc();
	double ymax    =  domain->Ymax();
	double dyLoc   = (ymaxLoc - yminLoc) / (mloc - 1);
	double y       =  yminLoc;
	
	double value = 0;
	int idx      = 0;
	for (int j = 0; j != mloc; ++j)
	{
		y = yminLoc + j * dyLoc;
		for (int i = 0; i != nloc; ++i)
		{
			value = cos(PI_ * y / ymax);
			(*atmos_loc)[idx] = value;
			++idx;
		}
	}
	double *atmos_loc_array;
	atmos_loc->ExtractView(&atmos_loc_array);
	F90NAME(m_inserts, insert_atmosphere)(atmos_loc_array);
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

	int len = local_row_scaling->MyLength();
	double *rowscal = new double[len];
	double *colscal = new double[len];

	// compute local average diagonal block
	F90NAME(m_scaling,average_block)(ldb);

	Comm->SumAll(ldb,gdb,_NUN_*_NUN_);

	for (int i=0;i<_NUN_*_NUN_;i++)
    {
		gdb[i]/=Comm->NumProc();
    }
	// compute row- and column scaling
	F90NAME(m_scaling,compute)(gdb,rowscal,colscal);

	// put them back in assembly vectors:
	// note: the scaling matrices in Trilinos are
	// defined as the inverse of those in THCM
	for (int i=0;i<len;i++)
    {
		(*local_row_scaling)[i] = 1.0/rowscal[i];
		(*local_col_scaling)[i] = 1.0/colscal[i];
    }

	// kick out the ghost nodes:
	domain->Assembly2Solve(*local_row_scaling,*row_scaling);
	domain->Assembly2Solve(*local_col_scaling,*col_scaling);

	// make sure T and S are scaled the same way in each cell
	// we need this because of our special block scaling for the
	// ATS matrix in the preconditioner.
	for (int i=TT-1;i<row_scaling->MyLength();i+=_NUN_)
    {
		double mean = 0.5*((*row_scaling)[i]+(*row_scaling)[i+1]);
		(*row_scaling)[i] = mean;
		(*row_scaling)[i+1] = mean;
		mean = 0.5*((*col_scaling)[i]+(*col_scaling)[i+1]);
		(*col_scaling)[i] = mean;
		(*col_scaling)[i+1] = mean;
    }

	delete [] rowscal;
	delete [] colscal;

	// (2) diagonal row scaling for T and S (obsolete!)
	if (row_scaling_TS!=Teuchos::null)
    {
		CHECK_ZERO(Jac->ExtractDiagonalCopy(*row_scaling_TS));
		for (int i=0;i<row_scaling_TS->MyLength();i+=_NUN_)
		{
			for (int j=i;j<i+4;j++)
			{
				(*row_scaling_TS)[j] = 1.0; //u,v,w,p
			}
			for (int j=i+4;j<i+6;j++)
			{
				(*row_scaling_TS)[j] = 1.0/(*row_scaling_TS)[j]; //T,S
			}
		}
    }// additional T/S scaling (obsolete)
}

//=============================================================================
void THCM::normalizePressure(Epetra_Vector& soln) const
{

	int i = n/2-1;
	int j = 6*m/8-1;
	int k = l-1;
	int ref_gid = FIND_ROW2(_NUN_,n,m,l,i,j,k,PP);
	int ref_lid, ref_host;
	double ref_value;
	soln.Map().RemoteIDList(1, &ref_gid, &ref_host, &ref_lid);

	DEBUG("+++ NORMALIZE P +++\n");
	DEBVAR(ref_gid);
	DEBVAR(ref_lid);
	DEBVAR(ref_host);
	DEBUG("+++++++++++++++++++\n");

	if (ref_lid>=0) ref_value = soln[ref_lid];
	soln.Comm().Broadcast (&ref_value, 1, ref_host);

	//subtract reference value from all 'P' points except land cells
	// TODO: 1) we do not handle land cells correctly here, yet!
	//       2) the whole thing seems to go wrong...
    //  for (int i=PP;i<=soln.MyLength();i+=_NUN_) soln[i-1] -= ref_value;
}

//=============================================================================
// Timing functionality
void THCM::startTiming(std::string fname)
{
	Teuchos::RCP<Epetra_Time> T=Teuchos::rcp(new Epetra_Time(*Comm));
	timerList.sublist("timers").set(fname,T);
}


//=============================================================================
void THCM::stopTiming(std::string fname,bool print)
{
	Teuchos::RCP<Epetra_Time> T = Teuchos::null;
	T=timerList.sublist("timers").get(fname,T);
	double elapsed=0;
	if (T!=Teuchos::null)
    {
		elapsed=T->ElapsedTime();
    }
	int ncalls=timerList.sublist("number of calls").get(fname,0);
	double total_time=timerList.sublist("total time").get(fname,0.0);
	timerList.sublist("number of calls").set(fname,ncalls+1);
	timerList.sublist("total time").set(fname,total_time+elapsed);
	if (print)
    {
		(std::cout) << "### timing: "<<fname<<" "<<elapsed<<std::endl;
    }
}

//=============================================================================
void THCM::printTiming(std::ostream& os)
{
	os << "================= TIMING RESULTS ====================="<<std::endl;
	os << "     Description                              ";
	os << " # Calls \t Cumulative Time \t Time/call\n";
	os << "======================================================"<<std::endl;

	Teuchos::ParameterList& ncallsList=timerList.sublist("number of calls");
	Teuchos::ParameterList& elapsedList=timerList.sublist("total time");
	for (Teuchos::ParameterList::ConstIterator i=ncallsList.begin();i!=ncallsList.end();i++)
    {
		const std::string& fname = i->first;
		int ncalls = ncallsList.get(fname,0);
		double elapsed = elapsedList.get(fname,0.0);
		os << fname << "\t" <<ncalls<<"\t"<<elapsed<<"\t"
		   << ((ncalls>0)? elapsed/(double)ncalls : 0.0) <<std::endl;
    }
	os << "====================================================="<<std::endl;
	DEBUG(timerList);
}

//=============================================================================
// convert parameter name to integer
int THCM::par2int(std::string label)
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
	int MKAP   = 29; int SPL2   = 30; int EXPO   = 31; int SEAS   = 32;
	int SEASW  = 33; int SEAST  = 34; int SEASS  = 35; int MASS   = 36;

	if (label == "Time")                   return TIME;
	else if (label == "AL_T")                     return AL_T;
	else if (label == "Rayleigh-Number")          return RAYL;
	else if (label == "Rossby-Number")          return ROSB;
	else if (label == "Vertical Ekman-Number")    return EK_V;
	else if (label == "Horizontal Ekman-Number")  return EK_H;
	else if (label == "MIXP")                     return MIXP;
	else if (label == "RESC")                     return RESC;
	else if (label == "SPL1")                     return SPL1;
	else if (label == "Homotopy")                 return HMTP;
	else if (label == "Sun")                      return SUNP;
	else if (label == "Vertical Peclet-Number")   return PE_V;
	else if (label == "Horizontal Peclet-Number") return PE_H;
	else if (label == "P_VC")                     return P_VC;
	else if (label == "LAMB")                     return LAMB;
	else if (label == "ARCL")                     return ARCL;
	else if (label == "Salinity Forcing")         return SALT;
	else if (label == "Wind Forcing")             return WIND;
	else if (label == "Temperature Forcing")      return TEMP;
	else if (label == "Nonlinear Factor")         return BIOT;
	else if (label == "Combined Forcing")         return COMB;
	else if (label == "CONT")                     return CONT;
	else if (label == "IFRICB")                   return IFRICB;
	else if (label == "NLES")                     return NLES;
	else if (label == "CMPR")                     return CMPR;
	else if (label == "ALPC")                     return ALPC;
	else if (label == "Energy")                   return ENER;
	else if (label == "SPER")                     return SPER;
	else if (label == "FPER")                     return FPER;
	else if (label == "MKAP")                     return MKAP;
	else if (label == "SPL2")                     return SPL2;
	else if (label == "Exponent")                 return EXPO;
	else if (label == "Seasonal Forcing")         return SEAS;// combination of T,S and Wind
	else if (label == "Seasonal Forcing (Temperature)")  return SEAST;
	else if (label == "Seasonal Forcing (Salinity)")     return SEASS;
	else if (label == "Seasonal Forcing (Wind)")         return SEASW;
	else if (label == "Mass")                            return MASS;
	else
    {
		INFO("Invalid continuation parameter label: '"<<label<<"'");
		ERROR("Parameter label is invalid!",__FILE__,__LINE__);
    }
	return -1;
}

//=============================================================================
// convert parameter name to integer
std::string THCM::int2par(int index)
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
	int MKAP   = 29; int SPL2   = 30; int EXPO   = 31; int SEAS   = 32;
	int SEASW  = 33; int SEAST  = 34; int SEASS  = 35; int MASS   = 36;

	if (index==TIME) label = "Time";
	else if (index==AL_T) label = "AL_T";
	else if (index==RAYL) label = "Rayleigh-Number";
	else if (index==EK_V) label = "Vertical Ekman-Number";
	else if (index==EK_H) label = "Horizontal Ekman-Number";
	else if (index==ROSB) label = "Rossby-Number";
	else if (index==MIXP) label = "MIXP";
	else if (index==RESC) label = "RESC";
	else if (index==SPL1) label = "SPL1";
	else if (index==HMTP) label = "Homotopy";
	else if (index==SUNP) label = "Sun";
	else if (index==PE_V) label = "Vertical Peclet-Number";
	else if (index==PE_H) label = "Horizontal Peclet-Number";
	else if (index==SALT) label = "Salinity Forcing";
	else if (index==WIND) label = "Wind Forcing";
	else if (index==TEMP) label = "Temperature Forcing";
	else if (index==BIOT) label = "Nonlinear Factor";
	else if (index==COMB) label = "Combined Forcing";
	else if (index==NLES)   label = "NLES";
	else if (index==ARCL)   label = "ARCL";
	else if (index==IFRICB)   label = "IFRICB";
	else if (index==CONT)   label = "CONT";
	else if (index==P_VC)   label = "P_VC";
	else if (index==LAMB) label = "LAMB";
	else if (index==CMPR) label = "CMPR";
	else if (index==ALPC) label = "ALPC";
	else if (index==ENER) label = "Energy";
	else if (index==MKAP) label = "MKAP";
	else if (index==SPL2) label = "SPL2";
	else if (index==FPER) label = "FPER";
	else if (index==SPER) label = "SPER";
	else if (index==EXPO) label = "Exponent";
	else if (index==SEAS) label = "Seasonal Forcing";
	else if (index==SEASW) label = "Seasonal Forcing (Wind)";
	else if (index==SEAST) label = "Seasonal Forcing (Temperature)";
	else if (index==SEASS) label = "Seasonal Forcing (Salinity)";
	else if (index==MASS)  label = "Mass";
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
	if (param>0 && param<=_NPAR_) // time (0) and exp/seas (31/32) are not passed to THCM
    {
		FNAME(setparcs)(&param,&value);
    }
	else if (param<0)
    {
		ERROR("Invalid Parameter",__FILE__,__LINE__);
    }
	else if (param==0) // 0 is non-dimensional time
    {
		// set monthly forcing data
		bool time_dep_forcing = paramList.get("Time Dependent Forcing",false);
		if ((value>=0.0) && time_dep_forcing)
		{
			double gamma=1.0,gammaT=1.0,gammaS=1.0,gammaW=1.0;//default values
			this->getParameter("Seasonal Forcing",gamma);
			this->getParameter("Seasonal Forcing (Wind)",gammaW);
			this->getParameter("Seasonal Forcing (Temperature)",gammaT);
			this->getParameter("Seasonal Forcing (Salinity)",gammaS);
			gammaT*=gamma; gammaS*=gamma; gammaW*=gamma;
			INFO("Set THCM time to "<<value);
			INFO("Seasonal forcing parameters (W,T,S): "<<gammaW<<", "<<gammaT<<", "<<gammaS);
			F90NAME(m_monthly,update_forcing)(&value,&gammaW,&gammaT,&gammaS);
			if (internal_forcing)
			{
				F90NAME(m_monthly,update_internal_forcing)(&value,&gammaT,&gammaS);
			}

		}
		else if (value<0.0) //! reset to constant forcing (used in 4D FFT solver)
		{
			double gammaT=0.0,gammaS=0.0,gammaW=0.0,val=0.0;
			INFO("Set THCM forcing to constant");
			F90NAME(m_monthly,update_forcing)(&val,&gammaW,&gammaT,&gammaS);
			if (internal_forcing)
			{
				F90NAME(m_monthly,update_internal_forcing)(&val,&gammaT,&gammaS);
			}
		}
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

	// create an non-overlapping distributed map
	int i0 = domain->FirstRealI()+1; // 'grid-style' indexing is 1-based
	int i1 = domain->LastRealI()+1;
	int j0 = domain->FirstRealJ()+1;
	int j1 = domain->LastRealJ()+1;
	int k0 = domain->FirstRealK()+1;
	int k1 = domain->LastRealK()+1;

	// add global boundary cells
	if (i0==1) i0--; if (i1==n) i1++;
	if (j0==1) j0--; if (j1==m) j1++;
	if (k0==1) k0--; if (k1==l+la) k1++;

	int I0 = 0; int I1 = n+1;
	int J0 = 0; int J1 = m+1;
	int K0 = 0; int K1 = l+la+1;

	DEBUG("create landmap without overlap...");
	Teuchos::RCP<Epetra_Map> landmap_loc0 = Utils::CreateMap(i0,i1,j0,j1,k0,k1,
															 I0,I1,J0,J1,K0,K1,*Comm);

	// create an overlapping distributed map
	i0 = domain->FirstI()+1; // 'grid-style' indexing is 1-based
	i1 = domain->LastI()+1;
	j0 = domain->FirstJ()+1;
	j1 = domain->LastJ()+1;
	k0 = domain->FirstK()+1;
	k1 = domain->LastK()+1;

	//add the boundary cells i=0,n+1 etc (this is independent of overlap)
	i0--; i1++; j0--; j1++; k0--; k1++;

	DEBUG("create landmap with overlap...");
	Teuchos::RCP<Epetra_Map> landmap_loc = Utils::CreateMap(i0,i1,j0,j1,k0,k1,
															I0,I1,J0,J1,K0,K1,*Comm);

	DEBUG("Create local vectors...");

	// distributed non-overlapping version of landm
	Teuchos::RCP<Epetra_IntVector> landm_loc0 = Teuchos::rcp(new Epetra_IntVector(*landmap_loc0));
	// distributed overlapping version of landm
	Teuchos::RCP<Epetra_IntVector> landm_loc = Teuchos::rcp(new Epetra_IntVector(*landmap_loc));

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
// implement integral condition for S in Jacobian and B-matrix
void THCM::intcond_S(Epetra_CrsMatrix& A, Epetra_Vector& B)
{
    int N=domain->GlobalN();
    int M=domain->GlobalM();
    int L=domain->GlobalL();

    int lastrow = rowintcon_;

    int root = Comm->NumProc()-1;

    Teuchos::RCP<Epetra_MultiVector> intcond_glob =
        Utils::Gather(*intcond_coeff,root);

    if (A.MyGRID(lastrow))
	{
		if (Comm->MyPID()!=root)
        {
			ERROR("S-integral condition should be on last processor!",__FILE__,__LINE__);
        }
		int lid = B.Map().LID(lastrow);
		B[lid] = 0.0; // no more time-dependence for this S-point
		int len = N*M*L;
		double *values = new double[len];
		int *indices = new int[len];


		int pos=0;
		for (int i=0;i<N;i++)
			for (int j=0;j<M;j++)
				for (int k=0;k<L;k++)
				{
					int gid = FIND_ROW2(_NUN_,N,M,L,i,j,k,SS);
					indices[pos] = gid;
					values[pos] = (*intcond_glob)[0][gid];
					pos++;
				}

		/*
		  len=1;
		  indices[0]=lastrow;
		  values[0]=1.0;
		*/
		if (A.Filled())
        {
			CHECK_NONNEG(A.ReplaceGlobalValues(lastrow,len,values,indices));
        }
		else
        {
			CHECK_NONNEG(A.InsertGlobalValues(lastrow,len,values,indices));
        }

		delete []  values;
		delete []  indices;
	}
    else if (Comm->MyPID()==root)
	{
		ERROR("S-integral condition should be on last processor!",__FILE__,__LINE__);
	}
}

//=============================================================================
void THCM::fixPressurePoints(Epetra_CrsMatrix& A, Epetra_Vector& B)
{
	for (int i=1;i<=2;i++) {
		int row = (i==1)? rowPfix1: rowPfix2;
		if (A.MyGRID(row))
		{
			int lid = B.Map().LID(row);
			B[lid] = 0.0; // no more time-dependence for this P-point
			double values;
			int indices;

			int len=1;
			indices=row;
			values=1.0;

			if (A.Filled())
			{
				CHECK_NONNEG(A.ReplaceGlobalValues(row,len,&values,&indices));
			}
			else
			{
				CHECK_NONNEG(A.InsertGlobalValues(row,len,&values,&indices));
			}
		}
	}//for two pressure dirichlet values
}

//=============================================================================
Teuchos::RCP<Epetra_CrsGraph> THCM::CreateMaximalGraph()
{
	DEBUG("Constructing maximal matrix graph...");
	int n=domain->LocalN();
	int m=domain->LocalM();
	int l=domain->LocalL();
	int ndim = StandardMap->NumMyElements();
	int *numEntriesPerRow = new int[ndim];
	//int *landm = new int[n*m*l];
	//F90NAME(m_thcm_utils,get_landm)(landm);
	//const int LAND=1;
	for (int k=1;k<=l;k++)
		for (int j=1;j<=m;j++)
			for (int i=1;i<=n;i++)
			{
				int lidU = FIND_ROW2(_NUN_,n,m,l,i-1,j-1,k-1,UU);
				int gidU = AssemblyMap->GID(lidU);
				if (StandardMap->MyGID(gidU)) // otherwise: ghost cell, not in Jacobian
				{
					int lid0 = StandardMap->LID(gidU)-1;
					numEntriesPerRow[lid0+UU] = 24;
					numEntriesPerRow[lid0+VV] = 22;
					numEntriesPerRow[lid0+WW] = 7;
					numEntriesPerRow[lid0+PP] = 11;
					numEntriesPerRow[lid0+TT] = 20;
					numEntriesPerRow[lid0+SS] = 20;
				}
			}

	Teuchos::RCP<Epetra_CrsGraph> graph
		= Teuchos::rcp(new Epetra_CrsGraph(Copy,*StandardMap,numEntriesPerRow,false));

	DEBVAR(rowintcon_);
	DEBVAR(sres);

	int indices[24];
	int N = domain->GlobalN();
	int M = domain->GlobalM();
	int L = domain->GlobalL();

	int I0 = domain->FirstRealI();
	int J0 = domain->FirstRealJ();
	int K0 = domain->FirstRealK();
	int I1 = domain->LastRealI();
	int J1 = domain->LastRealJ();
	int K1 = domain->LastRealK();
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
#define DEBUG_GRAPH_ROW(var)											\
				std::cout << "graph row "<<i<<" "<<j<<" "<<k<<" "<<var;	\
				std::cout << " (gid "<<(gid0+var)<<"):"<<std::endl;		\
				std::cout << "predicted length: "<<numEntriesPerRow[StandardMap->LID(gid0+var)]; \
				std::cout <<", actual length: "<<pos<<std::endl;		\
				for (int pp=0;pp<pos;pp++) std::cout << (indices)[pp]<<" ";	\
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
				if ((sres==0)&&((gid0+SS)==rowintcon_)) break;
#endif
				DEBUG_GRAPH_ROW(SS)
					CHECK_ZERO(graph->InsertGlobalIndices(gid0+SS,pos,indices));
			}

#ifndef NO_INTCOND
	if (sres==0)
    {
		int grid = rowintcon_;
		if (StandardMap->MyGID(grid))
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
	if (domain->IsPeriodic())
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
// set vmix_fix
void THCM::fixMixing(int value)
{
    if (vmix_GLB)
	{
		INFO("setting vmix_fix to "<<value);
		F90NAME(m_mix,set_vmix_fix)(&value);
	}
}

//=============================================================================
void THCM::SetupMonthlyForcing()
{
	DEBUG("Initialize monthly Levitus...");

	// maps for 2D fields (surface forcing)
	Teuchos::RCP<Epetra_Map> lev_map_dist = domain->CreateStandardMap(1,true);
	Teuchos::RCP<Epetra_Map> lev_map_root = Utils::Gather(*lev_map_dist,0);

	// maps for 3D fields (internal forcing)
	Teuchos::RCP<Epetra_Map> intlev_map_dist = domain->CreateStandardMap(1,false);
	Teuchos::RCP<Epetra_Map> intlev_map_root = Utils::Gather(*intlev_map_dist,0);

	// create sequential and parallel vectors to hold the data
	Teuchos::RCP<Epetra_MultiVector> temp_glob =
		Teuchos::rcp(new Epetra_Vector(*intlev_map_root));
	Teuchos::RCP<Epetra_MultiVector> salt_glob =
		Teuchos::rcp(new Epetra_Vector(*intlev_map_root));

	Teuchos::RCP<Epetra_MultiVector> tatm_glob =
		Teuchos::rcp(new Epetra_Vector(*lev_map_root));
	Teuchos::RCP<Epetra_MultiVector> emip_glob =
		Teuchos::rcp(new Epetra_Vector(*lev_map_root));

	Teuchos::RCP<Epetra_MultiVector> taux_glob =
		Teuchos::rcp(new Epetra_Vector(*lev_map_root));
	Teuchos::RCP<Epetra_MultiVector> tauy_glob =
		Teuchos::rcp(new Epetra_Vector(*lev_map_root));

	// get raw pointers to the data:
	double *tatm_g, *emip_g, *taux_g, *tauy_g, *temp_g, *salt_g;
	CHECK_ZERO((*tatm_glob)(0)->ExtractView(&tatm_g));
	CHECK_ZERO((*emip_glob)(0)->ExtractView(&emip_g));
	CHECK_ZERO((*temp_glob)(0)->ExtractView(&temp_g));
	CHECK_ZERO((*salt_glob)(0)->ExtractView(&salt_g));
	CHECK_ZERO((*taux_glob)(0)->ExtractView(&taux_g));
	CHECK_ZERO((*tauy_glob)(0)->ExtractView(&tauy_g));

	// now create distributed maps
	Teuchos::RCP<Epetra_Map> lev_map_loc    = domain->CreateAssemblyMap(1,true);
	Teuchos::RCP<Epetra_Map> intlev_map_loc = domain->CreateAssemblyMap(1,false);

	// and distributed vectors
	Teuchos::RCP<Epetra_Vector> tatm_loc = Teuchos::rcp(new Epetra_Vector(*lev_map_loc));
	Teuchos::RCP<Epetra_Vector> emip_loc = Teuchos::rcp(new Epetra_Vector(*lev_map_loc));
	Teuchos::RCP<Epetra_Vector> temp_loc = Teuchos::rcp(new Epetra_Vector(*intlev_map_loc));
	Teuchos::RCP<Epetra_Vector> salt_loc = Teuchos::rcp(new Epetra_Vector(*intlev_map_loc));
	Teuchos::RCP<Epetra_Vector> taux_loc = Teuchos::rcp(new Epetra_Vector(*lev_map_loc));
	Teuchos::RCP<Epetra_Vector> tauy_loc = Teuchos::rcp(new Epetra_Vector(*lev_map_loc));

	// extract pointers to the distributed data
	double* ctatm,*cemip,*ctaux,*ctauy,*ctemp,*csalt;
	CHECK_ZERO(tatm_loc->ExtractView(&ctatm));
	CHECK_ZERO(emip_loc->ExtractView(&cemip));
	CHECK_ZERO(temp_loc->ExtractView(&ctemp));
	CHECK_ZERO(salt_loc->ExtractView(&csalt));
	CHECK_ZERO(taux_loc->ExtractView(&ctaux));
	CHECK_ZERO(tauy_loc->ExtractView(&ctauy));

	// to import overlap
	Teuchos::RCP<Epetra_Import> loc2dist =
		Teuchos::rcp(new Epetra_Import(*lev_map_loc,*lev_map_dist));
	Teuchos::RCP<Epetra_Import> int_loc2dist =
		Teuchos::rcp(new Epetra_Import(*intlev_map_loc,*intlev_map_dist));

	for (int month=1;month<=12;month++)
    {
		if (Comm->MyPID()==0)
		{
			F90NAME(m_global,get_monthly_forcing)(tatm_g,emip_g,taux_g,tauy_g,&month);
			if (internal_forcing)
			{
				F90NAME(m_global,get_monthly_internal_forcing)(temp_g,salt_g,&month);
			}
		}

		// distribute levitus and wind fields
		Teuchos::RCP<Epetra_MultiVector> tatm_dist = Utils::Scatter(*tatm_glob,*lev_map_dist);
		Teuchos::RCP<Epetra_MultiVector> emip_dist = Utils::Scatter(*emip_glob,*lev_map_dist);
		Teuchos::RCP<Epetra_MultiVector> temp_dist = Utils::Scatter(*temp_glob,*intlev_map_dist);
		Teuchos::RCP<Epetra_MultiVector> salt_dist = Utils::Scatter(*salt_glob,*intlev_map_dist);
		Teuchos::RCP<Epetra_MultiVector> taux_dist = Utils::Scatter(*taux_glob,*lev_map_dist);
		Teuchos::RCP<Epetra_MultiVector> tauy_dist = Utils::Scatter(*tauy_glob,*lev_map_dist);

		CHECK_ZERO(tatm_loc->Import(*tatm_dist,*loc2dist,Insert));
		CHECK_ZERO(emip_loc->Import(*emip_dist,*loc2dist,Insert));
		CHECK_ZERO(temp_loc->Import(*temp_dist,*int_loc2dist,Insert));
		CHECK_ZERO(salt_loc->Import(*salt_dist,*int_loc2dist,Insert));
		CHECK_ZERO(taux_loc->Import(*taux_dist,*loc2dist,Insert));
		CHECK_ZERO(tauy_loc->Import(*tauy_dist,*loc2dist,Insert));

		double tatmmax, emipmax;
		double tatmmin, emipmin;
		CHECK_ZERO(tatm_dist->MaxValue(&tatmmax));
		CHECK_ZERO(emip_dist->MaxValue(&emipmax));
		CHECK_ZERO(tatm_dist->MinValue(&tatmmin));
		CHECK_ZERO(emip_dist->MinValue(&emipmin));

		INFO("Month: " << month);
		INFO("Temperature-forcing range: [" << tatmmin << ".." << tatmmax << "]");
		INFO("Salinity-forcing range: [" << emipmin << ".." << emipmax << "]");

		F90NAME(m_monthly,set_forcing)(ctatm,cemip,ctaux,ctauy,&month);
		if (internal_forcing)
		{
			F90NAME(m_monthly,set_internal_forcing)(ctemp,csalt,&month);
		}
    }
}


//=============================================================================
extern "C" {

// this is a cheat for the fortran routine fsint from forcing.F90
#ifdef HUYGENS
	void thcm_forcing_integral(double* qfun2, double* y, int* landm, double* fsint)
#else
		void thcm_forcing_integral_(double* qfun2, double* y, int* landm, double* fsint)
#endif
	{
		Teuchos::RCP<Epetra_Comm> comm = THCM::Instance().GetComm();
		Teuchos::RCP<TRIOS::Domain> domain = THCM::Instance().GetDomain();

		int n = domain->LocalN();
		int m = domain->LocalM();
		int l = domain->LocalL();

		double lsint = 0.0, lfsint=0.0, sint;

		int i0 = domain->FirstRealI()-domain->FirstI();
		int j0 = domain->FirstRealJ()-domain->FirstJ();

		int i1 = domain->LastRealI()-domain->FirstI();
		int j1 = domain->LastRealJ()-domain->FirstJ();

		for (int j=j0; j<=j1; j++)
		{
			for (int i=i0;i<=i1;i++)
			{
				//note: the fortran landm array is 0-based, so we add a 1 to i,j,k
				int pl = FIND_ROW2(1,n+2,m+2,l+2,i+1,j+1,l,1);
				int pq = FIND_ROW2(1,n,m,1,i,j,0,1);
				lfsint = qfun2[pq] * cos(y[j]) * (1-landm[pl]) + lfsint;
				lsint = cos(y[j]) * (1-landm[pl]) + lsint;
			}
		}
		CHECK_ZERO(comm->SumAll(&lfsint,fsint,1));
		CHECK_ZERO(comm->SumAll(&lsint,&sint,1));
		*fsint = *fsint/sint;
		INFO("Flux correction equals "<<*fsint);
	}

	Teuchos::RCP<const Epetra_MultiVector> THCM::getNullSpace()
	{
		if (nullSpace==Teuchos::null)
		{
			nullSpace = Teuchos::rcp(new Epetra_MultiVector
									 (*StandardMap,2,true) );

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
			for (int k=domain->FirstRealK();k<=domain->LastRealK();k++)
				for (int j=domain->FirstRealJ();j<=domain->LastRealJ();j++)
					for (int i=domain->FirstRealI();i<=domain->LastRealI();i++)
					{
						if ((i+j)%2)
						{
							(*(*nullSpace)(0))[pos] = 1;
						}
						else
						{
							(*(*nullSpace)(1))[pos] = 1;
						}
						pos+=_NUN_;
					}
			double nrm1,nrm2;
			CHECK_ZERO((*nullSpace)(0)->Norm2(&nrm1));
			CHECK_ZERO((*nullSpace)(1)->Norm2(&nrm2));
			CHECK_ZERO((*nullSpace)(0)->Scale(1.0/nrm1));
			CHECK_ZERO((*nullSpace)(1)->Scale(1.0/nrm2));
		}
		return nullSpace;
	}
};
