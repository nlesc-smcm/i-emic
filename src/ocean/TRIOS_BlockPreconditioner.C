/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#include "Teuchos_Utils.hpp"
#include <sstream>
#include "Epetra_Map.h"

#include "TRIOS_Macros.H"

#include "TRIOS_BlockPreconditioner.H"
#include "TRIOS_Saddlepoint.H"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrixTransposer.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StrUtils.hpp>
#include "AztecOO.h"
#include "Epetra_Time.h"
#include "AztecOO_string_maps.h"
#include <iomanip>
#include "Teuchos_oblackholestream.hpp"

#include "Utils.H"

#include "TRIOS_SolverFactory.H"
#include "TRIOS_Domain.H"

#include "TRIOS_Static.H"

#include "THCMdefs.H"
#include "GlobalDefinitions.H"

/// define this to set P=I, just remove checkerboard pressure modes
//#define DUMMY_PREC 1

using Teuchos::rcp_dynamic_cast;

// this reordering of the T/S system was intended
// for testing direct solvers like SuperLU_dist,
// but it has not really been tested in practice
//#define LINEAR_ARHOMU_MAPS

// here are some macros controling which parts
// of the preconditioner should be employed.
// This is only for testing purposes!

#ifdef TESTING
// error tolerance when checking i.e. result of linear solves etc.
#define _TESTTOL_ 1e-14
#endif


namespace TRIOS {

///////////////////////////////////////////////////////////////////////////////
// CLASS OCEANPREC
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////

    BlockPreconditioner::BlockPreconditioner(Teuchos::RCP<Epetra_CrsMatrix> jac,
                                             Teuchos::RCP<Domain> domain,
                                             Teuchos::ParameterList &List)
        :
        label_("Ocean Preconditioner"),
        jacobian(jac),
        domain(domain),
        needs_setup(true),
        IsComputed_(false)
    {
        INFO("Create new ocean preconditioner...");
        comm = domain->GetComm();
        this->SetParameters(List);
        this->Setup1();
    }

///////////////////////////////////////////////////////////////////////////////
// destructor
///////////////////////////////////////////////////////////////////////////////

    BlockPreconditioner::~BlockPreconditioner()
    {
        if (verbose>5)
        {
            INFO("Destroy ocean preconditioner...");
        }
        // not sure this loop is necessry but it certainly doesn't hurt

        for (int i =0; i < _NUMSUBM; i++)
        {
            SubMatrix[i]=Teuchos::null;
            SubMatrixRowMap[i]=Teuchos::null;
            SubMatrixColMap[i]=Teuchos::null;
            SubMatrixRangeMap[i]=Teuchos::null;
            SubMatrixDomainMap[i]=Teuchos::null;
            Importer[i]=Teuchos::null;
        }

        if (is_dummyP!=NULL) delete [] is_dummyP;
        if (is_dummyW!=NULL) delete [] is_dummyW;

        // rest should be handled automatically by Teuchos::rcp's
    }

    // to be called in constructors:
    void BlockPreconditioner::Setup1()
    {
        if (verbose > 5) INFO("Enter Setup1()...");
        // set pointer to domain map (map of vector x in y=A*x)
        // this preconditioner acts on standard vectors, i.e.
        // with 'nun' (6) unknowns per node
        domainMap = domain->GetSolveMap();

        // the map for y is the same
        rangeMap = domainMap;

        Teuchos::RCP<Epetra_Map> RowMap = domain->GetSolveMap();

        // column maps are maps with every grid point
        // present on every processor. in Epetra, the
        // column map of a matrix indicates which
        // column indices can possibly occur on a
        // subdomain. The reason why we have to make
        // column maps at all is that some of the
        // blocks map into different variables, i.e.
        // BwTS: TS -> etc.
        Teuchos::RCP<Epetra_Map> ColMap = domain->GetColMap();

        // split up row map:
        if (verbose>5) INFO("$   Split main map into uv|w|p|TS maps...");

        //  DEBVAR(*RowMap);

        const int labelUV[2] = {UU, VV};
        const int labelTS[2] = {TT, SS};

        mapUV = Utils::CreateSubMap(*RowMap, dof_, labelUV);
        mapW  = Utils::CreateSubMap(*RowMap, dof_, WW);
        mapP  = Utils::CreateSubMap(*RowMap, dof_, PP);
        mapTS = Utils::CreateSubMap(*RowMap, dof_, labelTS);

        // split up column map:

        if (verbose>5)
        {
            INFO("$   Build corresponding column maps...");
        }

        colmapUV = Utils::CreateSubMap(*ColMap, dof_, labelUV);
        colmapW  = Utils::CreateSubMap(*ColMap, dof_, WW);
        colmapP  = Utils::CreateSubMap(*ColMap, dof_, PP);
        colmapTS = Utils::CreateSubMap(*ColMap, dof_, labelTS);

        if (verbose>5)
        {
            INFO("$   Create Importers ...");
        }

        // construct import objects for row maps:

        // note: the terminology of import objects is very confusing:
        // Epetra_Import(target_map,source_map) means that we wish to
        // import values from the target map into the source map. The
        // corresponding function call is A.Export(B,Importer,...),
        // which means that the values of A are replaced by corresponding
        // values in B. In this preconditioner class we use importers not
        // for getting non-local data but mainly to extract certain rows/columns
        // from matrices (a permutation operation which doesn't involve any
        // communication)

        // Note that the same operation could be done with Export objects, and/or by
        // calling 'Import' instead of 'Export'. To be honest, there is no clear line
        // in this file, it is basically trial and error for every new import/export.

        importUV = Teuchos::rcp(new Epetra_Import(*RowMap, *mapUV) );
        importW  = Teuchos::rcp(new Epetra_Import(*RowMap,  *mapW) );
        importP  = Teuchos::rcp(new Epetra_Import(*RowMap,  *mapP) );
        importTS = Teuchos::rcp(new Epetra_Import(*RowMap, *mapTS) );

        Mzp1 = Teuchos::null;
        Mzp2 = Teuchos::null;

        // will be constructed when the preconditioner is computed
        AuvSolver  = Teuchos::null;
        ATSSolver  = Teuchos::null;
        SppSolver  = Teuchos::null;
        AuvPrecond = Teuchos::null;
        ATSPrecond = Teuchos::null;
        SppPrecond = Teuchos::null;

        is_dummyP=NULL;
        is_dummyW=NULL;

        QTS = Teuchos::null;

        if (verbose>5)
        {
            // this is called at the end of any constructor, so this message makes sense:
            INFO("BlockPreconditioner constructor completed.");
        }

// note: Once the Jacobian is available we still have to
// do some more setup, that's why needs_setup == true at this point
    }

// build preconditioner. This function is called by NOX/LOCA,
// I am not sure what the arguments mean but we do not use them
// anyway
    bool BlockPreconditioner::computePreconditioner(const Epetra_Vector& x,
                                                    Epetra_Operator& Prec,
                                                    Teuchos::ParameterList* p)
    {
        int ierr = this->Compute();
        return (ierr == 0);
    }


///////////////////////////////////////////////////////////////////////////////
// private setup function. Can only be done once the actual Jacobian
// is there (or at least it's pattern)
///////////////////////////////////////////////////////////////////////////////

    void BlockPreconditioner::Setup2()
    {
        if (verbose>5)
        {
            INFO("Enter Setup2()...");
        }
        // create some more maps we need for the submatrices and set
        // pointer arrays, create the matrices.

        // detect dummy w and p points

        int NmapW = mapW->NumMyElements();
        int NmapP = mapP->NumMyElements();
        is_dummyW = new bool[NmapW];
        is_dummyP = new bool[NmapP];

        INFO("mapP->NumMyElements = " << NmapP);
        INFO("mapW->NumMyElements = " << NmapW);

        // dummy points are defined to be such rows of the matrix that
        // a) belong to the corresponding eqn (p or w, depending on the map)
        // b) have exactly one entry, namely a 1 on the diagonal
        //
        // Physically these rows occur
        // - in land cells for P and W
        // - in ocean cells below atmosphere (or land) cells for w
        //
        // so that in practice there are always m*n more dummy w's than dummy p's
        // (the top ocean layer)
        dump = detect_dummies(*jacobian, *mapP, is_dummyP);
        dumw = detect_dummies(*jacobian, *mapW, is_dummyW);

        INFO(" Remove dummy W and P points ...");
        INFO(" dummy ws: " << dumw);
        INFO(" dummy ps: " << dump);

        INFO(" Create maps P1, W1, P^ ...");

        // create maps without dummy points: All matrices and vectors
        // will be based on these maps, not the original P and W maps
        // (which means that dummy points are completely ignored by
        // the preconditioner)
        mapP1 = Utils::CreateSubMap(*mapP,is_dummyP);
        colmapP1 = Utils::AllGather(*mapP1);

        mapW1 = Utils::CreateSubMap(*mapW,is_dummyW);
        colmapW1 = Utils::AllGather(*mapW1);

        // this map contains all P points corresponding to non-dummy W's.
        // That is, the P's in the top ocean layer are removed, in addition
        // to P's in land cells (original dummy P's)         ^
        // (in the paper the corresponding vector is denoted u_p).
        mapPhat = Utils::CreateSubMap(*mapP,is_dummyW);

        // create importers
        importW1 = Teuchos::rcp(new Epetra_Import(*(domain->GetSolveMap()), *mapW1) );
        importP1 = Teuchos::rcp(new Epetra_Import(*(domain->GetSolveMap()), *mapP1) );

        // note: The target map for this importer is P1, not the standard solve map.
        //       It is used to extract the relevant part of a P-vector.
        importPhat = Teuchos::rcp(new Epetra_Import(*mapPhat, *mapP1) );


        // create the map for \bar{p}, the surface (or 'depth-averaged') pressure
        // in ocean cells of the top layer:
        int l = domain->LocalL();

        int np = mapP->NumMyElements();
        int npbar = np/l;
        // pointer to 'top layer' in dummy array (k=l)
        bool *is_dummy_Pbar = &(is_dummyP[np-npbar]);

        int NumMyElements = npbar;
        int end = npbar;
        for (int i=0; i<end;i++)
        {
            if (is_dummy_Pbar[i]) NumMyElements--;
        }
        mapPbar = Teuchos::rcp(new Epetra_Map(-1,NumMyElements,0,*comm));

        // Set arrays with pointers. This makes accessing the objects
        // more convenient later on.
        INFO(" Set Pointer Arrays ...");

        // diagonal blocks:

        // Auv: uv->uv
        SubMatrixRowMap[_Auv]    = mapUV;
        SubMatrixColMap[_Auv]    = colmapUV;
        SubMatrixRangeMap[_Auv]  = mapUV;
        SubMatrixDomainMap[_Auv] = mapUV;
        Importer[_Auv]           = importUV;
        SubMatrixLabel[_Auv]     = "Auv";

        // Dw is not really needed, only its 'square part' Aw
        // Dw: W1->P1
        // Aw: W1->W1
        SubMatrixRowMap[_Dw]    = mapP1;
        SubMatrixColMap[_Dw]    = colmapW1;
        SubMatrixRangeMap[_Dw]  = mapP1;
        SubMatrixDomainMap[_Dw] = mapW1;
        Importer[_Dw]           = importP1;
        SubMatrixLabel[_Dw]     = "Dw";

        // ATS: TS->TS
        SubMatrixRowMap[_ATS]    = mapTS;
        SubMatrixColMap[_ATS]    = colmapTS;
        SubMatrixRangeMap[_ATS]  = mapTS;
        SubMatrixDomainMap[_ATS] = mapTS;
        Importer[_ATS]           = importTS;
        SubMatrixLabel[_ATS]     = "ATS";

        // B matrices

        // B_{uv} in TS equation: UV->TS
        SubMatrixRowMap[_BTSuv]    = mapTS;
        SubMatrixColMap[_BTSuv]    = colmapUV;
        SubMatrixRangeMap[_BTSuv]  = mapTS;
        SubMatrixDomainMap[_BTSuv] = mapUV;
        Importer[_BTSuv]           = importTS;
        SubMatrixLabel[_BTSuv]     = "BTSuv";

        // B_{TS} in w equation TS->W1
        SubMatrixRowMap[_BwTS]    = mapW1;
        SubMatrixColMap[_BwTS]    = colmapTS;
        SubMatrixRangeMap[_BwTS]  = mapW1;
        SubMatrixDomainMap[_BwTS] = mapTS;
        Importer[_BwTS]           = importW1;
        SubMatrixLabel[_BwTS]     = "BwTS";

        // B_w in TS eqn: W1->TS
        SubMatrixRowMap[_BTSw] = mapTS;
        SubMatrixColMap[_BTSw] = colmapW1;
        SubMatrixRangeMap[_BTSw] = mapTS;
        SubMatrixDomainMap[_BTSw] = mapW1;
        Importer[_BTSw] = importTS;
        SubMatrixLabel[_BTSw]="BTSw";

        //G and D matrices

        // pressure gradient, Guv: P1->UV
        SubMatrixRowMap[_Guv] = mapUV;
        SubMatrixColMap[_Guv] = colmapP1;
        SubMatrixRangeMap[_Guv] = mapUV;
        SubMatrixDomainMap[_Guv] = mapP1;
        Importer[_Guv] = importUV;
        SubMatrixLabel[_Guv]="Guv";

        // Gw without dummies: P1->W1
        SubMatrixRowMap[_Gw] = mapW1;
        SubMatrixColMap[_Gw] = colmapP1;
        SubMatrixRangeMap[_Gw] = mapW1;
        SubMatrixDomainMap[_Gw] = mapP1;
        Importer[_Gw] = importW1;
        SubMatrixLabel[_Gw]="Gw";

        // divergence in p-eqn Duv: UV->P1
        SubMatrixRowMap[_Duv]    = mapP1;
        SubMatrixColMap[_Duv]    = colmapUV;
        SubMatrixRangeMap[_Duv]  = mapP1;
        SubMatrixDomainMap[_Duv] = mapUV;
        Importer[_Duv]           = importP1;
        SubMatrixLabel[_Duv]     = "Duv";


        // allocate memory for submatrices
        for (int i=0;i<_NUMSUBM;i++)
        {
            SubMatrix[i] =
                Teuchos::rcp(new Epetra_CrsMatrix(Copy,
                                                  *(SubMatrixRowMap[i]),
                                                  *(SubMatrixColMap[i]),0));

            SubMatrix[i]->SetLabel(SubMatrixLabel[i].c_str());
        }

        INFO(" Build singular vectors of the pressure...");

        build_svp();

        INFO(" BlockPreconditioner setup done.");

        // needs_setup is still kept 'true'. After the submatrices have been extracted
        // for the first time their column maps will be adjusted, THEN needs_setup will
        // be false.

    }//Setup2()



///////////////////////////////////////////////////////////////////////////////
// create singular vectors of the pressure field.
// These are "checkerboard" modes which have to be
// explicitly removed from the pressure solution
///////////////////////////////////////////////////////////////////////////////

    void BlockPreconditioner::build_svp()
    {
        svp1 = Teuchos::rcp(new Epetra_Vector(*mapP1) );
        svp2 = Teuchos::rcp(new Epetra_Vector(*mapP1) );

        DEBUG("Build pressure vectors svp1/2...");

        // temporary vectors used for assembly

        // rectangular domain without overlap, only 1 DoF
        DEBUG("Calling domain->CreateStandardMap");
        Teuchos::RCP<Epetra_Map> stdMap = domain->CreateStandardMap(1,false);
        Teuchos::RCP<Epetra_Vector> std1 = Teuchos::rcp(new Epetra_Vector(*stdMap));
        Teuchos::RCP<Epetra_Vector> std2 = Teuchos::rcp(new Epetra_Vector(*stdMap));

        // load-balanced map, only 1 DoF
        DEBUG("Calling domain->CreateSolveMap");
        Teuchos::RCP<Epetra_Map> slvMap = domain->CreateSolveMap(1,false);

        // transfer operator
        Teuchos::RCP<Epetra_Import> import =
            Teuchos::rcp(new Epetra_Import(*stdMap,*slvMap));

        int pos = 0;

        // the svp's are fairly easy to construct, they are
        // so-called 'checkerboard' modes' in the x-y planes.
        // we first construct them for the standard rectan-
        // gular subdomains and then export them to the load-
        // balanced 'solve' map.

        // loop over all non-ghost subdomain cells:
        for (int k=domain->FirstRealK();k<=domain->LastRealK();k++)
            for (int j=domain->FirstRealJ();j<=domain->LastRealJ();j++)
                for (int i=domain->FirstRealI();i<=domain->LastRealI();i++)
                {
                    if ((i+j)%2)
                    {
                        (*std1)[pos++] = 1;
                    }
                    else
                    {
                        (*std2)[pos++] = 1;
                    }
                }

        // redistribute to solve map
        Teuchos::RCP<Epetra_Vector> slv1 = Teuchos::rcp(new Epetra_Vector(*slvMap));
        Teuchos::RCP<Epetra_Vector> slv2 = Teuchos::rcp(new Epetra_Vector(*slvMap));
        CHECK_ZERO(slv1->Export(*std1,*import,Insert));
        CHECK_ZERO(slv2->Export(*std2,*import,Insert));

        // replace our temporary map by the pressure map
        // (with global indices [3, 9, 15, ...]).
        CHECK_ZERO(slv1->ReplaceMap(*mapP));
        CHECK_ZERO(slv2->ReplaceMap(*mapP));

        // remove dummy points (i.e. land cells)
        import = Teuchos::rcp(new Epetra_Import(*mapP1,*mapP));
        CHECK_ZERO(svp1->Import(*slv1,*import,Zero));
        CHECK_ZERO(svp2->Import(*slv2,*import,Zero));


//  DEBVAR(*svp1);
//  DEBVAR(*svp2);

        // normalize the global svp's
        double nrm1, nrm2;
        CHECK_ZERO(svp1->Norm2(&nrm1));
        CHECK_ZERO(svp2->Norm2(&nrm2));
        CHECK_ZERO(svp1->Scale(1.0/nrm1));
        CHECK_ZERO(svp2->Scale(1.0/nrm2));

    }

#ifdef TESTING //[

// check svp's: They must satisfy Guv*svp=Gw*svp=0.
    bool BlockPreconditioner::test_svp()
    {
        INFO("Entering test_svp()");
        bool result = true;
        double norm;
        INFO("Test if svp1 and svp2 are singular vectors of P...");
        INFO("Check if Gw*svp1=0...");
        Epetra_Vector testW(*mapW1);
        CHECK_ZERO(SubMatrix[_Gw]->Multiply(false,*svp1,testW));
        CHECK_ZERO(testW.Norm2(&norm));
        if (norm>_TESTTOL_)
        {
            result = false;
            INFO("WARNING: ||Gw*svp1||_2 = " << norm);
            DEBVAR(*SubMatrix[_Gw]);
            DEBVAR(*svp1);
            DEBVAR(testW);
        }
        else
        {
            INFO("Test successful.");
        }
        INFO("Check if Gw*svp2=0...");
        CHECK_ZERO(SubMatrix[_Gw]->Multiply(false,*svp2,testW));
        CHECK_ZERO(testW.Norm2(&norm));
        if (norm>_TESTTOL_)
        {
            result = false;
            INFO("WARNING: ||Gw*svp2||_2 = "<< norm);
            DEBVAR(*SubMatrix[_Gw]);
            DEBVAR(*svp2);
            DEBVAR(testW);
        }
        else
        {
            INFO("Test successful.");
        }
        INFO("Check if Guv*svp1=0...");
        Epetra_Vector testUV(*mapUV);
        CHECK_ZERO(SubMatrix[_Guv]->Multiply(false,*svp1,testUV));
        CHECK_ZERO(testUV.Norm2(&norm));
        if (norm>_TESTTOL_)
        {
            result = false;
            INFO("WARNING: ||Guv*svp1||_2 = "<< norm );
            DEBVAR(*SubMatrix[_Guv]);
            DEBVAR(*svp1);
            DEBVAR(testUV);
        }
        else
        {
            INFO("Test successful.");
        }

        INFO("Check if Guv*svp2=0...");
        CHECK_ZERO(SubMatrix[_Guv]->Multiply(false,*svp2,testUV));
        CHECK_ZERO(testUV.Norm2(&norm));
        if (norm>_TESTTOL_)
        {
            result = false;
            INFO("WARNING: ||Gw*svp2||_2 = "<< norm );
            DEBVAR(*SubMatrix[_Guv]);
            DEBVAR(*svp2);
            DEBVAR(testUV);
        }
        else
        {
            INFO("Test successful.");
        }
        return result;
    }//test_svp

#endif //]TESTING

///////////////////////////////////////////////////////////////////////////////
// splits the jacobian into the various subblocks
///////////////////////////////////////////////////////////////////////////////

    void BlockPreconditioner::extract_submatrices(const Epetra_CrsMatrix& Jac)
    {
        if (verbose>5)
        {
            INFO("Extract submatrices..." << std::endl);
        }
        {

            // construct all submatrices
            for (int i = 0; i < _NUMSUBM; i++)
            {
                DEBUG("extract submatrix " << i);
                CHECK_ZERO(SubMatrix[i]->Export(Jac, *Importer[i], Zero));

                // the first time we do this we have to adjust the column maps
                if (needs_setup)
                {
                    // at first the col maps contain all possible nodes. Once we know
                    // what the submatrices look like (they always have the same structure)
                    // replace the col map by one that contains only cols actually present
                    // on the processor:
                    CHECK_ZERO(SubMatrix[i]->FillComplete(*SubMatrixDomainMap[i],
                                                          *SubMatrixRangeMap[i]));
                    
                    SubMatrixColMap[i] = Utils::CompressColMap(*SubMatrix[i]);
                    SubMatrix[i]=Utils::ReplaceColMap(SubMatrix[i],*SubMatrixColMap[i]);
                }
                CHECK_ZERO(SubMatrix[i]->FillComplete(*SubMatrixDomainMap[i],
                                                      *SubMatrixRangeMap[i]));
                CHECK_ZERO(SubMatrix[i]->OptimizeStorage());
                //DEBVAR(*SubMatrix[i]);
            }
        }



        // all operations that have to be done once the Jacobian
        // is available have now been performed:
        needs_setup = false;

        // change sign of Duv and Dw (the whole continuity equation is negated)
        CHECK_ZERO(SubMatrix[_Duv]->Scale(-1.0));
        CHECK_ZERO(SubMatrix[_Dw]->Scale(-1.0));

        // since we replace the matrices Auv and ATS by new ones (see next comment/commands),
        // there occurs a problem in ML (as of Trilinos 10), a segfault if we do not delete the
        // solver before the matrix. This is a bug in Trilinos and will probably be fixed soon.
        AuvPrecond=Teuchos::null;
        ATSPrecond=Teuchos::null;

        // Auv/ATS have to be Ifpack-safe. This is definitely
        // the case if we don't give any col map.
        // (I believe the 'LocalFilter' in Ifpack is buggy)
        Auv = Utils::RemoveColMap(SubMatrix[_Auv]);
        CHECK_ZERO(Auv->FillComplete(*mapUV,*mapUV));

        DEBUG("Adjust diagonal block ATS...");
        ATS = Utils::RemoveColMap(SubMatrix[_ATS]);
        CHECK_ZERO(ATS->FillComplete(*mapTS,*mapTS));

        // construct the matrices that are not in the arrays
        // (because they are not extracted directly from the Jacobian)

        DEBUG("Construct Aw: W1->W1 (square part of Dw: W1->P1)");

        //DEBVAR(*SubMatrix[_Dw])

        // Aw is the 'square part' of Dw, in practice it maps
        // to w vectors (which have less points) instead of p vectors:
        Aw = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *mapPhat, SubMatrix[_Dw]->ColMap(),
                                               SubMatrix[_Dw]->MaxNumEntries()) );

        CHECK_ZERO(Aw->Import(*SubMatrix[_Dw],*importPhat,Zero));

        // make it map to W1 instead of Phat:
// in Trilinos < 7.0.5, 'ReplaceRowMap' is buggy. Also, we need
// to make sure that row map = col map, so we do it by hand for now:
#if 0
        CHECK_ZERO(Aw->ReplaceRowMap(*mapW1));
#else
        CHECK_ZERO(Aw->FillComplete());
        Aw = Utils::ReplaceBothMaps(Aw,*mapW1,*mapW1);
#endif

        Aw->SetLabel("Aw");
        CHECK_ZERO(Aw->FillComplete(*mapW1,*mapW1));
        CHECK_ZERO(Aw->OptimizeStorage());

        //DEBVAR(*Aw);

        DEBUG("Construct Duv1: UV->W1 ('square part' of Duv)");

        // Duv1 is the 'square part' of Duv, in practice it maps
        // to w1 vectors (which have less points) instead of p1 vectors:

        Duv1 =
            Teuchos::rcp(new Epetra_CrsMatrix(Copy,
                                              *mapPhat,
                                              SubMatrix[_Duv]->ColMap(),
                                              SubMatrix[_Duv]->MaxNumEntries()) );
        CHECK_ZERO(Duv1->Import(*SubMatrix[_Duv],*importPhat,Zero));

        // make it map to W1 instead of Phat:
#if 0
        CHECK_ZERO(Duv1->ReplaceRowMap(*mapW1));
#else
        Duv1->FillComplete();
        Duv1 = Utils::ReplaceRowMap(Duv1,*mapW1);
#endif
        CHECK_ZERO(Duv1->FillComplete(*mapUV,*mapW1));
        Duv1->SetLabel("Duv1");

//    DEBVAR(*Duv1);
        DEBUG("All submatrices have been extracted");

#ifdef TESTING
        // The pressure vectors svp1,2 are built at the end of Setup2() and
        // they must fullfill the conditions Guv*svp=Gw*svp=0. Test this:
        INFO("Check singular vectors...");
        if (!test_svp())
            INFO("WARNING: Singular vectors did not pass test!");
#endif

    }


///////////////////////////////////////////////////////////////////////////////
// another setup function, find dummy cells w.r.t. p or w
///////////////////////////////////////////////////////////////////////////////


// count and label dummy P or W points (rows of A where there is only a 1 on the
// diagonal).
    int BlockPreconditioner::detect_dummies( const Epetra_CrsMatrix& A,
                                             const Epetra_Map& M, bool *is_dummy) const
    {
        int dim = M.NumMyElements();
        int ndummies = 0;
        int row;
        int len;
        //int col;
        //double val;
        int maxlen     = A.MaxNumEntries();
        int *indices   = new int[maxlen];
        double *values = new double[maxlen];

        len = 1;
        for (int i = 0; i < dim; i++)
        {
            row = M.GID(i);
            is_dummy[i] = false;

            CHECK_ZERO(A.ExtractGlobalRowCopy(row, maxlen, len, values, indices));

            for (int p = 0; p < len; p++)
            {
                if (indices[p] == row && std::abs(values[p]-1.0) < 1e-8)
                {
                    is_dummy[i] = true;
                }
                else if (values[p] != 0.0)
                {
                    is_dummy[i] = false;
                    break;
                }
            }

            if (is_dummy[i])
            {
                ndummies++;
            }
        }

        delete [] indices;
        delete [] values;
        return ndummies;
    }//detect_dummies






///////////////////////////////////////////////////////////////////////////////
// builds depth-averaging operator Mzp
///////////////////////////////////////////////////////////////////////////////

    Teuchos::RCP<Epetra_CrsMatrix> BlockPreconditioner::build_singular_matrix(Teuchos::RCP<Epetra_CrsMatrix> Gw)
    {

        //           _
        // Mzp: P -> P. The columns of Mzp are the singular vectors of Gw
        //
        // Since Gw is a bidiagonal matrix (discrete gradient) of the form
        //
        //  -k      k
        //    -k      k
        //      -k      k
        //
        // we can construct the Teuchos::nullvectors (singular vectors) explicitly
        // by traversing the matrix in a zickzack pattern:
        //
        //  -k  ->  k
        //    -k    | k
        //      -k  |   k
        //        -kv     k
        //          -k  ->  k

        // if the matrix is correctly built we can do this in parallel, as
        // there should only be connections in z-direction, which is not
        // split up among processors


#ifdef TESTING
        if (!Gw->Filled())
        {
            ERROR("build_singular_matrix: Matrix Gw is not ready yet!",__FILE__,__LINE__);
        }
#endif

        // start by creating a CRS matrix (adaptation of the Fortran routine
        // build_Mzp from m_bgskit)
        INFO(" build Mzp...");

        // determine dimension
        int nGw = Gw->NumMyRows();

        // upper bound for nonzeros in Mzp:
        int nnz = mapUV->NumGlobalElements();
        int np1 = mapP1->NumMyElements();
        DEBVAR(nGw);

        // allocate temporary CRS graph
        int* begM = new int[np1+1];
        int* jcoM = new int[nnz];
        double* coM = new double[nnz]; // not really essential, but convenient later on

        begM[nGw] = nnz;

        bool *track = new bool[np1]; // keep track of vectors that have been built
        for (int i = 0; i < np1; i++) track[i]=false;

        int svs = 0; // count singular vectors
        int j = 0;   // count nnz
        int k;
        int len;  // length of current row of Gw

        int *jcoGw; // column indices in current Gw row
        double *coGw; // value array of current Gw row

        begM[0]=0;

        int kglob = -1;
        double sf; //scaling factor

        // build all Teuchos::null-vectors of Gw. Note that Gw does not
        // contain the top-layer, so we might miss some: those
        // water-columns that are only one cell thick.
        for (int i = 0; i < nGw; i++)
        {
            //DEBVAR(i)
            if (!track[i]) // new singular vector
            {
                sf = 1.0;
                svs = svs+1;
                //DEBVAR(svs)
                k = i; // start in row i
                while (k < nGw) // jump through matrix
                {
                    // look at row k from the matrix Gw
                    CHECK_ZERO(Gw->ExtractMyRowView(k, len, coGw, jcoGw));
                    // there may be a single 0-entry in the top layer due to
                    // the fixed pattern of the Jacobian
                    if (len < 2) break;
#if 1
                    if (len!=2)
                    {
                        INFO("found unexpected local row in matrix Gw:");
                        for (int ii=0;ii<len;ii++)
                        {
                            INFO("("<<k<<", "<<jcoGw[ii]<<")\t"<<coGw[ii]);
                        }
                        ERROR("Gw should be a bi-diagonal matrix!",__FILE__,__LINE__);
                    }
#endif
                    //DEBVAR(k)

                    j = j+1; // counts nnz in M
                    //DEBVAR(j)
                    track[k] = true;
                    jcoM[j-1] = Gw->GCID(jcoGw[0]); // transform to global P index
                    coM[j-1] = sf;

                    // adjust scaling factor
                    sf = sf*(std::abs(coGw[0])/std::abs(coGw[1]));

                    // jump to the row of the superdiagonal column (>,V)

                    // this process is a little awkward because of the different indexing schemes:
                    // we want k to be the local index in mapW corresponding to the local index
                    // of the column in mapP. As the cell might be a dummy W but _not_ a dummy P,
                    // we store the P index seperately in kglob as it is otherwise lost.
                    kglob = Gw->GCID(jcoGw[1]); // kglob is the global index of the P node of the
                    // superdiagonal column
                    k = kglob - (PP-WW); // now k is the global index of the W-node
                    k = mapW1->LID(k); // now k is the local row index, or -1 if it's a dummy W cell
                    if (k<0) // this means that it is not in the map. If Gw is correct,
                    {      // this occurs only if the cell has a dummy w but not dummy p (top layer cell)
                        //mark the surface cell as tracked:
                        track[mapP1->LID(kglob)]=true;
                        break;   // go to next row i
                    }
                }
                j = j+1;
                //DEBVAR(j)
                jcoM[j-1] = kglob;
                coM[j-1] = sf;
                begM[svs] = j;
            }//track
        }//i

        // extra-loop over the top layer: Gw is not defined here,
        // but any columns that have not been tracked yet are
        // water columns of thickness 1 cell, where \bar{p}=p holds.
        INFO("enter extra-loop, ["<<nGw<<":"<<np1-1<<"]");
        for (int i=nGw;i<np1;i++)
        {
            if (!track[i]) // new singular vector
            {
                //DEBVAR(i)
                svs = svs+1;
                //DEBVAR(svs)
                j = j+1; // counts nnz in M
                //DEBVAR(j)
                track[i] = true;
                INFO("shallow layer node: ("<<i<<", "<<mapP1->GID(i)<<")");
                jcoM[j-1] = mapP1->GID(i); // transform to global P index
                coM[j-1] = 1.0;
                begM[svs] = j;
            }//track
        }//i

        INFO("final svs: "<<svs);

        Teuchos::RCP<Epetra_CrsMatrix> Mzp =
            Teuchos::rcp(new Epetra_CrsMatrix(Copy,*mapPbar,*colmapP1,domain->LocalL()) );

        int ipb; // row index in Pbar indexing

        //normalize rows and insert in Trilinos matrix:
        for (int i=0; i<svs; i++)
        {
            j = begM[i];
            k = begM[i+1];
            double scale = 0.0;
            int len = k-j;
            for (int jj=j;jj<k;jj++)
            {
                scale += (coM[jj])*(coM[jj]);
            }
            scale = 1.0/sqrt(scale);
            for (int jj=j;jj<k;jj++)
            {
                coM[jj]*=scale;
            }
            ipb = mapPbar->GID(i);

            if (ipb >= 0)
            {
                CHECK_ZERO(Mzp->InsertGlobalValues(ipb, len, coM+begM[i], jcoM+begM[i]));
            }
        }


        // pass in a domain map and a range map:
        CHECK_ZERO(Mzp->FillComplete(*mapP1,*mapPbar));

#if     1
        comm->Barrier();
        // make sure that Gw*Mzp'=0
        INFO("Depth-averaging matrix Mzp has been built.");

        INFO("Test if Gw*Mzp'=0 is satisfied...");
        Teuchos::RCP<Epetra_CrsMatrix> C = Teuchos::rcp(new Epetra_CrsMatrix(Copy,Gw->RangeMap(),10) );
        CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(*Gw, false, *Mzp,true,*C));
        if (C->NormInf()>1e-10)
        {
            INFO("WARNING: ||Gw*Mzp'||_inf = "<<C->NormInf());
            INFO("This should be very close to 0!");
            DEBUG("WARNING: Gw*Mzp' is NOT close to 0, here it comes:");
            DEBVAR(*C);
        }
        else
        {
            INFO("Test successful!");
        }
#endif

        delete [] track;
        delete [] begM;
        delete [] jcoM;
        delete [] coM;
        //DEBVAR(*Mzp);
        DEBUG("### leave BlockPreconditioner::build_singular_matrix ###");
        return Mzp;
    }


///////////////////////////////////////////////////////////////////////////////
// another setup function: build block systems, preconditioners and solvers
///////////////////////////////////////////////////////////////////////////////

    void BlockPreconditioner::build_preconditioner(void)
    {

        if (verbose>5)
        {
            INFO("Prepare preconditioner...");
        }
        {
            char ApType = lsParams.sublist("Ap Solver").get("Full or square", 'S');
            Ap = Teuchos::rcp(new ApMatrix(*SubMatrix[_Gw],
                                           *Mzp1, mapW1, mapP1,
                                           mapPhat, comm, ApType) );
        }
        if (Spp == Teuchos::null)
        {
            Spp = Teuchos::rcp(new SppDAMatrix(*Mzp1,*Mzp2,*Auv,*SubMatrix[_Guv],
                                               *SubMatrix[_Duv], comm));
        }
        else
        {
            //note: Guv and Duv are constant and so we can keep them
            Teuchos::rcp_dynamic_cast<SppDAMatrix>(Spp)->Update(*Auv);
        }

        bool rho_mixing = true; // TODO: for the moment this is hard-coded here
        bool rhomu = lsParams.get("ATS: rho/mu Transform", rho_mixing);
        if (rhomu) this->setup_rhomu();
        if (verbose>5)
        {
            INFO("*** Construct Krylov solvers...");
        }

        Teuchos::ParameterList& AuvSolverList = lsParams.sublist("Auv Solver");

        if (AuvSolver==Teuchos::null)
        {
            AuvSolver=SolverFactory::CreateKrylovSolver(AuvSolverList,verbose);
        }
        if (AuvSolver!=Teuchos::null) //may want to use only a prec (option "Method"="None")
        {
            CHECK_ZERO(AuvSolver->SetUserMatrix(Auv.get()));
            //CHECK_ZERO(AuvSolver->SetUserOperator(Auv.get()));
        }

#ifdef DEBUGGING
        comm->Barrier();
#endif

        if (SppSolver==Teuchos::null)
        {
            Teuchos::ParameterList& SppSolverList=lsParams.sublist("Saddlepoint Solver");
            SppSolver=SolverFactory::CreateKrylovSolver(SppSolverList,verbose);
        }
        if (SppSolver!=Teuchos::null)
        {
            CHECK_ZERO(SppSolver->SetUserOperator(Spp.get()));
        }

        if (verbose>5)
        {
            INFO("*** Create Preconditioners...");
        }
        Teuchos::ParameterList& AuvPrecList = lsParams.sublist("Auv Precond");

// Currently the Auv Preconditioner has to be reconstructed because
// the Auv pointer is no longer valid (it is a new one because of the
// call to Utils::RemoveColMap(...))
        {
            DEBUG("Create Auv Preconditioner...");
            AuvPrecond = SolverFactory::CreateAlgebraicPrecond(*Auv,AuvPrecList,verbose);

            DEBUG("Compute Auv Preconditioner...");
            SolverFactory::ComputeAlgebraicPrecond(AuvPrecond,AuvPrecList);

        }


        // currently has to be reconstructed, I think
        {
            DEBUG("Create SppSimplePrecond...");
            Teuchos::ParameterList& SimpleList = lsParams.sublist("Saddlepoint Preconditioner");
            // also pass on info about Auv solver (like tol, maxit etc)
            SimpleList.sublist("Auv Solver")  = lsParams.sublist("Auv Solver");
            SimpleList.sublist("Auv Precond") = lsParams.sublist("Auv Precond");

            // note: the parameters in "Simple: Auv Precond" are ignored as we already have
            // a preconditioner for Auv (and likewise for "Simple: Auv Solver")
            SppPrecond = Teuchos::rcp(new SppSimplePrec(Spp,SimpleList,comm,
                                                        AuvSolver,AuvPrecond, true) );

            bool test_spp = SimpleList.get("Analyze Preconditioned Spectrum",false);
            if (test_spp)
            {
                Teuchos::ParameterList testList;
                Teuchos::RCP<const Epetra_Operator> op = Spp;
                testList.set("Eigen-Analysis: Operator",op);
                SolverFactory::AnalyzeSpectrum(testList,SppPrecond);
            }

        }

        if (ATSSolver==Teuchos::null)
        {
            Teuchos::ParameterList& solverlist = lsParams.sublist("ATS Solver");
            ATSSolver = SolverFactory::CreateKrylovSolver(solverlist,verbose);
        }

        if (ATSSolver!=Teuchos::null)
        {
            if (rhomu)
            {
                CHECK_ZERO(ATSSolver->SetUserMatrix(Arhomu.get()));
            }
            else
            {
                CHECK_ZERO(ATSSolver->SetUserMatrix(ATS.get()));
            }
        }

        // ATS Precond has to be rebuilt. See comment for Auv Precond.
        {

            if (rhomu)
            {
                DEBUG("Create Preconditioner for Arhomu");
                ATSPrecond =
                    SolverFactory::CreateAlgebraicPrecond(*Arhomu,lsParams.sublist("ATS Precond"));
            }
            else
            {
                DEBUG("Create Preconditioner for ATS");
                ATSPrecond =
                    SolverFactory::CreateAlgebraicPrecond(*ATS,lsParams.sublist("ATS Precond"));
            }
            DEBUG("Compute ATSPrecond...");
            SolverFactory::ComputeAlgebraicPrecond(ATSPrecond,lsParams.sublist("ATS Precond"));
        }

        // tell the solvers which preconditioners to use:
        DEBUG("Set Preconditioner Operators...");
        if (AuvSolver!=Teuchos::null)
        {
            // set outer stream for Aztec solver
            AuvSolver->SetOutputStream(*OuterStream);
            AuvSolver->SetErrorStream(*OuterErrorStream);
            CHECK_ZERO(AuvSolver->SetPrecOperator(AuvPrecond.get()));
        }

        if (SppSolver!=Teuchos::null)
        {
            // set outer stream for Aztec solver
            SppSolver->SetOutputStream(*OuterStream);
            SppSolver->SetErrorStream(*OuterErrorStream);
            CHECK_ZERO(SppSolver->SetPrecOperator(SppPrecond.get()));
        }

        if (ATSSolver!=Teuchos::null)
        {
            // set outer stream for Aztec solver
            ATSSolver->SetOutputStream(*OuterStream);
            ATSSolver->SetErrorStream(*OuterErrorStream);
            if (lsParams.sublist("ATS Precond").get("Method","None")!="None")
            {
                CHECK_ZERO(ATSSolver->SetPrecOperator(ATSPrecond.get()));
            }
        }

        DEBUG("leave build_preconditioner");
    }//build_preconditioner


///////////////////////////////////////////////////////////////////////////////
// Apply preconditioner matrix (not available)
///////////////////////////////////////////////////////////////////////////////

    int BlockPreconditioner::
    Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
        // Not implemented: Throw an error!
        ERROR("BlockPreconditioner::Apply() not implemented",__FILE__,__LINE__);
        return -1;
    }

///////////////////////////////////////////////////////////////////////////////
// apply inverse preconditioner step
///////////////////////////////////////////////////////////////////////////////

    int BlockPreconditioner::
    ApplyInverse(const Epetra_MultiVector& input,
                 Epetra_MultiVector& result) const
    {
        bool noisy = (verbose>=8);
        if (noisy) INFO("Apply Block-Preconditioner...");

//  DEBVAR(input);

        if ((input.NumVectors()>1) || (input.NumVectors()!=result.NumVectors()))
        {
            ERROR("Ocean Preconditioner not implemented for multiple RHS, yet!",__FILE__,__LINE__);
        }

// check if input vectors are multivectors or standard vectors
        const Epetra_Vector *input_vec = dynamic_cast<const Epetra_Vector*>(&input);
        Epetra_Vector *result_vec = dynamic_cast<Epetra_Vector*>(&result);

        if (input_vec==NULL)
        {
            input_vec = input(0);
            result_vec= result(0);
        }

        // cast blockvectors into vectors b=input, x=outpu
        const Epetra_Vector& b = *(input_vec);
        Epetra_Vector& x       = *(result_vec);

        // make the solvers report to our own files
        // (note that Aztec uses a static stream
        // because it is based on old C code)

        // note: it is not 100% clear if we have to
        // set it for all solvers, so we just do it
        if (AuvSolver!=Teuchos::null)
        {
            AuvSolver->SetOutputStream(*InnerStream);
            AuvSolver->SetErrorStream(*InnerErrorStream);
        }
        if (SppSolver!=Teuchos::null)
        {
            SppSolver->SetOutputStream(*InnerStream);
            SppSolver->SetErrorStream(*InnerErrorStream);
        }
        if (ATSSolver!=Teuchos::null)
        {
            ATSSolver->SetOutputStream(*InnerStream);
            ATSSolver->SetErrorStream(*InnerErrorStream);
        }

        // (0) PREPROCESSING

        if (noisy)  INFO("(0) Split rhs vector ...");

        // split b = [buv,bw,bp,bTS]' and x = [xuv,xw,xp,xTS]'  // ++scales++
        Epetra_Vector buv(*mapUV);
        Epetra_Vector bw(*mapW1);
        Epetra_Vector bp(*mapP1);
        Epetra_Vector bTS(*mapTS);

        Epetra_Vector xuv(*mapUV);
        Epetra_Vector xw(*mapW1);
        Epetra_Vector xp(*mapP1);
        Epetra_Vector xTS(*mapTS);

        CHECK_ZERO(buv.Export(b,*importUV,Zero));
        CHECK_ZERO(bw.Export(b,*importW1,Zero));
        CHECK_ZERO(bp.Export(b,*importP1,Zero));
        CHECK_ZERO(bTS.Export(b,*importTS,Zero));

        CHECK_ZERO(xuv.Export(x,*importUV,Zero));
        CHECK_ZERO(xw.Export(x,*importW1,Zero));
        CHECK_ZERO(xp.Export(x,*importP1,Zero));
        CHECK_ZERO(xTS.Export(x,*importTS,Zero));

        // set bp = -bp (the sign of the cont. eqn. has been changed)
        CHECK_ZERO(bp.Scale(-1.0));

        Epetra_Vector yuv(*mapUV);
        Epetra_Vector yw(*mapW1);
        Epetra_Vector yp(*mapP1);
        Epetra_Vector yTS(*mapTS);


        // We try to include the buoyancy based on x_init. Apparantly,
        // based on the number of max iterations, x may contain bad
        // stuff... So let's switch this off.
        if (false) 
        {
            CHECK_ZERO(SubMatrix[_BwTS]->Multiply(false,xTS,yw));
            CHECK_ZERO(bw.Update(-1.0,yw,1.0));
            CHECK_ZERO(yw.PutScalar(0.0));
        }

        if (zero_init)
        {
            CHECK_ZERO(xuv.PutScalar(0.0));
            CHECK_ZERO(xw.PutScalar(0.0));
            CHECK_ZERO(xp.PutScalar(0.0));
            CHECK_ZERO(xTS.PutScalar(0.0));
        }

        if (scheme=="ILU")
        {
            ERROR("Block-ILU is no longer supported!!!",__FILE__,__LINE__);
            //solve y = L\b
            /*
              if (noisy) INFO("(1) Solve Ly=b...");
              SolveLower(buv,bw,bp,bTS,yuv, yw,yp,yTS);
            */
            // solve x = U\y
            /*
              if (noisy) INFO("(2) Solve Ux=y...");
              SolveUpper(yuv,yw,yp,yTS,xuv, xw,xp,xTS);
            */
        }
        else if (scheme=="Gauss-Seidel")
        {
            //solve y = (D+wL)\b
            if (noisy) INFO("(1) Solve (D+wL)x=b...");
            if (permutation==1)
                SolveLower1(buv,bw,bp,bTS,xuv, xw,xp,xTS);
            else if (permutation==2)
                SolveLower2(buv,bw,bp,bTS,xuv, xw,xp,xTS);
            else if (permutation==3)
                SolveLower3(buv,bw,bp,bTS,xuv, xw,xp,xTS);
            else
                ERROR("Invalid choice of parameter 'Permutation' (should be 1, 2 or 3)",__FILE__,__LINE__);
        }
        else if (scheme=="symmetric Gauss-Seidel")
        {
            ERROR("symmetric GS is no longer supported!!!",__FILE__,__LINE__);
            /* this method has not proved useful and is no longer supported
            //solve x = (D+wL)\b
            if (noisy) INFO("(1) Solve (D+wL)x=b...");
            SolveLower(buv,bw,bp,bTS,xuv, xw,xp,xTS);
            if (noisy) INFO("(2) apply (D+wU)\\D (BwTS correction)...");
            CHECK_ZERO(SubMatrix[_BwTS]->Multiply(false,xTS,yw));
            yp.PutScalar(0.0);
            Ap->ApplyInverse(yw,yp);
            CHECK_ZERO(xp.Update(-DampingFactor,yp,1.0));
            if (noisy) INFO("(4) Scale with w(2-w)...");
            double fac = DampingFactor*(2-DampingFactor);
            CHECK_ZERO(xuv.Scale(fac));
            CHECK_ZERO(xw.Scale(fac));
            CHECK_ZERO(xp.Scale(fac));
            CHECK_ZERO(xTS.Scale(fac));
            */
        }
        else
        {
            ERROR("Unsupported Scheme: "+scheme,__FILE__,__LINE__);
        }

        // (3) Postprocessing

        if (noisy) INFO("(3) Postprocess: construct final result...");

        // (3.1) fill result vector x
        CHECK_ZERO(x.PutScalar(0.0));
        CHECK_ZERO(x.Import(xuv,*importUV,Add));
        CHECK_ZERO(x.Import(xw,*importW1,Add));
        CHECK_ZERO(x.Import(xp,*importP1,Add));
        CHECK_ZERO(x.Import(xTS,*importTS,Add));

        // reset the static Aztec stream for the outer iteration
        // as it may happen that there is no Auv Solver, we do
        // it for all three solvers to make sure the stream is
        // reset
        if (AuvSolver!=Teuchos::null)
        {
            AuvSolver->SetOutputStream(*OuterStream);
            AuvSolver->SetErrorStream(*OuterErrorStream);
        }
        else if (SppSolver!=Teuchos::null)
        {
            SppSolver->SetOutputStream(*OuterStream);
            SppSolver->SetErrorStream(*OuterErrorStream);
        }
        else if (ATSSolver!=Teuchos::null)
        {
            ATSSolver->SetOutputStream(*OuterStream);
            ATSSolver->SetErrorStream(*OuterErrorStream);
        }
        return 0;
    }

///////////////////////////////////////////////////////////////////////////////
// compute orthogonal transform QTS so that QTS*ATS*QTS is easier to solve   //
///////////////////////////////////////////////////////////////////////////////

    void BlockPreconditioner::setup_rhomu()
    {
        if (QTS==Teuchos::null)
        {
            QTS=Teuchos::rcp(new Epetra_CrsMatrix(Copy,*mapTS,*mapTS,2,true));
            int indices[2];
            double values[2];
            double alphaT = 1.8e-4; //TODO: hard-coded for the moment
            double alphaS = 7.6e-4;
            double lambda = alphaS/alphaT; // par(LAMB) in THCM
            //double idet = 1.0/sqrt(1+lambda*lambda);
            double idet = 1.0/sqrt(2.0);
            // we choose Q such that it is orthonormal, Q^2=I
            int imax = QTS->NumMyRows();

            for (int i=0;i<imax;i+=2)
            {
                int gidT=QTS->GRID(i);
                int gidS=QTS->GRID(i+1);
                indices[0]=gidT; indices[1]=gidS;
                values[0] =-idet; values[1]=lambda*idet;
                CHECK_ZERO(QTS->InsertGlobalValues(gidT,2,values,indices));
                values[0] = idet/lambda; values[1]=idet;
                CHECK_ZERO(QTS->InsertGlobalValues(gidS,2,values,indices));
            }

            for (int i=imax;i<QTS->NumMyRows();i+=2)
            {
                int gidT=QTS->GRID(i);
                int gidS=QTS->GRID(i+1);

                indices[0]=gidT;
                values[0] =1.0;
                CHECK_ZERO(QTS->InsertGlobalValues(gidT,1,values,indices));
                indices[0]=gidS;
                CHECK_ZERO(QTS->InsertGlobalValues(gidS,1,values,indices));
            }

            CHECK_ZERO(QTS->FillComplete());
        }

        Arhomu = Utils::TripleProduct(false,*QTS,false,*ATS,false,*QTS);

        Arhomu->SetLabel("A_(rho,mu)");
#ifdef LINEAR_ARHOMU_MAPS
        //{ use linear map for Arhomu (to allow use of SuperLU)
        Arhomu_linearmap = Teuchos::rcp(new Epetra_Map(Arhomu->NumGlobalRows(),0,Arhomu->Comm()));
        Teuchos::RCP<Epetra_Map> lincolmap = Utils::AllGather(*Arhomu_linearmap);

        int maxlen = Arhomu->MaxNumEntries();
        int len;
        int *ind = new int[maxlen];
        double *val = new double[maxlen];
        int nloc = Arhomu->NumMyRows();
        int *row_lengths = new int[nloc];
        for (int i=0;i<nloc;i++) row_lengths[i]=Arhomu->NumMyEntries(i);
        Teuchos::RCP<Epetra_CrsMatrix> tmpmat;
        tmpmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*Arhomu_linearmap,*lincolmap,row_lengths) );

        int rowA,rowNew;

        for (int i=0;i<Arhomu->NumMyRows();i++)
        {
            rowA = Arhomu->GRID(i);
            rowNew = Arhomu_linearmap->GID(i);
            CHECK_ZERO(Arhomu->ExtractGlobalRowCopy(rowA,maxlen,len,val,ind));
            for (int j=0;j<len;j++)
            {
                int rem=MOD(ind[j],_NUN_);
                int newind=(ind[j]-rem)/_NUN_ + MOD(ind[j],2);
                ind[j] = newind;
            }
            CHECK_ZERO(tmpmat->InsertGlobalValues(rowNew, len, val, ind));
        }
        //Arhomu=Utils::RemoveColMap(tmpmat);
        Arhomu=tmpmat;
        CHECK_ZERO(Arhomu->FillComplete());
#endif
#ifdef STORE_MATRICES
        Utils::Dump(*QTS,"QTS");
        Utils::Dump(*Arhomu,"Arhomu");
#endif
    }

    ///////////////////////////////////////////////////////////////////////////////
    // inf-norm (n/a)
    ///////////////////////////////////////////////////////////////////////////////


    double BlockPreconditioner::NormInf() const
    {
        // Not implemented: Throw an error!
        std::cout << "ERROR: BlockPreconditioner::NormInf() - "
                  << "method is NOT implemented!!  " << std::endl;
        throw "Error: method not implemented";
    }




    //////////////////////////////////////////////////////////////////////////////
    // solve Ly = b for y:                                                      //
    //////////////////////////////////////////////////////////////////////////////
    void BlockPreconditioner::SolveLower1(const Epetra_Vector& buv,
                                          const Epetra_Vector& bw,
                                          const Epetra_Vector& bp,
                                          const Epetra_Vector& bTS,
                                          Epetra_Vector& yuv,
                                          Epetra_Vector& yw,
                                          Epetra_Vector& yp,
                                          Epetra_Vector& yTS) const
    {
#ifdef DUMMY_PREC
        if (DoPresCorr)
        {
            yuv=buv;
            yw=bw;
            yp=bp;
            yTS=bTS;
            double fac1,fac2;
            CHECK_ZERO(yp.Dot(*svp1,&fac1));
            CHECK_ZERO(yp.Dot(*svp2,&fac2));
            CHECK_ZERO(yp.Update(-fac1,*svp1,-fac2,*svp2,1.0));
        }
#else

        // Compute the pressure (yp)
        // Compute ytilp = Ap\[bw,0]'
        Epetra_Vector ytilp(*mapP1);
        Ap->ApplyInverse(bw,ytilp);

        TIMER_START("BlockPrec: solve depth-av Spp");
        // Solve the depth-averaged Saddlepoint problem
        // (a) depth-average bzp = Mzp*bp
        Epetra_Vector bzp(*mapPbar);
        CHECK_ZERO(Mzp2->Multiply(false,bp,bzp));

        // (b) construct 'uv' rhs for Spp

        // yuv = buv-Guv*ytilp
        CHECK_ZERO(SubMatrix[_Guv]->Multiply(false,ytilp,yuv));
        CHECK_ZERO(yuv.Update(1.0,buv,-DampingFactor));
        // (c) construct vector bzuvp = [bzuv,bzp]'
        //     or [buv,bzp]', respectively
        Epetra_Vector bzuvp(Spp->OperatorRangeMap());
        Epetra_Vector yzuvp(Spp->OperatorDomainMap());

        int nzp = bzp.MyLength();

        Teuchos::RCP<Epetra_Vector> bzuv;
        bzuv = Teuchos::rcp(&yuv,false);

        int nzuv = bzuv->MyLength();
        for (int i=0;i<nzuv;i++) bzuvp[i] = (*bzuv)[i];
        for (int i=0;i<nzp ;i++) bzuvp[nzuv+i] = bzp[i];

        yzuvp = bzuvp;

        if (zero_init)
            CHECK_ZERO(yzuvp.PutScalar(0.0));

        { 
            if (SppSolver!=Teuchos::null)
            {
                // (d) solve Saddlepoint problem yzuvp = Spp\bzuvp
                //     using Krylov method
                //     with our own preconditioner
                CHECK_ZERO(SppSolver->SetRHS(&bzuvp));
                CHECK_ZERO(SppSolver->SetLHS(&yzuvp));
                CHECK_NONNEG(SppSolver->Iterate(nitSpp,tolSpp));
            }
            else
                CHECK_ZERO(SppPrecond->ApplyInverse(bzuvp,yzuvp));
        }
        TIMER_STOP("BlockPrec: solve depth-av Spp");

        // Construct the pressure
        // a) yp = ytilp + Mzp1'*yzp
        Epetra_Vector yzp(*mapPbar);
        for (int i=0; i<nzp; i++)
        {
            yzp[i]=yzuvp[nzuv+i];
        }
        CHECK_ZERO(Mzp1->Multiply(true,yzp,yp));
        CHECK_ZERO(yp.Update(1.0,ytilp,1.0));

        // (b)  pressure correction: xp = xp - <xp,svp1>*svp1
        //                                   - <xp,svp2>*svp2
        if (DoPresCorr)
        {
            double fac1,fac2;
            CHECK_ZERO(yp.Dot(*svp1,&fac1));
            CHECK_ZERO(yp.Dot(*svp2,&fac2));
            CHECK_ZERO(yp.Update(-fac1,*svp1,-fac2,*svp2,1.0));
        }
        // Solve the velocity field yuv
        for (int i=0;i<nzuv;i++) yuv[i] = yzuvp[i];

        // Solve vertical velocity field
        // yw = bp(1:nw) - Duv1*yuv
        // note that the sign of Duv has been changed!
        // Duv maps to P1, but Duv1 has been shifted
        // and cropped so that it maps to W1 instead
        CHECK_ZERO(Duv1->Multiply(false,yuv,yw));

        // can't 'Update' because bp lives in the wrong space:
        for (int i=0;i<yw.MyLength();i++) yw[i]=bp[i]-DampingFactor*yw[i];

        // yw = Aw\yw (lower tri-solve)
        Epetra_Vector rhsw = yw;

        // taking care of a no diagonal case
        bool unitDiag = (Aw->NoDiagonal()) ? true : false;

        CHECK_ZERO(Aw->Solve(false, false, unitDiag, rhsw, yw));
        //     Utils::TriSolve(*Aw,rhsw,yw);

        // temperature and salinity equations

        // yTS = BTSuv*yuv
        CHECK_ZERO(SubMatrix[_BTSuv]->Multiply(false,yuv,yTS));

        // yTS2 = BTSw*yw
        Epetra_Vector yTS2 = yTS;
        CHECK_ZERO(SubMatrix[_BTSw]->Multiply(false,yw,yTS2));

        // yTS2 = bTS - yTS - yTS2
        CHECK_ZERO(yTS2.Update(1.0,bTS,-DampingFactor,yTS,-DampingFactor));
        { 
            this->SolveATS(yTS2,yTS,tolATS,nitATS);
        }
#endif

        return;

    } //SolveLower1

    void BlockPreconditioner::SolveLower2(const Epetra_Vector& buv, const Epetra_Vector& bw,
                                          const Epetra_Vector& bp, const Epetra_Vector& bTS,
                                          Epetra_Vector& yuv, Epetra_Vector& yw,
                                          Epetra_Vector& yp, Epetra_Vector& yTS) const
    {

        // Solve the depth-averaged Saddlepoint problem

        // (a) depth-average bzp = Mzp*bp
        Epetra_Vector bzp(*mapPbar);
        CHECK_ZERO(Mzp2->Multiply(false,bp,bzp));

        // (b) construct vector bzuvp = [buv,bzp]'
        Epetra_Vector bzuvp(Spp->OperatorRangeMap());
        Epetra_Vector yzuvp(Spp->OperatorDomainMap());

        int nzp = bzp.MyLength();

        int nuv = buv.MyLength();
        for (int i=0;i<nuv;i++) bzuvp[i] = buv[i];
        for (int i=0;i<nzp ;i++) bzuvp[nuv+i] = bzp[i];

        yzuvp = bzuvp;

        if (zero_init)
        {
            CHECK_ZERO(yzuvp.PutScalar(0.0));
        }
        {

            if (SppSolver!=Teuchos::null)
            {
                // (d) solve Saddlepoint problem yzuvp = Spp\bzuvp using Krylov method
                // with our own preconditioner
                CHECK_ZERO(SppSolver->SetRHS(&bzuvp));
                CHECK_ZERO(SppSolver->SetLHS(&yzuvp));

                CHECK_NONNEG(SppSolver->Iterate(nitSpp,tolSpp));
            }
            else
            {
                CHECK_ZERO(SppPrecond->ApplyInverse(bzuvp,yzuvp));
            }

        }
        // Extract the velocity field yuv
        for (int i=0;i<nuv;i++) yuv[i] = yzuvp[i];

        // Diagnose vertical velocity field from conti-equation

        // yw = bp(1:nw) - Duv1*yuv
        // note that the sign of Duv has been changed!
        // Duv maps to P1, but Duv1 has been shifted
        // and cropped so that it maps to W1 instead
        CHECK_ZERO(Duv1->Multiply(false,yuv,yw));

        // can't 'Update' because bp lives in the wrong space:
        for (int i=0;i<yw.MyLength();i++) yw[i]=bp[i]-DampingFactor*yw[i];

        // yw = Aw\yw (lower tri-solve)
        Epetra_Vector rhsw = yw;
        CHECK_ZERO(Aw->Solve(false,false,false,rhsw,yw));


        // temperature and salinity equantions

        // yTS = BTSuv*yuv
        CHECK_ZERO(SubMatrix[_BTSuv]->Multiply(false,yuv,yTS));

        // yTS2 = BTSw*yw
        Epetra_Vector yTS2 = yTS;
        CHECK_ZERO(SubMatrix[_BTSw]->Multiply(false,yw,yTS2));

        // yTS2 = bTS - yTS - yTS2
        CHECK_ZERO(yTS2.Update(1.0,bTS,-DampingFactor,yTS,-DampingFactor));
        {

            this->SolveATS(yTS2,yTS,tolATS,nitATS);
        }
        // Compute the pressure (yp)

        // a) ytilp = Ap\(bw - BTS*yTS)
        CHECK_ZERO(SubMatrix[_BwTS]->Multiply(false,yTS,rhsw));
        CHECK_ZERO(rhsw.Update(1.0,bw,-1.0));
        Epetra_Vector ytilp(*mapP1);
        Ap->ApplyInverse(rhsw,ytilp);

        Epetra_Vector yzp(*mapPbar);
        for (int i=0; i<nzp; i++)
        {
            yzp[i]=yzuvp[nuv+i];
        }
        CHECK_ZERO(Mzp1->Multiply(true,yzp,yp));
        CHECK_ZERO(yp.Update(1.0,ytilp,1.0));


        // (b)  pressure correction: xp = xp - <xp,svp1>*svp1
        //                                   - <xp,svp2>*svp2
        if (DoPresCorr)
        {
            double fac1,fac2;
            CHECK_ZERO(yp.Dot(*svp1,&fac1));
            CHECK_ZERO(yp.Dot(*svp2,&fac2));
            CHECK_ZERO(yp.Update(-fac1,*svp1,-fac2,*svp2,1.0));
        }

    }//SolveLower2

    void BlockPreconditioner::SolveLower3(const Epetra_Vector& buv, const Epetra_Vector& bw,
                                          const Epetra_Vector& bp, const Epetra_Vector& bTS,
                                          Epetra_Vector& yuv, Epetra_Vector& yw,
                                          Epetra_Vector& yp, Epetra_Vector& yTS) const
    {

        // yw = Aw\bw (lower tri-solve)
        CHECK_ZERO(Aw->Solve(false,false,false,bp,yw));

        // temperature and salinity equantions

        // yTS2 = BTSw*yw
        Epetra_Vector yTS2 = yTS;
        CHECK_ZERO(SubMatrix[_BTSw]->Multiply(false,yw,yTS2));

        // yTS2 = bTS - yTS2
        CHECK_ZERO(yTS2.Update(1.0,bTS,-DampingFactor));
        {

            this->SolveATS(yTS2,yTS,tolATS,nitATS);
        }
        // hydrostatic balance

        // Compute ytilp = Ap\[bw,0]'
        Epetra_Vector rhsw = yw;
        CHECK_ZERO(SubMatrix[_BwTS]->Multiply(false,yTS,rhsw));
        CHECK_ZERO(rhsw.Update(1.0,bw,-1.0));
        Epetra_Vector ytilp(*mapP1);
        CHECK_ZERO(Ap->ApplyInverse(rhsw,ytilp));

        // Saddle point problem

        // (a) depth-average bzp = Mzp*bp
        Epetra_Vector bzp(*mapPbar);
        CHECK_ZERO(Mzp2->Multiply(false,bp,bzp));

        // (b) construct vector bzuvp = [buv-Guv yp,bzp]'
        CHECK_ZERO(SubMatrix[_Guv]->Multiply(false,ytilp,yuv));
        Epetra_Vector bzuvp(Spp->OperatorRangeMap());
        Epetra_Vector yzuvp(Spp->OperatorDomainMap());

        int nzp = bzp.MyLength();
        int nuv = buv.MyLength();

        for (int i=0;i<nuv;i++) bzuvp[i] = buv[i]-yuv[i];
        for (int i=0;i<nzp ;i++) bzuvp[nuv+i] = bzp[i];

        yzuvp = bzuvp;

        if (zero_init)
        {
            CHECK_ZERO(yzuvp.PutScalar(0.0));
        }
        {

            if (SppSolver!=Teuchos::null)
            {
                // (d) solve Saddlepoint problem yzuvp = Spp\bzuvp using Krylov method
                // with our own preconditioner
                CHECK_ZERO(SppSolver->SetRHS(&bzuvp));
                CHECK_ZERO(SppSolver->SetLHS(&yzuvp));

                CHECK_NONNEG(SppSolver->Iterate(nitSpp,tolSpp));
            }
            else
            {
                CHECK_ZERO(SppPrecond->ApplyInverse(bzuvp,yzuvp));
            }

        }
        // Construct the pressure

        // a) yp = ytilp + Mzp1'*yzp
        Epetra_Vector yzp(*mapPbar);
        for (int i=0; i<nzp; i++)
        {
            yzp[i]=yzuvp[nuv+i];
        }
        CHECK_ZERO(Mzp1->Multiply(true,yzp,yp));
        CHECK_ZERO(yp.Update(1.0,ytilp,1.0));

        // (b)  pressure correction: xp = xp - <xp,svp1>*svp1
        //                                   - <xp,svp2>*svp2
        if (DoPresCorr)
        {
            double fac1,fac2;
            CHECK_ZERO(yp.Dot(*svp1,&fac1));
            CHECK_ZERO(yp.Dot(*svp2,&fac2));
            CHECK_ZERO(yp.Update(-fac1,*svp1,-fac2,*svp2,1.0));
        }

    }//SolveLower3

    // apply x=U\y
    void BlockPreconditioner::SolveUpper(const Epetra_Vector& yuv, const Epetra_Vector& yw,
                                         const Epetra_Vector& yp, const Epetra_Vector& yTS,
                                         Epetra_Vector& xuv, Epetra_Vector& xw,
                                         Epetra_Vector& xp, Epetra_Vector& xTS) const

    {
        // temporary vectors
        Epetra_Vector zuv1 = yuv;
        Epetra_Vector zuv = yuv;
        Epetra_Vector zw1 = yw;
        Epetra_Vector zw = yw;
        Epetra_Vector zp = yp;

        // (2) Apply x = U\y
        DEBUG("(3) Solve Ux=y for x");

        // (2.1) compute xTS

        // xTS = yTS
        xTS = yTS;

        // (2.2) compute xw

        // apply zw1 = BwTS*yTS
        CHECK_ZERO(SubMatrix[_BwTS]->Multiply(false,yTS,zw1));

        // apply zp=Ap\(BwTS*yTS)
        Ap->ApplyInverse(zw1,zp);

        //note: zp will be used later on to correct xp

#ifdef TESTING
        {
// check if the Ap solve worked out:
// Ap = [Gw;Mzp]'
            Epetra_Vector vw(*mapW1);
            CHECK_ZERO(SubMatrix[_Gw]->Multiply(false,zp,vw));
            vw.Update(-1.0,zw1,1.0);
            double nrm,nrmb;
            CHECK_ZERO(vw.Norm2(&nrm));
            CHECK_ZERO(zw1.Norm2(&nrmb));
            if (nrm/nrmb>_TESTTOL_)
            {
                INFO("WARNING: ||Ap*(Ap\\zw1)-zw1||_2 = "<<nrm<<"!");
                INFO("        (||zw1||_2 = "<<nrmb<<")");
                INFO("("<<__FILE__<<", line "<<__LINE__<<")");
            }
        }
#endif


        // apply zuv1 = Guv*Ap\BwTS*yTS
        CHECK_ZERO(SubMatrix[_Guv]->Multiply(false,zp,zuv1));

// it is sufficient to apply the preconditioner of Auv once here

        // zuv ~= Auv\zuv1 (apply one preconditioning step)
        CHECK_ZERO(AuvPrecond->ApplyInverse(zuv1,zuv));

        //note: zuv will be used again to compute xuv

        // zw1 = Duv1*zuv
        CHECK_ZERO(Duv1->Multiply(false,zuv,zw1));

        // zw = Aw \ zw1 (lower tri-solve)
        CHECK_ZERO(Aw->Solve(false,false,false,zw1,zw));
//   Utils::TriSolve(*Aw,zw1,zw);

#ifdef TESTING
        {
// check if the Aw solve worked out:
            Epetra_Vector vw(*mapW1);
            CHECK_ZERO(Aw->Multiply(false,zw,vw));
            vw.Update(-1.0,zw1,1.0);
            double nrm,nrmb;
            CHECK_ZERO(vw.Norm2(&nrm));
            CHECK_ZERO(zw1.Norm2(&nrmb));
            if (nrm/nrmb>_TESTTOL_)
            {
                INFO("WARNING: ||Aw*(Aw*zw)-zw1||_2 = "<<nrm<<"!");
                INFO("        (||zw1||_2 = "<<nrmb<<")");
                INFO("("<<__FILE__<<", line "<<__LINE__<<")");
                DEBUG("WARNING: ||Aw*(Aw\\zw1)-zw1||_2 = "<<nrm<<"!");
                DEBUG("        (||zw1||_2 = "<<nrmb<<")");
                DEBVAR(*Aw);
                DEBVAR(zw1);
                DEBVAR(zw);
            }
        }
#endif

        // construct final solution components of x:

        // xw = yw - zw
        CHECK_ZERO(xw.Update(1.0,yw,-1.0,zw,0.0));

        // xuv = yuv + zuv
        CHECK_ZERO(xuv.Update(1.0,yuv,1.0,zuv,0.0));

        // set xp = yp - Ap\BwTS*yTS=yp-zp
        CHECK_ZERO(xp.Update(1.0,yp,-1.0,zp,0.0));

    }//SolveUpper

    void BlockPreconditioner::SolveATS(Epetra_Vector& rhs,
                                       Epetra_Vector& sol,
                                       double tol, int maxit) const
    {
        if (zero_init)
        {
            CHECK_ZERO(sol.PutScalar(0.0));
        }
        Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcp(&rhs,false);
        Teuchos::RCP<Epetra_Vector> sol_ptr = Teuchos::rcp(&sol,false);
        if (QTS!=Teuchos::null)
        {
            rhs_ptr = Teuchos::rcp(new Epetra_Vector(*mapTS));
            sol_ptr = Teuchos::rcp(new Epetra_Vector(*mapTS));
            CHECK_ZERO(QTS->Multiply(false,sol,*sol_ptr));
            CHECK_ZERO(QTS->Multiply(false,rhs,*rhs_ptr));
        }
// TODO: This is for direct solvers for A_(rho/mu) and irrelevant in practice
#ifdef LINEAR_ARHOMU_MAPS
        if (Arhomu_linearmap!=Teuchos::null)
        {
            CHECK_ZERO(rhs_ptr->ReplaceMap(*Arhomu_linearmap));
            CHECK_ZERO(sol_ptr->ReplaceMap(*Arhomu_linearmap));
        }
#endif

        // Solve system with ATS or Arhomu iteratively or apply preconditioner once
        if (ATSSolver!=Teuchos::null)
        {
            TIMER_START("BlockPrec: solve ATS");
            ATSSolver->SetRHS(rhs_ptr.get());
            ATSSolver->SetLHS(sol_ptr.get());
            CHECK_NONNEG(ATSSolver->Iterate(maxit,tol));
            TIMER_STOP("BlockPrec: solve ATS");
        }
        else
        {
            CHECK_ZERO(ATSPrecond->ApplyInverse(*rhs_ptr,*sol_ptr));
        }
#ifdef LINEAR_ARHOMU_MAPS
        if (Arhomu_linearmap!=Teuchos::null)
        {
            CHECK_ZERO(rhs_ptr->ReplaceMap(*Arhomu_linearmap));
            CHECK_ZERO(sol_ptr->ReplaceMap(*Arhomu_linearmap));
        }
#endif
        if (QTS!=Teuchos::null)
        {
            CHECK_ZERO(QTS->Multiply(false,*sol_ptr,sol));
        }
    }

// we need a simple search for column indices since it seems that in parallel
// the notion of 'Sorted()' is different from the serial case (?)
    bool BlockPreconditioner::find_entry(int col, int* indices, int numentries,int& pos)
    {
        pos = 0;
        while (pos<numentries)
        {
            if (indices[pos]==col) break;
            pos++;
        }
        return (pos<numentries);
    }

// store Jacobian, rhs, start guess and all the preconditioner 'hardware'
// (i.e. depth-averaging operators etc) in an HDF5 file
    void BlockPreconditioner::dumpLinSys(const Epetra_Vector& x, const Epetra_Vector& b) const
    {
#ifndef HAVE_XDMF
        INFO("WARNING: cannot dump linear system, hdf5 is not available!");
#else
        Teuchos::RCP<EpetraExt::HDF5> hdf5 = Teuchos::rcp(new EpetraExt::HDF5(*comm));
        hdf5->Create("linsys.h5");
// write linear system:
        hdf5->Write("x",x);
        hdf5->Write("rhs",b);
        hdf5->Write("jacobian",*jacobian);

// write preconditioner hardware:
        hdf5->Write("mzp1",*Mzp1);
        hdf5->Write("mzp2",*Mzp2);
        hdf5->Write("mapuv",*mapUV);
        hdf5->Write("mapw",*mapW1);
        hdf5->Write("mapp",*mapP1);
        hdf5->Write("mappbar",*mapPbar);
        hdf5->Write("mapts",*mapTS);

        hdf5->Write("svp1",*svp1);
        hdf5->Write("svp2",*svp2);

        hdf5->Close();
        comm->Barrier();

// it seems that the HDF5 lib on my notebook doesn't like our singular vectors...
// we store them in ASCII for the moment:
/*
  if (comm->NumProc()>1) ERROR("dump linsys doesn't support parallel runs presently!",__FILE__,__LINE__);
  std::ofstream ofs1("svp1.txt");
  std::ofstream ofs2("svp2.txt");
  ofs1 << *svp1;
  ofs2 << *svp2;
  ofs1.close();
  ofs2.close();
*/
#endif
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION of Ifpack_Preconditioner interface. This allows us to use an                                                 //
// BlockPreconditioner inside an Ifpack_AdditiveSchwarz, i.e. for the Seasonal Cycle problem.                                  //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BlockPreconditioner::BlockPreconditioner(Epetra_RowMatrix* RowMat)
        : label_("Ocean Preconditioner"),
          needs_setup(true), IsComputed_(false)
    {
        INFO("BlockPreconditioner, Ifpack constructor");
        Epetra_CrsMatrix* CrsMat = dynamic_cast<Epetra_CrsMatrix*>(RowMat);
        if (CrsMat==NULL) ERROR("BlockPreconditioner needs CrsMatrix!",__FILE__,__LINE__);
        jacobian = Teuchos::rcp(CrsMat,false);

        // we also need a 'domain' object. The only way to get that right now is
        // from the global THCM instance:
        domain = TRIOS::Static::GetDomain();
        if (domain==Teuchos::null)
        {
            ERROR("could not get static domain pointer, required for Ifpack constructor.",__FILE__,__LINE__);
        }

        comm = domain->GetComm();

        // no params given: set defaults by passing in an empty list
        Teuchos::RCP<Teuchos::ParameterList> List =
            Teuchos::rcp(new Teuchos::ParameterList);
        this->SetParameters(*List);

        // finish constructor
        this->Setup1();
    }

    int BlockPreconditioner::SetParameters(Teuchos::ParameterList &List)
    {
        lsParams = List;
        Teuchos::RCP<std::ostream> defaultStream =
            Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));

        OuterStream      = lsParams.get("Outer Output Stream", defaultStream);
        defaultStream    = Teuchos::rcp_dynamic_cast<Teuchos::FancyOStream>(OuterStream);
        OuterErrorStream = lsParams.get("Outer Error Stream", defaultStream);
        InnerStream      = lsParams.get("Inner Output Stream",defaultStream);
        InnerErrorStream = lsParams.get("Inner Error Stream", defaultStream);

        nitAuv = lsParams.sublist("Auv Solver").get("Max Num Iter",1);
        nitSpp = lsParams.sublist("Saddlepoint Solver").get("Max Num Iter", 5);

        tolAuv = lsParams.sublist("Auv Solver").get("Tolerance",1e-4);
        DEBVAR(tolAuv);
        tolSpp = lsParams.sublist("Saddlepoint Solver").get("Tolerance",1e-8);

        scheme = lsParams.get("Scheme","Gauss-Seidel");
        permutation = lsParams.get("Permutation",1);
        verbose = lsParams.get("Verbosity",10);
        zero_init = lsParams.get("Zero Initial Guess",true);

        DampingFactor = lsParams.get("Relaxation: Damping Factor",1.0);

        // for B-grid
        DoPresCorr = lsParams.get("Subtract Spurious Pressure Modes", true);

        if (scheme=="ILU") //need to solve Schur-complement instead of ATS
        {
            ERROR("BILU Preconditioner is no longer supported!",__FILE__,__LINE__);
        }

        nitATS = lsParams.sublist("ATS Solver").get("Max Num Iter",25);
        tolATS = lsParams.sublist("ATS Solver").get("Tolerance",1e-10);
        return 0;
    }


    int BlockPreconditioner::Initialize()
    {
        // this concept is not clearly implemented here,
        // the ocean preconditioner keeps track of its
        // state using the needs_setup flag
        return 0;
    }

    bool BlockPreconditioner::IsInitialized() const
    {
        return true;
    }

    int BlockPreconditioner::Compute()
    {
        INFO("  Compute Ocean Preconditioner for " << jacobian->Label());

        if (needs_setup) Setup2(); // allocate memory, build submaps...
        // This has to be done exactly once,
        // but not before the Jacobian is there.

        // Extract Submatrices:
        extract_submatrices(*jacobian);

        // construct depth-averaging operators Mzp1/2
        // this has to be done only once as they depend only
        // on the topography and the grid:
        // construct depth-averaging operators Mzp1/2
        // this has to be done only once as they depend only
        // on the topography and the grid:
        if (Mzp1 == Teuchos::null)
        {
            INFO(" build Mzp1...");
            //Mzp1 is the Teuchos::null-space of Gw
            Mzp1 = build_singular_matrix(SubMatrix[_Gw]);
        }

        if (Mzp2 == Teuchos::null)
        {
            //Mzp2 is the Teuchos::null-space of Dw^T
            if (verbose>5)
            {
                INFO("transpose Dw...");
            }
            Epetra_RowMatrixTransposer Trans(SubMatrix[_Dw].get());
            Epetra_CrsMatrix* tmpGw;
            CHECK_ZERO(Trans.CreateTranspose(true, tmpGw, SubMatrixRowMap[_Gw].get()));

            if (verbose>5)
            {
                INFO("build Mzp2...");
            }
            Mzp2 = build_singular_matrix(Teuchos::rcp(tmpGw));
            //TODO: should we or should we not delete this?
            // delete tmpGw;
        }


#ifdef STORE_MATRICES // this is only for debugging, and only for moderate dimensions
        for (int i=0;i<_NUMSUBM;i++)
        {
            Utils::Dump(*SubMatrix[i],SubMatrixLabel[i]);
        }
        Utils::Dump(*Mzp1,"Mzp1");
        Utils::Dump(*Mzp2,"Mzp2");
        Utils::Dump(*Aw,    "Aw");
        Utils::Dump(*Duv1,"Duv1");
#endif

        // build blocksystems, preconditioners and solvers
        build_preconditioner();
        IsComputed_=true;
        return 0;
    }

    bool BlockPreconditioner::IsComputed() const
    {
        // if needs_setup==false, it has certainly been computed once.
        // but not necessarily for this matrix, so I'm not sure here.
        return !needs_setup;
    }

    double BlockPreconditioner::Condest() const
    {
        // we can't compute that right now
        return -1.0;
    }

    double BlockPreconditioner::Condest(const Ifpack_CondestType CT,
                                        const int MaxIters,
                                        const double Tol,
                                        Epetra_RowMatrix* Matrix)
    {
        return -1.0;
    }



    const Epetra_RowMatrix& BlockPreconditioner::Matrix() const
    {
        return *jacobian;
    }

    int BlockPreconditioner::NumInitialize() const
    {
        return -1; // no idea, we don't need this counting business
        // (THCM class takes care of the important counting/timing)
    }

    int BlockPreconditioner::NumCompute() const
    {
        return -1;
    }

    int BlockPreconditioner::NumApplyInverse() const
    {
        return -1;
    }

    double BlockPreconditioner::InitializeTime() const
    {
        return 0.0;
    }
    double BlockPreconditioner::ComputeTime() const
    {
        return 0.0;
    }

    double BlockPreconditioner::ApplyInverseTime() const
    {
        return 0.0;
    }

    double BlockPreconditioner::InitializeFlops() const
    {
        return 0.0;
    }

    double BlockPreconditioner::ComputeFlops() const
    {
        return 0.0;
    }

    double BlockPreconditioner::ApplyInverseFlops() const
    {
        return 0.0;
    }

    std::ostream& BlockPreconditioner::Print(std::ostream& os) const
    {
        os << Label()<<std::endl;
        os << "Scheme: " << scheme << ", ordering: " << permutation << std::endl;
        return os;
    }

#ifdef TESTING
    void BlockPreconditioner::Test()
    {
        Setup2();
        extract_submatrices(*jacobian);
        test_svp();
        Teuchos::RCP<Epetra_CrsMatrix> T =
            Teuchos::rcp(new Epetra_CrsMatrix(Copy, (*SubMatrix[_Guv]).RowMap(),
                                              (*SubMatrix[_Guv]).MaxNumEntries()));
        INFO("Build Mzp1");
        Mzp1 = build_singular_matrix(SubMatrix[_Gw]);

        INFO("Multiply Guv with Mzp1");
        CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(*SubMatrix[_Guv], false,
                                                     *Mzp1, true, *T));
        INFO("transpose Dw...");
        Epetra_RowMatrixTransposer Trans(SubMatrix[_Dw].get());
        Epetra_CrsMatrix* tmpGw;
        CHECK_ZERO(Trans.CreateTranspose(true, tmpGw, SubMatrixRowMap[_Gw].get()));

        INFO("build Mzp2...");
        Mzp2 = build_singular_matrix(Teuchos::rcp(tmpGw));
        Teuchos::RCP<Epetra_CrsMatrix> T2 =
            Teuchos::rcp(new Epetra_CrsMatrix(Copy, (*Mzp2).RowMap(),
                                              (*SubMatrix[_Duv]).ColMap(),
                                              (*Mzp2).MaxNumEntries()));

        DEBVAR(*Mzp2);
        DEBVAR(*SubMatrix[_Duv]);
        DEBVAR((*SubMatrix[_Duv]).Importer());
        DEBVAR((*SubMatrix[_Duv]).ColMap().SameAs(T2->ColMap()));
        INFO("Multiply Mzp2 with Duv");
        CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(*Mzp2, false,
                                                     *SubMatrix[_Duv], false, *T2));
    }
#endif

/////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION of auxiliary operator classes                            //
/////////////////////////////////////////////////////////////////////////////

// This is a matrix class that represents the matrix Ap.
// it needs its own class because we want it to map from
// P to W and still be invertible.

// class ApMatrix

// public:

    // constructor
    ApMatrix::ApMatrix(const Epetra_CrsMatrix &Gw,
                       const Epetra_CrsMatrix &Mp,
                       Teuchos::RCP<Epetra_Map> mapW1_,
                       Teuchos::RCP<Epetra_Map> mapP1_,
                       Teuchos::RCP<Epetra_Map> mapPhat,
                       Teuchos::RCP<Epetra_Comm> comm_,
                       char ApType_ )
        :
        rangeMap(mapW1_), domainMap(mapP1_),
        mapW1(mapW1_), mapP1(mapP1_), ApType(ApType_)
    {
        INFO("ApMatrix constructor...");

        INFO(" building new Ap matrix, Ap = Gw(W1,W1)");
        const Epetra_Map &RowMapGw = Gw.RowMap();
        const Epetra_Map &ColMapGw = Gw.ColMap();
        const Epetra_Map &DomMapGw = Gw.DomainMap();
        const Epetra_Map &RowMapMp = Mp.RowMap();
        const Epetra_Map &ColMapMp = Mp.ColMap();
        const Epetra_Map &DomMapMp = Mp.DomainMap();

        // extract column maps for the blocks of Gw and Mzp
        Teuchos::RCP<Epetra_Map> ColMapGw1 = Teuchos::null;
        Teuchos::RCP<Epetra_Map> DomMapGw1 = Teuchos::null;
        Teuchos::RCP<Epetra_Map> ColMapGw2 = Teuchos::null;
        Teuchos::RCP<Epetra_Map> DomMapGw2 = Teuchos::null;
        Teuchos::RCP<Epetra_Map> ColMapMp1 = Teuchos::null;
        Teuchos::RCP<Epetra_Map> DomMapMp1 = Teuchos::null;
        Teuchos::RCP<Epetra_Map> ColMapMp2 = Teuchos::null;
        Teuchos::RCP<Epetra_Map> DomMapMp2 = Teuchos::null;

        INFO("  Split column maps of Gw...");
        // Gw = [Gw1 Gw2] where Gw1 is square.
        // note: Gw: P1->W1
        // Gw1 contains all P cells from 0 to nrowsGw
        int minGID = 0;
        int maxGID = RowMapGw.MaxAllGID()+(PP-WW);// must adjust from W to P index

        ColMapGw1 = Utils::ExtractRange(ColMapGw,minGID,maxGID);
        DomMapGw1 = Utils::ExtractRange(DomMapGw,minGID,maxGID);
        ColMapMp1 = Utils::ExtractRange(ColMapMp,minGID,maxGID);
        DomMapMp1 = Utils::ExtractRange(DomMapMp,minGID,maxGID);

        // The column map of Gw2 contains the remaining P cells
        minGID    = RowMapGw.MaxAllGID()+(PP-WW)+_NUN_;
        maxGID    = ColMapGw.MaxAllGID();

        ColMapGw2 = Utils::ExtractRange(ColMapGw,minGID,maxGID);
        DomMapGw2 = Utils::ExtractRange(DomMapGw,minGID,maxGID);

        ColMapMp2 = Utils::ExtractRange(ColMapMp,minGID,maxGID);
        DomMapMp2 = Utils::ExtractRange(DomMapMp,minGID,maxGID);

        INFO("  Split matrix...");
        Gw1 = Teuchos::rcp(new Epetra_CrsMatrix(Copy, RowMapGw,     *ColMapGw1, Gw.MaxNumEntries()));
        Gw2 = Teuchos::rcp(new Epetra_CrsMatrix(Copy, RowMapGw,     *ColMapGw2, Gw.MaxNumEntries()));
        Mp1 = Teuchos::rcp(new Epetra_CrsMatrix(Copy, RowMapMp, *ColMapMp1, Mp.MaxNumEntries()));
        Mp2 = Teuchos::rcp(new Epetra_CrsMatrix(Copy, RowMapMp, *ColMapMp2, Mp.MaxNumEntries()));

        // we use dummy importers to do the actual splitting.
        Teuchos::RCP<Epetra_Import> importGw =
            Teuchos::rcp(new Epetra_Import(RowMapGw, RowMapGw) );
        Teuchos::RCP<Epetra_Import> importMp =
            Teuchos::rcp(new Epetra_Import(RowMapMp, RowMapMp) );

        INFO("  Import matrix entries...");
        CHECK_ZERO(Gw1->Import(Gw, *importGw, Zero));
        CHECK_ZERO(Gw2->Import(Gw, *importGw, Zero));
        CHECK_ZERO(Mp1->Import(Mp, *importMp, Zero));
        CHECK_ZERO(Mp2->Import(Mp, *importMp, Zero));

        CHECK_ZERO(Gw1->FillComplete(*DomMapGw1, Gw.RangeMap()));
        CHECK_ZERO(Gw2->FillComplete(*DomMapGw2, Gw.RangeMap()));
        CHECK_ZERO(Mp1->FillComplete(*DomMapMp1, Mp.RangeMap()));
        CHECK_ZERO(Mp2->FillComplete(*DomMapMp2, Mp.RangeMap()));

        Gw1->SetLabel("Gw1"); Gw2->SetLabel("Gw2");
        Mp1->SetLabel("Mp1"); Mp2->SetLabel("Mp2");

#if 0
        DUMPMATLAB("Gw.ascii", Gw);
        DUMPMATLAB("Mp.ascii", Mp);
        DUMPMATLAB("Mp1.ascii", *Mp1);
        DUMPMATLAB("Mp2.ascii", *Mp2);
        DUMPMATLAB("Mp1.ascii", *Mp1);
        DUMPMATLAB("Gw.ascii", Gw);
        DUMPMATLAB("Mp.ascii", Mp);
#endif 

        // create the importers we need in applyinverse
        importPhat = Teuchos::rcp(new Epetra_Import(*mapP1, Mp1->DomainMap()));
        importPbar = Teuchos::rcp(new Epetra_Import(*mapP1, Mp2->DomainMap()));

        INFO("  Testing the splitting...");
        Teuchos::RCP<Epetra_CrsMatrix> R =
            Teuchos::rcp(new Epetra_CrsMatrix(Copy, Gw.RangeMap(), 10) );
        CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(*Gw1, false, *Mp1, true, *R));
        Teuchos::RCP<Epetra_CrsMatrix> S =
            Teuchos::rcp(new Epetra_CrsMatrix(Copy, Gw.RangeMap(), 10) );

        CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(*Gw2, false, *Mp2, true, *S));
        CHECK_ZERO(EpetraExt::MatrixMatrix::Add(*R, false, 1.0, *S, 1.0));
        INFO("   || Gw1*Mp1' + Gw2*Mp2' ||_inf = " << S->NormInf());

        INFO("  replace maps Gw1...");
        Gw1 = Utils::ReplaceBothMaps(Gw1, *mapPhat, *mapPhat);

        CHECK_ZERO(Gw1->FillComplete(*mapP1, *mapP1));

        Gw1->OptimizeStorage();
        Mp1->OptimizeStorage();
        Mp2->OptimizeStorage();

        INFO("ApMatrix constructor: done");
    }

    // destructor
    ApMatrix::~ApMatrix()
    {
        // handled by Teuchos Teuchos::rcp's
    }

    //
    // apply inverse operator x=Ap\b.
    // Ap: P1->W1 => x lives in P1 and b lives in W1.
    // In the cells were w is 'dummy' but p is not, w is assumed to be 0.
    //
    // we have the following block-factorization of Ap:
    //
    //            | I  0  |  | I -M1'| -1  | G1 0 | -1   | I -M1'| |inv(G1) 0|
    //  inv(Ap) = | 0  M2'|  | 0   I |     | M1 I |.   = | 0  M2'| |   M1   I|
    //
    //                                                    inv(AU)   inv(AL)
    //
    // Applying inv(Ap) is done in two distinct steps:
    //
    // 1) lower block triangular solve y = AL\b.     y=[y1,y2]' where y1 lives in P1 and
    //                                               y2 lives in W2.
    //
    // 2) upper block triangular solve z = AU\y.     z=[z1,z2]' where z1 lives in P1 and
    //                                               z2 lives in P2
    //
    // note: alternatively we can just treat Ap as the square part of Gw (Gw1), this approach
    // is now implemented instead
    int ApMatrix::ApplyInverse (const Epetra_Vector &b, Epetra_Vector &x) const
    {

        // DUMP_VECTOR("b.ascii", b);
#ifdef TESTING
        if (!b.Map().SameAs(*rangeMap))
        {
            ERROR("bad rhs vector for solve with Ap!",__FILE__,__LINE__);
        }
        if (!x.Map().SameAs(*domainMap))
        {
            ERROR("bad lhs vector for solve with Ap!",__FILE__,__LINE__);
        }
#endif

        // b is based on the W1 map, x on the P1 map
        // we convert b to a P vector first:
        Epetra_Vector bhat(*mapP1, true);

        for (int i = 0; i < b.MyLength(); i++)
        {
            bhat[i] = b[i];
        }

        // taking care of a no diagonal case
        bool unitDiag = (Gw1->NoDiagonal()) ? true : false;

        if (ApType == 'S') // Only Square part of Gw
        {
            CHECK_ZERO(Gw1->Solve(true, false, unitDiag, bhat, x));
        }
        else if (ApType == 'F') // Full Ap solve
        {
            // Create the support vectors
            Epetra_Vector utmp(Mp1->RangeMap(),  true);
            Epetra_Vector vtmp(Mp2->DomainMap(), true);
            Epetra_Vector wtmp(Mp1->DomainMap(), true);
            Epetra_Vector ztmp(Mp1->DomainMap(), true);

            CHECK_ZERO(Gw1->Solve(true, false, unitDiag, bhat, wtmp));

            CHECK_ZERO(Mp1->Multiply(false, wtmp, utmp));

            CHECK_ZERO(Mp1->Multiply(true, utmp, ztmp));

            CHECK_ZERO(wtmp.Update(-1.0, ztmp, 1.0));

            CHECK_ZERO(Mp2->Multiply(true, utmp, vtmp));
            CHECK_ZERO(vtmp.Scale(-1.0));

            CHECK_ZERO(x.Import(wtmp, *importPhat, Add));
            CHECK_ZERO(x.Import(vtmp, *importPbar, Add));

#if 0
            INFO("  testing ApplyInverse... ");
            Epetra_Vector tmp1(Mp1->RangeMap(), true);
            Epetra_Vector tmp2(Mp2->RangeMap(), true);
            Mp1->Multiply(false, wtmp, tmp1);
            Mp2->Multiply(false, vtmp, tmp2);
            tmp1.Update(1.0, tmp2, 1.0);
            double nrm;
            tmp1.Norm2(&nrm);
            INFO(" ||b2 - (M1*x1 + M2*x2)|| = " << nrm);

            Epetra_Vector tmp3(Gw1->RangeMap(), true);
            Epetra_Vector tmp4(Gw1->RangeMap(), true);

            Gw1->Multiply(false, x, tmp3);
            tmp3.Norm2(&nrm);
            INFO(" ||Gw1*x1|| = " << nrm);
            Gw2->Multiply(false, vtmp, tmp4);
            tmp4.Norm2(&nrm);
            INFO(" ||Gw2*x2|| = " << nrm);

            tmp3.Update(1.0, tmp4, 1.0);
            tmp3.Update(1.0, bhat, -1.0);

            tmp3.Norm2(&nrm);
            INFO(" ||b1 - (G1*x1 + G2*x2)|| = " << nrm);

            wtmp.Norm2(&nrm);
            INFO(" ||x1|| = " << nrm);
            vtmp.Norm2(&nrm);
            INFO(" ||x2|| = " << nrm);

            x.Norm2(&nrm);
            INFO(" ||x|| = " << nrm);
#endif

        }
        return 0;

    }//ApMatrix::ApplyInverse
}//namespace TRIOS
