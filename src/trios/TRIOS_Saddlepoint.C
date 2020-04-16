/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"

#include "AztecOO.h"

#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_VectorOut.h"

#include "Utils.H"

#include "TRIOS_Macros.H"
#include "TRIOS_Saddlepoint.H"
#include "TRIOS_SolverFactory.H"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FancyOStream.hpp"

//////////////////////////////////////////////////////////////////////////////////////////
// Saddlepoint base class. implements construction from blocks, inf norm and MVM
//////////////////////////////////////////////////////////////////////////////////////////

namespace TRIOS {

    //! private constructor
    SaddlepointMatrix::SaddlepointMatrix(Teuchos::RCP<Epetra_Comm> comm)
        :
        label_("Saddlepoint Matrix"),
        comm_(comm)
    {}

    //! constructor
    SaddlepointMatrix::SaddlepointMatrix(Teuchos::RCP<Epetra_CrsMatrix> a11,
                                         Teuchos::RCP<Epetra_CrsMatrix> a12,
                                         Teuchos::RCP<Epetra_CrsMatrix> a21,
                                         Teuchos::RCP<Epetra_Comm> comm)
        :
        label_("Saddlepoint Matrix"),
        comm_(comm)
    {
        SetBlocks(a11,a12,a21);
    }

    // destructor
    SaddlepointMatrix::~SaddlepointMatrix()
    {
        DEBUG("Destroy Saddlepoint Matrix");
    }

    void SaddlepointMatrix::SetBlocks(Teuchos::RCP<Epetra_CrsMatrix> a11,
                                      Teuchos::RCP<Epetra_CrsMatrix> a12,
                                      Teuchos::RCP<Epetra_CrsMatrix> a21)
    {
        A11_ = a11;
        A12_ = a12;
        A21_ = a21;

        // define range and domain of this operator:
        const Epetra_Map& rangemap1  = A11_->RangeMap();
        const Epetra_Map& rangemap2  = A21_->RangeMap();
        const Epetra_Map& domainmap1 = A11_->DomainMap();
        const Epetra_Map& domainmap2 = A12_->DomainMap();
        if (!(rangemap1.SameAs(domainmap1)&&rangemap2.SameAs(domainmap2)))
        {
            INFO("ERROR: bad components for saddlepoint matrix!");
            INFO("range of A11: "<<rangemap1);
            INFO("range of A21: "<<rangemap2);
            INFO("domain of A11: "<<domainmap1);
            INFO("domain of A12: "<<domainmap2);

            ERROR("bad components for saddlepoint matrix!",__FILE__,__LINE__);
        }
        int n1 = rangemap1.NumMyElements();
        int n2 = rangemap2.NumMyElements();

        int *MyGlobalElements = new int[n1+n2];

        // merge the uvbar and Pbar maps
        int offset = rangemap1.MaxAllGID()+1;
        for (int i=0;i<n1; i++) MyGlobalElements[i] = rangemap1.GID(i);
        for (int i=0;i<n2; i++) MyGlobalElements[n1+i] = offset+rangemap2.GID(i);

        rangeMap = Teuchos::rcp(new Epetra_Map(-1,n1+n2,
                                               MyGlobalElements, rangemap1.IndexBase(), *comm_) );

        domainMap = rangeMap;

        // create column maps
        delete [] MyGlobalElements;

        recompute_normInf();

    }
//!@}

//!\name Mathematical functions

    //! apply operator Y=Op*X
    int SaddlepointMatrix::Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
        if (X.NumVectors()>1) ERROR("Only one vector allowed right now!",__FILE__,__LINE__);
        // TODO: only implemented for standard vectors

        // this awkward procedure seems to be required to make it work for both
        // 'real' Epetra_Vectors and Epetra_MultiVectors. I have no clue why.
        // After all, an Epetra_Vector is just a special case of an Epetra_MV.
        const Epetra_Vector* x_ptr = dynamic_cast<const Epetra_Vector*>(&X);
        Epetra_Vector* y_ptr = dynamic_cast<Epetra_Vector*>(&Y);

        if (x_ptr==NULL) x_ptr = X(0);
        if (y_ptr==NULL) y_ptr = Y(0);

        const Epetra_Vector& x = *x_ptr;
        Epetra_Vector& y = *y_ptr;

        const Epetra_Map& map1 = A11_->RowMap();
        const Epetra_Map& map2 = A21_->RowMap();

        // split input and output vectors
        Epetra_Vector x1(map1);
        Epetra_Vector x2(map2);
        Epetra_Vector y1(map1);
        Epetra_Vector y2(map2);

        int n1 = x1.MyLength();
        int n2 = x2.MyLength();

        for (int i=0;i<n1;i++)
        {
            x1[i] = x[i];
        }
        for (int i=0;i<n2;i++)
        {
            x2[i] = x[n1+i];
        }

        this->Apply(x1,x2,y1,y2);

//   CHECK_ZERO(y.Update(1.0,yuv,1.0,yp,1.0)); //(Doesn't work because of yp!)
        for (int i=0;i<n1;i++)
        {
            y[i] = y1[i];
        }
        for (int i=0;i<n2;i++)
        {
            y[n1+i] = y2[i];
        }
        return 0;
    }//Apply


    //! apply operator to pre-split vector
    int SaddlepointMatrix::Apply(const Epetra_Vector& x1, const Epetra_Vector& x2,
                                 Epetra_Vector& y1,       Epetra_Vector& y2) const
    {

        Epetra_Vector tmp1(y1.Map());

        // DEBUG("set y1 = A11*x1...");
        CHECK_ZERO(A11_->Multiply(false,x1,y1));

        // DEBUG("set tmp1 = A12*x2...");
        CHECK_ZERO(A12_->Multiply(false,x2,tmp1));

        // DEBUG("set y2 = A21*x1...");
        CHECK_ZERO(A21_->Multiply(false,x1,y2));

        // construct final y1
        CHECK_ZERO(y1.Update(1.0,tmp1,1.0));
        return 0;
    }//Apply (private)

    //! apply inverse operator (n/a)
    int SaddlepointMatrix::ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
        INFO("WARNING: SaddlepointMatrix::ApplyInverse(...) is not implemented!");
        INFO("("<<__FILE__<<", line " << __LINE__<<")\n");
        return -1;
    }



    //! recompute infty-norm
    void SaddlepointMatrix::recompute_normInf()
    {
        normInf = 0.0;
        double rowsum;
        double *rowA11, *rowA12;
        int *dummyPtr,lenA11,lenA12;
        for (int i=0;i<A11_->NumMyRows();i++)
        {
            CHECK_ZERO(A11_->ExtractMyRowView(i,lenA11,rowA11,dummyPtr));
            CHECK_ZERO(A12_->ExtractMyRowView(i,lenA12,rowA12,dummyPtr));
            rowsum = 0;
            for (int j=0;j<lenA11;j++)
            {
                rowsum += std::abs(rowA11[j]);
            }
            for (int j=0;j<lenA12;j++)
            {
                rowsum += std::abs(rowA12[j]);
            }
            normInf = std::max(normInf,rowsum); // edited_15_10_2014_erik max -> std::max
        }
        // perform maximum over all rows
        double normInfGlobal;
        comm_->MaxAll(&normInf, &normInfGlobal, 1);
        // edited_15_10_2014_erik max -> std::max
        normInf = std::max(normInfGlobal, A21_->NormInf());
    }




////////////////////////////////////////////////////////////////////////////

// class for depth-averaged spp
    SppDAMatrix::SppDAMatrix(Epetra_CrsMatrix& Mzp1, Epetra_CrsMatrix& Mzp2,
                             Epetra_CrsMatrix& Auv,      Epetra_CrsMatrix& Guv,
                             Epetra_CrsMatrix& Duv,
                             Teuchos::RCP<Epetra_Comm> comm)
        : SaddlepointMatrix(comm)
    {
        Teuchos::RCP<Epetra_CrsMatrix> MAuv, MGuv, MDuv;

        const Epetra_Map& mapPbar = Mzp1.RowMap();
        //const Epetra_Map& mapP    = Duv.RowMap();
        const Epetra_Map& mapUV   = Auv.RowMap();

        // note: depth-averaging Auv is no longer implemented
        MAuv = Teuchos::rcp(&Auv,false);
        DEBUG("Depth-average subsystem Guv...");
        MGuv = Utils::MatrixProduct(false,Guv,true,Mzp1);
        DEBUG("Depth-average subsystem Duv...");
#if 1
        MDuv =
            Teuchos::rcp(new Epetra_CrsMatrix(Copy, Mzp2.RowMap(),
                                              Duv.ColMap(),
                                              Mzp2.MaxNumEntries()));
        INFO("Multiply Mzp2 with Duv");
        CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(Mzp2, false,
                                                     Duv, false, *MDuv, false));
#else
        MDuv = Utils::MatrixProduct(false, Mzp2, false, Duv, true);
#endif

        // ___         _
        // Duv:  uv -> p
        INFO("CHECK_ZERO(MDuv->FillComplete(mapUV,mapPbar))");
        CHECK_ZERO(MDuv->FillComplete(mapUV,mapPbar));

        // ___   _
        // Guv:  p -> uv_ZERO(MGuv->FillComplete(mapPbar,mapUV));

        label_ = "Depth-averaged Spp (U/V/p)";

        SaddlepointMatrix::SetBlocks(MAuv,MGuv,MDuv);
    }

    void SppDAMatrix::Update(Epetra_CrsMatrix& Auv)
    {
        A11_= Teuchos::rcp(&Auv,false);
        recompute_normInf();
    }

// destructor
    SppDAMatrix::~SppDAMatrix()
    {
        DEBUG("Destroy SppDAMatrix");
    }

//////////////////////////////////////////////////
// class SppSimplePrec : public Epetra_Operator //
//////////////////////////////////////////////////

// constructor
    SppSimplePrec::SppSimplePrec(Teuchos::RCP<SaddlepointMatrix> Spp_,
                                 Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Comm> comm_,
                                 Teuchos::RCP<AztecOO> A11Solver_,
                                 Teuchos::RCP<Epetra_Operator> A11Precond_,
                                 bool zero_init_)
        :
        comm          (comm_),
        zero_init     (zero_init_),
        Spp           (Spp_),
        A11Solver     (A11Solver_),
        A11Precond    (A11Precond_)
    {
        scheme            = params.get("Scheme","SR");
        scaleChat         = params.get("Scale Chat", false);
        fixSingularChat   = params.get("Fix singular Chat", false);
        fixChatTol        = params.get("Fix Chat tolerance", 1e3);

        printSingularChat = params.get("Print zero diagonal indices", false);
        fixSingularA11    = params.get("Fix singular A11", false);

        Teuchos::ParameterList& SpaIList = params.sublist("Approximate Inverse");
        std::string spai_scheme=SpaIList.get("Method","Block Diagonal");
        label_ = "Simple Preconditioner ("+scheme+", "+spai_scheme+")";

        // verify the scheme is valid
        if (scheme!="SI"&&scheme!="SL"&&scheme!="SR")
        {
            INFO("WARNING: invalid scheme for Simple preconditioner: "<<scheme);
            INFO("         currently supported schemes are SI,SL and SR");
            INFO("         Will use SI scheme.                         ");
            INFO("         ("<<__FILE__<<", line "<<__LINE__<<")\n");
            scheme="SI";
        }

        DEBUG("Build Spp Simple preconditioner...");

        rangeMap  = Spp->GetRangeMap();
        domainMap = Spp->GetDomainMap();

        {

#ifdef HAVE_PARASAILS
            if (spai_scheme=="ParaSails")
            {
                Teuchos::RCP<Epetra_Operator> spai =
                    SolverFactory::CreateAlgebraicPrecond(Spp->A11(),SpaIList);

                SolverFactory::ComputeAlgebraicPrecond(spai,SpaIList);
                Teuchos::RCP<ParaSailsPrecond> psPre =
                    Teuchos::rcp_dynamic_cast<ParaSailsPrecond>(spai);

                if (psPre== Teuchos::null) ERROR("failed to dynamic_cast SpaI",__FILE__,__LINE__);
                DEBUG("extract SpaI CRS matrix...");
                BlockDiagA11 = psPre->GetSpaI();
            }
            else
#endif
            {
                if (spai_scheme!="Block Diagonal")
                {
                    INFO("WARNING: Bad Sparse approximate inverse method: "+spai_scheme);
                    INFO("         using \"Block Jacobi\" instead!");
                }
                // create CRS matrix to store the block diagonal of A11:
                BlockDiagA11 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,Spp->A11().RowMap(),2,true) );

                DEBUG("extract (inverse) block diagonal from Auv...");
                ExtractInverseBlockDiagonal(Spp->A11(),*BlockDiagA11); //??

                BlockDiagA11->SetLabel("(inverse) 2x2 Block-diagonal of A11");
            }
        }
        {

            DEBUG("compute the Schur-complement...");
            DEBUG(" compute Chat = Duv*inv(diag(Auv))*Guv");

            DEBUG("  initialize TMP");
            Teuchos::RCP<Epetra_CrsMatrix> TMP =
                Teuchos::rcp(new Epetra_CrsMatrix(Copy, (Spp->A21()).RowMap(),
                                                  (Spp->A21()).MaxNumEntries()));
            DEBUG("  initialize AB");
            Teuchos::RCP<Epetra_CrsMatrix> AB =
                Teuchos::rcp(new Epetra_CrsMatrix(Copy, (Spp->A21()).RowMap(),
                                                  (Spp->A21()).MaxNumEntries()));



            DEBUG("  perform AB = Spp->A21*BlockDiagA11");
            EpetraExt::MatrixMatrix::Multiply(Spp->A21(),    false,
                                              *BlockDiagA11, false, *AB);
            DEBUG("  perform Chat = AB*Spp->A12");
            EpetraExt::MatrixMatrix::Multiply(*AB, false, Spp->A12(), false, *TMP );
            DEBUG("  finished MM's");
            Chat = TMP;

            /*
              Chat = Utils::TripleProduct(false,  Spp->A21(),
              false, *BlockDiagA11,
              false,  Spp->A12());
            */

            CHECK_ZERO(Chat->Scale(-1.0));
            Chat->SetLabel("Schur-Complement Chat of Simple Precond");
            // if AdjustChat does something to the local part of Chat,
            // these indicate which rows have been modified after the call:
            fixp1 = -1;
            fixp2 = -1;
            valp=0.0;
            this->AdjustChat(Chat);
        }

        /*
          Chat = MatrixUtils::ReadThcmMatrix("Chat",*comm,Chat->RowMap(),Chat->ColMap());
          BlockDiagA11 = MatrixUtils::ReadThcmMatrix("PDAuv",*comm,BlockDiagA11->RowMap(),BlockDiagA11->ColMap());
        */

#ifdef DEBUGGING
        if (comm->MyPID()==0) std::cout << " original Chat: "<<std::endl;
        comm->Barrier();
        std::cout << "PID "<<comm->MyPID()<<": Chat rows="<<Chat->NumMyRows()<<std::endl;
#endif
#ifdef HAVE_ZOLTAN
        bool repart = params.get("Repartition Chat",false);
        if (repart)
        {
            RepartChat = Teuchos::rcp(new Repart(*Chat));
            Chat = RepartChat->Redistribute(*Chat);
        }
#endif

        // -->unused parameters...?
        //Teuchos::ParameterList& AuvSolverList = params.sublist("Auv Solver");
        //Teuchos::ParameterList& AuvPrecList = params.sublist("Auv Precond");

        // note: the solver may be null if "None" is specified as Method
        if (A11Solver.get()!=NULL)
        {
            CHECK_ZERO(A11Solver->SetUserMatrix(&Spp->A11()));
            CHECK_ZERO(A11Solver->SetPrecOperator(A11Precond.get()));
        }

        Teuchos::ParameterList& ChatSolverList = params.sublist("Chat Solver");
        Teuchos::ParameterList& ChatPrecList   = params.sublist("Chat Precond");

        // make preconditioner for Chat
        DEBUG("Construct preconditioner for Chat...");
        {
            ChatPrecond = SolverFactory::CreateAlgebraicPrecond(*Chat, ChatPrecList);
            DEBUG("Compute preconditioner for Chat...");
            SolverFactory::ComputeAlgebraicPrecond(ChatPrecond, ChatPrecList);
        }

        ChatSolver = SolverFactory::CreateKrylovSolver(ChatSolverList);
        if (ChatSolver.get()!=NULL)
        {
            CHECK_ZERO(ChatSolver->SetUserMatrix(Chat.get()));
            std::string cHat_prec = ChatPrecList.get("Method","None");
            if (cHat_prec!="None")
            {
                CHECK_ZERO(ChatSolver->SetPrecOperator(ChatPrecond.get()));
            }
        }

#ifdef DEBUGGING
        if (comm->MyPID()==0) std::cout << " redistributed Chat: "<<std::endl;
        comm->Barrier();
        std::cout << "PID "<<comm->MyPID()<<": Chat rows="<<Chat->NumMyRows()<<std::endl;
#endif

        Teuchos::ParameterList& A11SolverList = params.sublist("Auv Solver");
        //Teuchos::ParameterList& A11PrecList = params.sublist("Auv Precond");

        nitA11 = A11SolverList.get("Max Num Iter",10);
        tolA11 = A11SolverList.get("Tolerance",1.0e-8);
        nitChat = ChatSolverList.get("Max Num Iter",10);
        tolChat = ChatSolverList.get("Tolerance",1.0e-8);

    }//SppSimplePrec


// Destructor
    SppSimplePrec::~SppSimplePrec()
    {
        DEBUG("Destroy SppSimplePrec");
        // handled by Teuchos Teuchos::rcp's
    }

// fix two points of the pressure to avoid singular Chat
    void SppSimplePrec::AdjustChat(Teuchos::RCP<Epetra_CrsMatrix> P)
    {
        int maxgid = Chat->Map().MaxAllGID();
        //int maxgid = Chat->Map().MaxMyGID();
        if (Chat->MyGRID(maxgid))
        {
            fixp1=Chat->LRID(maxgid-1);
            fixp2=Chat->LRID(maxgid);
            int maxlen=Chat->MaxNumEntries();
            int *indices = new int[maxlen];
            double *values = new double[maxlen];
            int len;
            int row = maxgid-1;
            CHECK_ZERO(Chat->ExtractGlobalRowCopy(row,maxlen,len,values,indices));
            for (int i=0;i<len;i++)
            {
                if (indices[i]==row)
                {
                    values[i]=1.0;
                }
                else
                {
                    values[i]=0.0;
                }
            }
            CHECK_ZERO(Chat->ReplaceGlobalValues(row,len,values,indices));
            row = maxgid;
            CHECK_ZERO(Chat->ExtractGlobalRowCopy(row,maxlen,len,values,indices));
            for (int i=0;i<len;i++)
            {
                if (indices[i]==row)
                {
                    values[i]=1.0;
                }
                else
                {
                    values[i]=0.0;
                }
            }
            CHECK_ZERO(Chat->ReplaceGlobalValues(row,len,values,indices));


            delete [] indices;
            delete [] values;
        }

        if (printSingularChat)
        {
            Epetra_Vector diagonal(Chat->RowMap());
            Epetra_Vector zeroRows(Chat->RowMap());
            CHECK_ZERO(Chat->ExtractDiagonalCopy(diagonal));
            INFO("  chat diagonal length = " << diagonal.GlobalLength());
            int numMyElements = diagonal.Map().NumMyElements();
            double tol = 1e-10;

            for (int i = 0; i != numMyElements; ++i)
                if (std::abs(diagonal[i]) < tol)
                    zeroRows[i] = 1;

            INFO("Printing zero rows in Chat");
            EpetraExt::VectorToMatlabFile("zerorows", zeroRows);
        }


        // =============================================================================
        // IMPROVE CONDITION NUMBER OF CHAT ----
        // I'm also going to try to fix zero diagonal elements due to the landmask -Erik
        if (fixSingularChat)
        {
            Epetra_Vector diagonal(Chat->RowMap());
            CHECK_ZERO(Chat->ExtractDiagonalCopy(diagonal));
            INFO("  Fixing singular Chat...");
            INFO("  tolerance = " << fixChatTol);
            INFO("  chat diagonal length = " << diagonal.GlobalLength());

            Teuchos::RCP<Epetra_CrsMatrix> TMP =
                Teuchos::rcp(new Epetra_CrsMatrix(Copy, Chat->RowMap(), 0));

            double values[1]  = {0.0};
            int colinds[1]    = {0};
            int numMyElements = diagonal.Map().NumMyElements();

            int *myGlobalElements = diagonal.Map().MyGlobalElements();
            int row;
            for (int i = 0; i != numMyElements; ++i)
            {
                if (std::abs(diagonal[i]) < fixChatTol)
                {
                    values[0] = (i > 0) ? diagonal[i-1] : -fixChatTol;
                    row = myGlobalElements[i];
                    colinds[0]  = row;
                    diagonal[i] = values[0];
                    CHECK_ZERO(TMP->InsertGlobalValues(row, 1, values, colinds));
                }
            }
            Epetra_Import Chat2TMP(TMP->RowMap(), Chat->RowMap());

            CHECK_ZERO(TMP->Import(*Chat, Chat2TMP, Insert));
            CHECK_ZERO(TMP->FillComplete());

            Chat = TMP;
            INFO("  Fixing singular Chat... done");

        }

        // Traditional scaling although I'm still not sure whether I'm doing it right
        // See also ApplyInverse()
        if (scaleChat)
        {
            scalingChat = Teuchos::rcp(new Epetra_Vector(Chat->RowMap()));

            Chat->InvRowSums(*scalingChat);
            Chat->LeftScale(*scalingChat);
            // Chat->RightScale(*scalingChat);
            // DUMPMATLAB("CHAT3", *Chat);
        }
    }


// extracts 2x2 block diagonal from a matrix, inverts the blocks and stores
// them as a new Crs matrix.
// note: this is _not_ very general
    void SppSimplePrec::ExtractInverseBlockDiagonal(const Epetra_CrsMatrix& A,
                                                    Epetra_CrsMatrix& D)
    {
        int *indAu, *indAv, indDu[2], indDv[2];
        double *valAu, *valAv, valDu[2], valDv[2];
        int len,maxlen;
        int grid;

        maxlen = A.MaxNumEntries();

        indAu = new int[maxlen];
        indAv = new int[maxlen];
        valAu = new double[maxlen];
        valAv = new double[maxlen];

        // this makes the assumption that the matrix entries
        // of A are sorted by ascending col index. This is
        // true as soon as FillComplete() has
        // been called (as of Trilinos 7.0.4)
        if (!A.Filled())
        {
            ERROR("Simple: Matrix A11 not filled!",__FILE__,__LINE__);
        }

        for (int i = 0; i < A.NumMyRows(); i += 2)
        {
            grid = A.GRID(i);
            // 'u' variable
            CHECK_ZERO(A.ExtractGlobalRowCopy(grid,maxlen,len,valAu,indAu) );
#if 0
            DEBUG("A11, global row " << i << ": \n");
            for (int j = 0; j < len; j++)
                DEBUG(grid << "\t" << indAu[j] << "\t" <<valAu[j] << "\n");
#endif
            int j = 0;
            while ((indAu[j]!=grid)&&(j<len)) j++;
            int k = 0;
            while ((indAu[k]!=grid+1)&&(k<len)) k++;
            if (j>=len)
            {
                ERROR("matrix has zero on diagonal!",__FILE__,__LINE__);
            }
            indDu[0]=grid;
            indDu[1]=grid+1;
            valDu[0] = valAu[j];
            valDu[1] = (k<len)? valAu[k]:0.0;

            // 'v' variable
            CHECK_ZERO(A.ExtractGlobalRowCopy(grid+1,maxlen,len,valAv,indAv) );
#ifdef DEBUGGING
            //(*debug) << "A11, global row "<<grid+1<<": \n";
            //for (int j=0;j<len;j++) (*debug) << grid+1<<"\t"<<indAv[j]<<"\t"<<valAv[j]<<"\n";
#endif
            j=0;
            while ((indAv[j]!=grid+1)&&(j<len)) j++;
            k=0;
            while ((indAv[k]!=grid)&&(k<len)) k++;
            if (j>=len)
            {
                ERROR("matrix has zero on diagonal!",__FILE__,__LINE__);
            }
            indDv[0]=grid;
            indDv[1]=grid+1;
            valDv[0] = (k<len)? valAv[k]:0.0;
            valDv[1] = valAv[j];

            // invert the block
            double detInv = 1.0/(valDu[0]*valDv[1]-valDu[1]*valDv[0]);

            double tmp = valDu[0];
            valDu[0]=valDv[1]*detInv; valDv[1]=tmp*detInv;
            valDu[1]*=-detInv; valDv[0]*=-detInv;
#ifdef DEBUGGING
            //(*debug) << "inverse block: \n";
            //(*debug) << "("<<indDu[0]<<") "<<valDu[0]<<"\t("<<indDu[1]<<") "<<valDu[1]<<"\n";
            //(*debug) << "("<<indDv[0]<<") "<<valDv[0]<<"\t("<<indDv[1]<<") "<<valDv[1]<<"\n";
#endif

            // insert the block
            CHECK_ZERO(D.InsertGlobalValues(grid,2,valDu,indDu));
            CHECK_ZERO(D.InsertGlobalValues(grid+1,2,valDv,indDv));
        }
        CHECK_ZERO(D.FillComplete());

        delete [] indAu;
        delete [] indAv;
        delete [] valAu;
        delete [] valAv;

    }//ExtractInverseBlockDiagonal


// Apply preconditioner operator (n/a)
    int SppSimplePrec::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
        INFO("WARNING: SppSimplePrec::Apply not implemented!");
        INFO("("<<__FILE__<< ", line "<<__LINE__<<")");
        return -1;
    }


// Apply preconditioner operator inverse
    int SppSimplePrec::ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
    {
        if (X.NumVectors()>1) ERROR("Only one vector allowed right now!",__FILE__,__LINE__);

        // TODO: only implemented for standard vectors
        const Epetra_Vector *b_ptr = dynamic_cast<const Epetra_Vector*>(&B);
        Epetra_Vector *x_ptr = dynamic_cast<Epetra_Vector*>(&X);


        if (x_ptr==NULL) x_ptr = X(0);
        if (b_ptr==NULL) b_ptr = B(0);

        const Epetra_Vector& b = *b_ptr;
        Epetra_Vector& x = *x_ptr;

        // DEBUG("Apply SppSimplePrec...");

        // temporary vector
        const Epetra_Map& map1 = Spp->A11().RowMap();
        const Epetra_Map& map2 = Spp->A21().RowMap();

        Teuchos::RCP<Epetra_Vector> x1 = Teuchos::rcp(new Epetra_Vector(map1));
        Teuchos::RCP<Epetra_Vector> x2 = Teuchos::rcp(new Epetra_Vector(map2));

        Teuchos::RCP<Epetra_Vector> b1 = Teuchos::rcp(new Epetra_Vector(map1));
        Teuchos::RCP<Epetra_Vector> b2 = Teuchos::rcp(new Epetra_Vector(map2));

        int n1 = b1->MyLength();
        int n2 = b2->MyLength();

        // split vector b = [b1;b2]
        for (int i=0;i<n1;i++) (*b1)[i] = b[i];
        for (int i=0;i<n2;i++) (*b2)[i] = b[n1+i];

        if (scheme=="SI")
        {
            CHECK_ZERO(this->ApplyInverse(*b1,*b2,*x1,*x2,false));
        }
        else if (scheme=="SL")
        {
            CHECK_ZERO(this->ApplyInverse(*b1,*b2,*x1,*x2,true));
        }
        else if (scheme=="SR"||scheme=="SPAI")
        {
            Teuchos::RCP<Epetra_Vector> xtmp1 = Teuchos::rcp(new Epetra_Vector(map1));
            Teuchos::RCP<Epetra_Vector> xtmp2 = Teuchos::rcp(new Epetra_Vector(map2));
            Teuchos::RCP<Epetra_Vector> btmp1 = Teuchos::rcp(new Epetra_Vector(map1));
            Teuchos::RCP<Epetra_Vector> btmp2 = Teuchos::rcp(new Epetra_Vector(map2));
            // apply SL step:
            CHECK_ZERO(this->ApplyInverse(*b1,*b2,*x1,*x2,true));

            // apply saddlepoint operator

            CHECK_ZERO(Spp->Apply(*x1,*x2,*btmp1,*btmp2));
            //CHECK_ZERO(Spp->A11().Multiply(false,x1,xtmp1));
            //CHECK_ZERO(Spp->A12().Multiply(false,x2,btmp1));

            CHECK_ZERO(btmp1->Update(1.0,*b1,-1.0));
            CHECK_ZERO(btmp2->Update(1.0,*b2,-1.0));

            // apply SI step
            CHECK_ZERO(this->ApplyInverse(*btmp1,*btmp2,*xtmp1,*xtmp2,false));

            // construct final result
            CHECK_ZERO(x1->Update(1.0,*xtmp1,1.0));
            CHECK_ZERO(x2->Update(1.0,*xtmp2,1.0));
        }

        // compose final vector x = [xuv;xp]
        for (int i=0;i<n1;i++) x[i] = (*x1)[i];
        for (int i=0;i<n2;i++) x[n1+i] = (*x2)[i];

        return 0;
    }

// apply standard Simple method (SI, if transp=false) or
// simple(L) (SL if transp=true);
    int SppSimplePrec::ApplyInverse(Epetra_Vector& b1, Epetra_Vector& b2,
                                    Epetra_Vector& x1, Epetra_Vector& x2,
                                    bool trans) const
    {
        Teuchos::RCP<Epetra_Vector> y1     =Teuchos::rcp(new Epetra_Vector(b1.Map()));
        Teuchos::RCP<Epetra_Vector> ytmp1  =Teuchos::rcp(new Epetra_Vector(b1.Map()));
        Teuchos::RCP<Epetra_Vector> y2     =Teuchos::rcp(new Epetra_Vector(b2.Map()));
        Teuchos::RCP<Epetra_Vector> ytmp2  =Teuchos::rcp(new Epetra_Vector(b2.Map()));
        Teuchos::RCP<Epetra_Vector> rhs,sol;

        if (!trans) // Simple
        {
            {
                if (zero_init)
                {
                    CHECK_ZERO(y1->PutScalar(0.0));
                }
                // apply inv(L):
                TIMER_START("BlockPrec: solve Auv");
                if (A11Solver.get()==NULL)
                {
                    CHECK_ZERO(A11Precond->ApplyInverse(b1,*y1));
                }
                else
                {
                    CHECK_ZERO(A11Solver->SetRHS(&b1));
                    CHECK_ZERO(A11Solver->SetLHS(y1.get()));
                    CHECK_NONNEG(A11Solver->Iterate(nitA11,tolA11));
                }
                TIMER_STOP("BlockPrec: solve Auv");
            }
            CHECK_ZERO(Spp->A21().Multiply(false,*y1,*y2));
            CHECK_ZERO(y2->Update(1.0,b2,-1.0));
            // fix pressure in two points (if they are on this subdomain)
            if (fixp1>=0) (*y2)[fixp1]=valp;
            if (fixp2>=0) (*y2)[fixp2]=valp;
            {
                if (zero_init)
                {
                    CHECK_ZERO(x2.PutScalar(0.0));
                }

                // apply inv(L):
                rhs=y2; sol=Teuchos::rcp(&x2,false);
                if (scaleChat) rhs->Multiply(1.0, *scalingChat,*rhs, 0.0); // scale rhs
#ifdef HAVE_ZOLTAN
                if (RepartChat!= Teuchos::null)
                {
                    rhs = Teuchos::rcp(new Epetra_Vector(Chat->RowMap()));
                    sol = Teuchos::rcp(new Epetra_Vector(Chat->RowMap()));
                    RepartChat->Redistribute(*y2,*rhs);
                }
#endif
                TIMER_START("BlockPrec: solve Chat");
                if (ChatSolver.get()==NULL)
                {
                    CHECK_ZERO(ChatPrecond->ApplyInverse(*rhs,*sol));
                }
                else
                {
                    CHECK_ZERO(ChatSolver->SetRHS(rhs.get()));
                    CHECK_ZERO(ChatSolver->SetLHS(sol.get()));
                    CHECK_NONNEG(ChatSolver->Iterate(nitChat,tolChat));
                }
                TIMER_STOP("BlockPrec: solve Chat");
#ifdef HAVE_ZOLTAN
                if (RepartChat!= Teuchos::null)
                {
                    RepartChat->Undistribute(*sol,x2);
                }
#endif
                // if (scaleChat) sol->Multiply(1.0, *scalingChat,*sol, 0.0); // scale sol
            }
            CHECK_ZERO(Spp->A12().Multiply(false,x2,*ytmp1));
            CHECK_ZERO(BlockDiagA11->Multiply(false,*ytmp1,x1));
            CHECK_ZERO(x1.Update(1.0,*y1,-1.0));
        }
        else  // Simple(L)
        {
            // apply inv(U')
            CHECK_ZERO(BlockDiagA11->Multiply(false,b1,*y1));
            CHECK_ZERO(Spp->A21().Multiply(false,*y1,*y2));
            CHECK_ZERO(y2->Update(1.0,b2,-1.0));

            if (fixp1>=0) (*y2)[fixp1]=valp;
            if (fixp2>=0) (*y2)[fixp2]=valp;
            {

                if (zero_init)
                {
                    CHECK_ZERO(x2.PutScalar(0.0));
                }
                rhs=y2; sol=Teuchos::rcp(&x2,false);
                if (scaleChat) rhs->Multiply(1.0, *scalingChat,*rhs, 0.0); // scale rhs
#ifdef HAVE_ZOLTAN
                if (RepartChat!= Teuchos::null)
                {
                    sol = Teuchos::rcp(new Epetra_Vector(Chat->RowMap()));
                    rhs = Teuchos::rcp(new Epetra_Vector(Chat->RowMap()));
                    RepartChat->Redistribute(*y2,*rhs);
                }
#endif
                TIMER_START("BlockPrec: solve Chat");
                if (ChatSolver.get()==NULL)
                {
                    CHECK_ZERO(ChatPrecond->ApplyInverse(*rhs,*sol));
                }
                else
                {
                    CHECK_ZERO(ChatSolver->SetRHS(rhs.get()));
                    CHECK_ZERO(ChatSolver->SetLHS(sol.get()));
                    CHECK_NONNEG(ChatSolver->Iterate(nitChat,tolChat));
                }
                TIMER_STOP("BlockPrec: solve Chat");
#ifdef HAVE_ZOLTAN
                if (RepartChat!= Teuchos::null)
                {
                    RepartChat->Undistribute(*sol,x2);
                }
#endif
                // if (scaleChat) sol->Multiply(1.0, *scalingChat,*sol, 0.0); // scale sol
            }

            // apply inv(L')
            CHECK_ZERO(Spp->A12().Multiply(false,x2,*y1));
            CHECK_ZERO(y1->Update(1.0,b1,-1.0));
            {


                if (zero_init)
                {
                    CHECK_ZERO(x1.PutScalar(0.0));
                }

                TIMER_START("BlockPrec: solve Auv");
                if (A11Solver.get()==NULL)
                {
                    CHECK_ZERO(A11Precond->ApplyInverse(*y1,x1));
                }
                else
                {
                    CHECK_ZERO(A11Solver->SetRHS(y1.get()));
                    CHECK_ZERO(A11Solver->SetLHS(&x1));
                    CHECK_NONNEG(A11Solver->Iterate(nitA11,tolA11));
                }
                TIMER_STOP("BlockPrec: solve Auv");
            }
        }
        return 0;
    }

// Computing infinity norm
    double SppSimplePrec::NormInf() const
    {
        INFO("WARNING: SppSimplePrec::NormInf not implemented!");
        INFO("("<<__FILE__<< ", line "<<__LINE__<<")");
        return -1;
    }


}//namespace TRIOS

////////////////////////////////////////////////////////////////



