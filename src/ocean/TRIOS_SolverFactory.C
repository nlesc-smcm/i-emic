/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#include "TRIOS_SolverFactory.H"
#include "Teuchos_Utils.hpp"
#include <sstream>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StrUtils.hpp>
#include "AztecOO.h"
#include "Epetra_Time.h"
#include "AztecOO_string_maps.h"
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_MRILU.h"
#include <iomanip>
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
// ML
#include "ml_MultiLevelPreconditioner.h"

#include "Utils.H"

#include "ml_struct.h"
#include "ml_operator.h"
#include "ml_epetra_utils.h"

// block preconditioner for THCM jacobian
#include "TRIOS_BlockPreconditioner.H"

// for the info stream
#include "GlobalDefinitions.H"

// don't know why...
#ifdef DEBUGGING
#ifdef HAVE_ANASAZI
#undef HAVE_ANASAZI
#endif
#endif

// this is for eigen-analysis using Anasazi BlockKrylovSchur:
#ifdef HAVE_ANASAZI
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#endif

using std::max;
using std::min;

namespace TRIOS {

// create an algebraic preconditioner (i.e. for a submatrix)
    Teuchos::RCP<Epetra_Operator> SolverFactory::CreateAlgebraicPrecond(Epetra_CrsMatrix& A,
                                                                        Teuchos::ParameterList& plist, int verbose)
    {
        std::string PrecType = plist.get("Method","Ifpack");
        DEBUG("Enter SolverFactory::CreateAlgebraicPrecond ("+PrecType+")");

        Teuchos::RCP<Epetra_Operator> prec;

        bool is_ifpack = (PrecType=="Ifpack");

        if (is_ifpack)
        {
            int OverlapLevel = plist.get("Ifpack Overlap Level",0);
            std::string SubType = plist.get("Ifpack Method","ILUT");
            if (SubType.find("MRILU") == 0)
            {
                int out = plist.sublist("MRILU").get("Output Level",0);
                if (verbose>5) out = max(out,verbose);
                if (verbose<5) out = min(out,verbose);
                plist.sublist("MRILU").set("Output Level",out);
            }

            Teuchos::RCP<Ifpack_Preconditioner> Prec;
            if (SubType == "MRILU")
                Prec = Teuchos::rcp(new Ifpack_AdditiveSchwarz<Ifpack_MRILU>(
                                        &A, OverlapLevel));
            else if (SubType == "MRILU stand-alone")
                Prec = Teuchos::rcp(new Ifpack_MRILU(&A));
            else
            {
                Ifpack PreconditionerFactory;
                Prec = Teuchos::rcp(PreconditionerFactory.Create(
                                        SubType, &A, OverlapLevel));
            }

            CHECK_ZERO(Prec->SetParameters(plist));
            CHECK_ZERO(Prec->Initialize());
            prec = Prec;
        }
        else if (PrecType=="ML")
        {
#ifndef NO_ML

            Teuchos::ParameterList& mllist = plist.sublist("ML");
            if (A.Comm().MyPID()==0)
            {
                int out = mllist.get("output",0);
                if (verbose>5) out = max(out,verbose);
                if (verbose<5) out = min(out,verbose);
                mllist.set("output",out);
                if (out>=10)
                {
                    plist.sublist("smoother: ifpack list").sublist("MRILU").set("Output Level",10);
                }

                if (out>0)
                {
                    std::cout << "Constructing Multi-Level Preconditioner for ";
                    std::cout << A.Label()<<"\n";
                }
            }
            if (mllist.get("smoother: type","Aztec")=="Aztec")
            {

                Teuchos::RCP<std::vector<int> > az_options =
                    Teuchos::rcp(new std::vector<int>(AZ_OPTIONS_SIZE));
                Teuchos::RCP<std::vector<double> > az_params =
                    Teuchos::rcp(new std::vector<double>(AZ_PARAMS_SIZE));

                AZ_defaults(&(*az_options)[0], &(*az_params)[0]);

                Teuchos::ParameterList& azlist = mllist.sublist("smoother: aztec list");

                // some reasonable default options for Krylov smoothers:
                (*az_options)[AZ_solver]=AZ_GMRESR;
                (*az_options)[AZ_scaling]=AZ_none;
                (*az_options)[AZ_precond]=AZ_dom_decomp;
                (*az_options)[AZ_subdomain_solve]=AZ_ilut;
                (*az_options)[AZ_max_iter]=3;
                (*az_options)[AZ_output]=0;
                (*az_options)[AZ_overlap]=0;
                (*az_options)[AZ_print_freq]=0;

                (*az_params)[AZ_tol]=0.0;
                (*az_params)[AZ_drop]=1.0e-12;
                (*az_params)[AZ_ilut_fill]=2.0;

                ExtractAztecOptions( azlist, &(*az_options)[0], &(*az_params)[0] );

                mllist.set("smoother: Aztec options",az_options);
                mllist.set("smoother: Aztec params",az_params);
            }

            if (mllist.isSublist("smoother: aztec list"))
            {
                mllist.remove("smoother: aztec list");
            }

            std::string smoo = mllist.get("smoother: type","Aztec");
            if (smoo=="IFPACK")
            {
                //Teuchos::ParameterList& ifp_list=mllist.sublist("smoother: ifpack list");
                std::string ifp_type=mllist.get("smoother: ifpack type","Amesos");

                // I experimented with Line relaxation in trilinos_thcm, but
                // it the DD-MRILU approach was far superior.
                if (ifp_type == "block relaxation")
                {
                    ERROR("line relaxation no longer supported",__FILE__,__LINE__);
                }
            }

            bool visualize = mllist.get("viz: enable",false);
            int repartition = mllist.get("repartition: enable",0);
            // I had something here in trilinos_thcm, but I don't know if it ever worked.
            if (visualize||repartition)
            {
                ERROR("visualization and repartitioning not supported",__FILE__,__LINE__);
            }
            prec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(A,mllist,false));
#else
            HYMLS::Tools::Error("ML is not available, choose another preconditioner!",__FILE__,__LINE__);
#endif
        }
        else if (PrecType=="ParaSails")
        {
#ifdef HAVE_PARASAILS
            prec = Teuchos::rcp(new ParaSailsPrecond(A,plist));
#else
            ERROR("ParaSails is not available, choose another preconditioner!",__FILE__,__LINE__);
#endif
        }
        else if (PrecType=="None")
        {
            prec=Teuchos::rcp(new IdentityOperator(A.RangeMap(),A.DomainMap(),A.Comm()));
        }
        else
        {
            ERROR("Bad Preconditioner Type: "+PrecType,__FILE__,__LINE__);
        }
        bool eigen_analysis = plist.get("Analyze Preconditioned Spectrum",false);
        // for the eigen-analysis after computing the precond,
        // we need access to the operator in a straight-forward way
        if (eigen_analysis)
        {
            Teuchos::RCP<const Epetra_Operator> eig_op = Teuchos::rcp(&A,false);
            plist.set("Eigen-Analysis: Operator", eig_op);
        }
        DEBUG("Leave SolverFactory::CreateAlgebraicPrecond ("+PrecType+")");
        return prec;
    }

// create an algebraic preconditinoer (i.e. for a submatrix)
    void SolverFactory::ComputeAlgebraicPrecond(Teuchos::RCP<Epetra_Operator> P, Teuchos::ParameterList& plist)
    {
        std::string PrecType = plist.get("Method","None");
        DEBUG("Enter SolverFactory::ComputeAlgebraicPrecond ("+PrecType+")");

        bool is_ifpack = (PrecType=="Ifpack");
        if (is_ifpack)
        {
            Teuchos::RCP<Ifpack_Preconditioner> Prec =
                Teuchos::rcp_dynamic_cast<Ifpack_Preconditioner>(P);

#ifdef DEBUGGING
            std::string SubType   = plist.get("Ifpack Method","None");
            int overlapLevel = plist.get("Ifpack Overlap Level", 0);
            DEBVAR(SubType);
            if (SubType=="ILU")
            {
                Teuchos::RCP<Ifpack_AdditiveSchwarz<Ifpack_ILU> > ilu =
                    Teuchos::rcp_dynamic_cast<Ifpack_AdditiveSchwarz<Ifpack_ILU> >(Prec, overlapLevel);
                const Epetra_RowMatrix& Aloc = dynamic_cast<const Epetra_RowMatrix&>(ilu->Inverse()->Matrix());
                //MatrixUtils::PrintRowMatrix(Aloc,*debug);
            }
            else if (SubType=="ILUT")
            {
                Teuchos::RCP<Ifpack_AdditiveSchwarz<Ifpack_ILUT> > ilu =
                    Teuchos::rcp_dynamic_cast<Ifpack_AdditiveSchwarz<Ifpack_ILUT> >(Prec, overlapLevel);
                const Epetra_RowMatrix& Aloc = dynamic_cast<const Epetra_RowMatrix&>(ilu->Inverse()->Matrix());
                //MatrixUtils::PrintRowMatrix(Aloc,*debug);
            }
            else if (SubType=="MRILU")
            {
                /*
                  Teuchos::RCP<Ifpack_AdditiveSchwarz<MRILU_Prec> > ilu =
                  Teuchos::rcp_dynamic_cast<Ifpack_AdditiveSchwarz<MRILU_Prec> >(Prec);
                */
            }
#endif
            // perform factorization...
            CHECK_ZERO(Prec->Compute());

#ifdef DEBUGGING
            if (SubType=="ILUT stand-alone")
            {
                Teuchos::RCP<Ifpack_ILUT> ilut = Teuchos::rcp_dynamic_cast<Ifpack_ILUT>(Prec);
                //DEBVAR(ilut->L());
                //DEBVAR(ilut->U());
            }
            else if (SubType=="ILU stand-alone")
            {
                Teuchos::RCP<Ifpack_ILU> ilu = Teuchos::rcp_dynamic_cast<Ifpack_ILU>(Prec);
                //DEBVAR(ilu->L());
                //DEBVAR(ilu->U());
            }
            else if (SubType=="ILU")
            {
                Teuchos::RCP<Ifpack_AdditiveSchwarz<Ifpack_ILU> > ilu =
                    Teuchos::rcp_dynamic_cast<Ifpack_AdditiveSchwarz<Ifpack_ILU> >(Prec);
                //DEBVAR(ilu->Inverse()->L());
                //DEBVAR(ilu->Inverse()->U());
            }
            else if (SubType=="ILUT")
            {
                Teuchos::RCP<Ifpack_AdditiveSchwarz<Ifpack_ILUT> > ilu =
                    Teuchos::rcp_dynamic_cast<Ifpack_AdditiveSchwarz<Ifpack_ILUT> >(Prec);
                //DEBVAR(ilu->Inverse()->L());
                //DEBVAR(ilu->Inverse()->U());
            }
#endif
        }
#ifndef NO_ML //[
        else if (PrecType=="ML")
        {

            Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> Prec =
                Teuchos::rcp_dynamic_cast<ML_Epetra::MultiLevelPreconditioner>(P);

            Teuchos::ParameterList& mllist = plist.sublist("ML");
            if (Prec->Comm().MyPID() == 0)
            {
                if (mllist.get("output", 0) > 0)
                {
                    std::cout << "Computing Multi-Level Preconditioner for ";
                    std::cout << Prec->RowMatrix().Label()<<"\n";
                }
            }

            Prec->ComputePreconditioner();

            // this has no effect unless the user sets "viz: enable" to true:
            if (mllist.get("viz: enable",false)==true)
            {
                Prec->VisualizeAggregates();
            }

//#ifdef TESTING //[
            bool analyze_cycle = plist.get("ML: Analyze Cycle",false);
            bool test_smoothers = plist.get("ML: Test Smoothers",false);
            if (test_smoothers || analyze_cycle)
            {
                // ML dumps all of its output to std::cout
                std::ostream& os = std::cout;
                //int NumSmoo = mllist.get("smoother: sweeps",1);
                int NumCycles = mllist.get("cycle applications",1);
                int NumPre=0, NumPost=0;
                std::string PreOrPost = mllist.get("smoother: pre or post","both");
                if ((PreOrPost=="pre")||(PreOrPost=="both"))
                {
                    NumPre=NumCycles;
                }
                if ((PreOrPost=="post")||(PreOrPost=="both"))
                {
                    NumPost=NumCycles;
                }


                if (analyze_cycle)
                {
                    if (Prec->Comm().MyPID()==0)
                    {
                        os << "Print hierarchy and analyze effect of MG cycle on a random vector\n";
                        os << "Matrix in question: " << Prec->RowMatrix().Label();
                    }
                    Prec->AnalyzeHierarchy(true,NumPre,NumPost,NumCycles);
                }
                if (test_smoothers)
                {
                    if (Prec->Comm().MyPID()==0)
                    {
                        os << "Test various Smoothers for ";
                        os << Prec->RowMatrix().Label()<<"\n";
                    }
                    Prec->TestSmoothers(plist.sublist("ML: Smoother Test"));
                }
            }
            bool dump_matrices = plist.get("ML: Dump Matrices",false);
            if (dump_matrices) DumpMLHierarchy(P);

//#endif //]
        }//ML
#endif //]
#ifdef HAVE_PARASAILS
        else if (PrecType=="ParaSails")
        {
            Teuchos::rcp_dynamic_cast<ParaSailsPrecond>(P)->Compute();
        }
#endif
        else if (PrecType=="None")
        {
            // ... //
        }
        else
        {
            ERROR("Bad preconditiner type: "+PrecType,__FILE__,__LINE__);
        }
        bool eigen_analysis = plist.get("Analyze Preconditioned Spectrum",false);
        if (eigen_analysis) AnalyzeSpectrum(plist,P);

        DEBUG("Leave SolverFactory::ComputeAlgebraicPrecond ("+PrecType+")");
    }

// note: we can currently only return the 'Teuchos::RCP<AztecOO>' type. Once Belos is
// available this should be redefined, but that means that Aztec will no longer
// be supported by our class.
    Teuchos::RCP<AztecOO> SolverFactory::CreateKrylovSolver(Teuchos::ParameterList& plist,int verbose)
    {
        Teuchos::RCP<AztecOO> Solver = Teuchos::null;
        std::string SolverType = plist.get("Method","AztecOO");
        if (SolverType == "AztecOO")
        {
            Solver = Teuchos::rcp(new AztecOO() );
            Solver->SetOutputStream(*outFile);
            Solver->SetErrorStream(*outFile);
            if (verbose >  5) plist.set("Output",1);
            if (verbose == 0) plist.set("Output",0);
            Solver->SetParameters(plist);
        }
        else if (SolverType != "None")
        {
            ERROR("Invalid Solver Method: "+SolverType,__FILE__,__LINE__);
        }
        return Solver;
    }



///////////////////////////////////////////////////////////////////////////////////////
// convert Teuchos::ParameterList entries to aztec options without an actual AztecOO              //
///////////////////////////////////////////////////////////////////////////////////////

    void SolverFactory::ExtractAztecOptions(Teuchos::ParameterList& azlist,int* options, double* params)
    {

// note: this is basically copied from class AztecOO in Trilinos 7.0.4

        Teuchos::ParameterList::ConstIterator
            pl_iter = azlist.begin(), pl_end  = azlist.end();

        for(; pl_iter != pl_end; ++pl_iter)
        {
            //create an upper-case copy of the entry's name and prepend AZ_ if necessary
//    std::string name = AztecOO_uppercase((*pl_iter).first);
            std::string name = Teuchos::StrUtils::allCaps((*pl_iter).first);

            if (!(name[0] == 'A' && name[1] == 'Z'))
            {
                std::string az_("AZ_");
                name=az_+name;
            }

            const Teuchos::ParameterEntry& entry = (*pl_iter).second;
            Teuchos::map<std::string,int>& azoo_key_map = AztecOO_key_map();
            Teuchos::map<std::string,int>::iterator result = azoo_key_map.find(name);
            bool entry_used = false;

            if (result != azoo_key_map.end())
            {
                int offset = (*result).second;
                if (offset < 0) return;

                int dummy_int;
                double dummy_double;
                std::string dummy_string;

                if (entry.isType<int>() || entry.isType<unsigned>())
                {
                    if (offset < AZ_FIRST_USER_OPTION)
                    {
                        int ival = entry.getValue(&dummy_int);
                        options[offset]=ival;
                        entry_used = true;
                    }
                }
                else if (entry.isType<std::string>())
                {
                    if (offset < AZ_FIRST_USER_OPTION)
                    {
                        std::string sname = Teuchos::StrUtils::allCaps(
                            entry.getValue(&dummy_string));
                        if (!(sname[0] == 'A' && sname[1] == 'Z'))
                        {
                            std::string az_("AZ_");
                            sname=az_+sname;
                        }
                        Teuchos::map<std::string,int>& val_map = AztecOO_value_map();
                        Teuchos::map<std::string,int>::iterator result = val_map.find(sname);
                        if (result != val_map.end())
                        {
                            options[offset] = (*result).second;
                            entry_used = true;
                        }
                    }
                }
                else if (entry.isType<double>())
                {
                    if (offset < AZ_FIRST_USER_PARAM)
                    {
                        double entry_value = entry.getValue(&dummy_double);
                        params[offset] = entry_value;
                        entry_used = true;
                    }
                }
            }
            if (!entry_used)
            {
                ERROR("your Aztec option '"+name+"' was not used!",__FILE__,__LINE__);
            }
        }
    }




    void SolverFactory::DumpMLHierarchy(Teuchos::RCP<Epetra_Operator> P)
    {
        Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> Prec =
            Teuchos::rcp_dynamic_cast<ML_Epetra::MultiLevelPreconditioner>(P);
        const ML* ml = Prec->GetML();
        int num_levels = ml->ML_num_actual_levels;
        Epetra_CrsMatrix* Ai;
        for (int i=0;i<num_levels;i++)
        {
            std::stringstream fname;
            fname<<"ML_Matrix_Level_"<<i<<".txt";
            CHECK_ZERO(ML_Operator2EpetraCrsMatrix(&(ml->Amat[i]), Ai));

            delete Ai;
        }
    }


// analyze spectrum of an operator inv(P)A-I, if P is a good preconditioner,
// this operator should have no dominant eigenmodes
    void SolverFactory::AnalyzeSpectrum(Teuchos::ParameterList& plist, Teuchos::RCP<const Epetra_Operator> Prec)
    {
#ifdef HAVE_ANASAZI
        Teuchos::RCP<const Epetra_Operator> A = null;
        A=plist.get("Eigen-Analysis: Operator", A);

        if (A->Comm().MyPID()==0)
        {
            std::cout << "##############################################################################\n";
            std::cout << "# Eigen-Analysis of preconditioner for "<<A->Label()<<std::endl;
            std::cout << "##############################################################################\n";
        }

        // the mass-matrix is not required. If it is non-null, the initial guess
        // is pre-multiplied by M, i.e. to get 0's into the conti-equation etc.
        Teuchos::RCP<Epetra_CrsMatrix> M = null;
        M=plist.get("Eigen-Analysis: Mass-Matrix", M);

        // this should have been set in 'CreateAlgebraicPrecond'
        if (A==null) ERROR("SolverFactory::AnalyzeSpectrum: matrix pointer not set",__FILE__,__LINE__);

        typedef double ScalarType;
        typedef Teuchos::ScalarTraits<ScalarType>          SCT;
        typedef SCT::magnitudeType               MagnitudeType;
        typedef Epetra_MultiVector                          MV;
        typedef Epetra_Operator                       OP;
        typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
        typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

        int ierr;

        // create operator P\A-I
        Teuchos::RCP<Epetra_Operator> MyOp;
        MyOp = Teuchos::rcp(new PAmI_Operator(A,Prec));

        // ************************************
        // Start the block Arnoldi iteration
        // ***********************************
        //
        //  Variables used for the Block Krylov Schur Method
        //
        bool boolret;
        int MyPID = Prec->Comm().MyPID();

        bool verbose = true;
        bool debug = false;
        std::string which("LM");

        int nev = 20;
        int blockSize = 1;
        int numBlocks = 50;
        int maxRestarts = 100;
        //int stepSize = 5;
        double tol = 1e-6;

        // Create a sort manager to pass into the block Krylov-Schur solver manager
        //  Make sure the reference-counted pointer is of type Anasazi::SortManager<>
        //  The block Krylov-Schur solver manager uses Anasazi::BasicSort<> by default,
        //  so you can also pass in the parameter "Which", instead of a sort manager.
        Teuchos::RCP<Anasazi::SortManager<ScalarType> > MySort =
            Teuchos::rcp( new Anasazi::BasicSort<ScalarType>( which ) );

        // Set verbosity level
        int verbosity = Anasazi::Errors + Anasazi::Warnings;
        if (verbose) {
            verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
        }
        if (debug) {
            verbosity += Anasazi::Debug;
        }
        //
        // Create parameter list to pass into solver manager
        //
        Teuchos::ParameterList MyPL;
        MyPL.set( "Verbosity", verbosity );
        MyPL.set( "Sort Manager", MySort );
        //MyPL.set( "Which", which );
        MyPL.set( "Block Size", blockSize );
        MyPL.set( "Num Blocks", numBlocks );
        MyPL.set( "Maximum Restarts", maxRestarts );
        //MyPL.set( "Step Size", stepSize );
        MyPL.set( "Convergence Tolerance", tol );

        // Create an Epetra_MultiVector for an initial vector to start the solver.
        // Note:  This needs to have the same number of columns as the blocksize.
        Teuchos::RCP<Epetra_MultiVector> ivec =
            Teuchos::rcp( new Epetra_MultiVector(A->OperatorRangeMap(), blockSize) );
        ivec->Random();

        if (M!=null)
        {
            Epetra_MultiVector tmp = *ivec;
            CHECK_ZERO(M->Multiply(false,tmp,*ivec));
        }

        // Create the eigenproblem.
        Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
            Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(MyOp, ivec) );

        // Inform the eigenproblem that the operator A is symmetric
        MyProblem->setHermitian(false);

        // Set the number of eigenvalues requested
        MyProblem->setNEV( nev );

        // Inform the eigenproblem that you are finishing passing it information
        boolret = MyProblem->setProblem();
        if (boolret != true) {
            if (verbose && MyPID == 0) {
                cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
            }
            ERROR("AnalyzeSpectrum: Could not create Anasazi Problem!",__FILE__,__LINE__);
        }

        // Initialize the Block Arnoldi solver
        Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

        // Solve the problem to the specified tolerances or length
        Anasazi::ReturnType returnCode;
        returnCode = MySolverMgr.solve();
        if (returnCode != Anasazi::Converged && MyPID==0 && verbose) {
            cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
        }

        // Get the Ritz values from the eigensolver
        std::vector<Anasazi::Value<double> > ritzValues = MySolverMgr.getRitzValues();


        // Output computed eigenvalues and their direct residuals
        if (verbose && MyPID==0) {
            int numritz = (int)ritzValues.size();
            cout.setf(std::ios_base::right, std::ios_base::adjustfield);
            cout<<endl<< "Computed Ritz Values"<< endl;
            if (MyProblem->isHermitian()) {
                cout<< std::setw(16) << "Real Part"
                    << endl;
                cout<<"-----------------------------------------------------------"<<endl;
                for (int i=0; i<numritz; i++) {
                    cout<< std::setw(16) << ritzValues[i].realpart
                        << endl;
                }
                cout<<"-----------------------------------------------------------"<<endl;
            }
            else {
                cout<< std::setw(16) << "Real Part"
                    << std::setw(16) << "Imag Part"
                    << endl;
                cout<<"-----------------------------------------------------------"<<endl;
                for (int i=0; i<numritz; i++) {
                    cout<< std::setw(16) << ritzValues[i].realpart
                        << std::setw(16) << ritzValues[i].imagpart
                        << endl;
                }
                cout<<"-----------------------------------------------------------"<<endl;
            }
        }
#endif
    }



}//namespace TRIOS
