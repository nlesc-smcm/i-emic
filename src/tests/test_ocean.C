#include "TestDefinitions.H"

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "NumericalJacobian.H"
#include "THCMdefs.H"
#include "Ocean.H"
#include "Continuation.H"

#include "TRIOS_Domain.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Teuchos::ParameterList> oceanParams;
    Teuchos::RCP<Ocean> ocean;
    Teuchos::RCP<Epetra_Comm> comm;
}

//------------------------------------------------------------------
TEST(Ocean, Initialization)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        oceanParams = Utils::obtainParams("ocean_params.xml", "Ocean parameters");
        oceanParams->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Ocean, RHSNorm)
{
    ocean->computeRHS();
    double stateNorm = Utils::norm(ocean->getState('V'));
    double rhsNorm   = Utils::norm(ocean->getRHS('V'));
    std::cout << "stateNorm = " << stateNorm << std::endl;
    std::cout << "RHSNorm   = " << rhsNorm   << std::endl;
    EXPECT_LT(rhsNorm, 1e-6);
}

//------------------------------------------------------------------
// Check mass matrix contents
TEST(Ocean, MassMat)
{
    Epetra_Vector v   = *ocean->getState('C');
    Epetra_Vector out = *ocean->getState('C');

    Utils::MaskStruct mask = ocean->getLandMask();

    v.PutScalar(1.0);
    ocean->applyMassMat(v, out);

    std::ofstream file;
    file.open("massmat");
    file << out;
    file.close();

    int dim = mask.global_borderless->size();

    EXPECT_EQ(dim * _NUN_, out.GlobalLength());

    double rosb = ocean->getPar("Rossby-Number");

    int gid = -1, lid = -1;
    for (int i = 0; i < dim; ++i)
    {
        gid = i * _NUN_;
        lid = out.Map().LID(gid);
        if ((lid >= 0) && ((*mask.global_borderless)[i] == 0))
        {
            if (std::abs(out[lid])>0) // UU
            {
                EXPECT_EQ(out[lid], -rosb);
            }

            if (std::abs(out[lid+1])>0) // VV
            {
                EXPECT_EQ(out[lid+1], -rosb);
            }

            EXPECT_EQ(out[lid+2], 0.0);  // WW
            EXPECT_EQ(out[lid+3], 0.0);  // PP
            EXPECT_EQ(out[lid+4], -1.0); // TT
            if (std::abs(out[lid+5])>0)  // SS
            {
                EXPECT_EQ(out[lid+5], -1.0);
            }
        }
    }

    // check integral condition
    int sres = oceanParams->sublist("THCM").get("Restoring Salinity Profile", 1);

    if (sres == 0)
    {
        Teuchos::RCP<TRIOS::Domain> domain = ocean->getDomain();
        int N = domain->GlobalN();
        int M = domain->GlobalM();
        int L = domain->GlobalL();
        int rowIntCon = FIND_ROW2(_NUN_,N,M,L,N-1,M-1,L-1,SS);
        if (out.Map().MyGID(rowIntCon))
        {
            int lid = out.Map().LID(rowIntCon);
            EXPECT_EQ(out[lid], 0.0);
        }
    }
}

//------------------------------------------------------------------
TEST(Ocean, ComputeJacobian)
{
    if (comm->NumProc() > 1)
    {
        WARNING("We are not going to do this in parallel...",
                __FILE__, __LINE__);
    }
    else
    {
        bool failed = false;
        try
        {
            ocean->setPar("Combined Forcing", 0.1);
            ocean->getState('V')->PutScalar(1.234);
            ocean->computeJacobian();
        }
        catch (...)
        {
            failed = true;
        }
        EXPECT_EQ(failed, false);
    }
}

//------------------------------------------------------------------
TEST(Ocean, NumericalJacobian)
{
    // only do this test for small problems in serial
    int nmax = 2.0e3;

    if ( (comm->NumProc() == 1) &&
         (ocean->getState('V')->GlobalLength() < nmax) )
    {
        bool failed = false;
        try
        {
            NumericalJacobian<Teuchos::RCP<Ocean>,
                              Teuchos::RCP<Epetra_Vector> > njmat;

            njmat.setTolerance(1e-10);
            njmat.seth(1);
            njmat.compute(ocean, ocean->getState('V'));

            std::string fname("ocean_numjac");

            INFO(" Printing numerical Jacbian " << fname);

            njmat.print(fname);
        }
        catch (...)
        {
            failed = true;
        }
        EXPECT_EQ(failed, false);
    }
    else
    {
        std::cout << "Problem size: " << ocean->getState('V')->GlobalLength()
                  << ", number of procs: " << comm->NumProc() << std::endl;

        WARNING("We are not going to do this for this setup...",
                __FILE__, __LINE__);
    }
}

//------------------------------------------------------------------
TEST(Ocean, Continuation)
{
    bool failed = false;
    try
    {
        ocean->setPar("Combined Forcing", 0.0);
        ocean->getState('V')->PutScalar(0.0);

        // Create continuation params
        Teuchos::RCP<Teuchos::ParameterList> continuationParams =
            Teuchos::rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params.xml",
                                    continuationParams.ptr());

        // Create contination
        Continuation<Teuchos::RCP<Ocean>> continuation(ocean, continuationParams);

        Teuchos::RCP<Epetra_CrsMatrix> mat = ocean->getJacobian();
        DUMPMATLAB("ocean_jac", *mat);

        Teuchos::RCP<Epetra_Vector> diagB = ocean->getMassMat();
        EXPECT_NE(Utils::norm(diagB), 0.0);
        DUMP_VECTOR("ocean_B", *diagB);

        // Run continuation
        int status = continuation.run();
        EXPECT_EQ(status, 0);

        mat  = ocean->getJacobian();
        DUMPMATLAB("ocean_jac", *mat);

        diagB = ocean->getMassMat();
        EXPECT_NE(Utils::norm(diagB), 0.0);
        DUMP_VECTOR("ocean_B", *diagB);

        Teuchos::RCP<Epetra_Vector> intcond_coeff = ocean->getIntCondCoeff();
        DUMP_VECTOR("intcond_coeff", *intcond_coeff);

    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//-------------------------------------------------------------------
TEST(Ocean, Integrals)
{
    double salt_advection = 0.0;
    double salt_diffusion = 0.0;

    Teuchos::RCP<Epetra_Vector> un = ocean->getState('C');

    ocean->integralChecks(un, salt_advection, salt_diffusion);

    EXPECT_NEAR(salt_advection, 0.0, 1e-10);

    // This is not such a great test as it does not check the actual
    // discretisation.
    // EXPECT_NEAR(salt_diffusion, 0.0, 1e-10);

    Teuchos::RCP<Epetra_Vector> icCoef = ocean->getIntCondCoeff();
    Teuchos::RCP<Epetra_Vector> x      = ocean->getSolution('V');
    Teuchos::RCP<Epetra_CrsMatrix> mat = ocean->getJacobian();


    Teuchos::RCP<Epetra_Vector> e   = Teuchos::rcp(new Epetra_Vector(x->Map()));
    Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(x->Map()));
    Teuchos::RCP<Epetra_Vector> col = Teuchos::rcp(new Epetra_Vector(x->Map()));

    Teuchos::RCP<TRIOS::Domain> domain = ocean->getDomain();
    int N = domain->GlobalN();
    int M = domain->GlobalM();
    int L = domain->GlobalL();

    int rowS, lid;

    int rowIntCon = ocean->getRowIntCon();

    int lidIntCon = e->Map().LID(rowIntCon);
    if (lidIntCon >= 0)
        (*icCoef)[lidIntCon] = 0;

    TIMER_START("Test ocean: integral method 1");
    std::vector<double> integrals;
    double dot;
    for (int k = 0; k != L; ++k)
        for (int j = 0; j != M; ++j)
            for (int i = 0; i != N; ++i)
            {
                rowS = FIND_ROW2(_NUN_,N,M,L,i,j,k,SS);

                lid  = e->Map().LID(rowS);

                if (lid >= 0)
                    (*e)[lid] = 1;

                col->PutScalar(0.0);
                mat->Apply(*e, *col);
                dot = Utils::dot(icCoef, col);

                if (lid >= 0)
                {
                    (*e)[lid]   = 0;
                    (*tmp)[lid] = dot;
                }

                if (k < L-2) // here we test everything except the top rows
                    integrals.push_back(dot);
            }

    for (auto &el: integrals)
    {
        EXPECT_NEAR(el, 0.0, 1e-7);
    }
    TIMER_STOP("Test ocean: integral method 1");

    TIMER_START("Test ocean: integral method 2");
    Teuchos::RCP<Epetra_Vector> integrals2 = ocean->getColumnIntegral(mat);
    EXPECT_NEAR(Utils::norm(tmp), Utils::norm(integrals2), 1e-7);
    TIMER_STOP("Test ocean: integral method 2");

    int badRows = ocean->analyzeJacobian1();
    std::cout << " bad P rows: " << badRows << std::endl;

    badRows = ocean->analyzeJacobian1();
    std::cout << " bad S ints: " << badRows << std::endl;
}

//------------------------------------------------------------------
TEST(Ocean, CreateASecondOcean)
{
    // First destroy the old Ocean since there can only be one
    // instance of THCM.
    ocean = Teuchos::null;

    bool failed = false;
    try
    {
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    comm = initializeEnvironment(argc, argv);
    if (outFile == Teuchos::null)
        throw std::runtime_error("ERROR: Specify output streams");

    ::testing::InitGoogleTest(&argc, argv);

    // -------------------------------------------------------
    // TESTING
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    // Get rid of possibly parallel objects for a clean ending.
    ocean = Teuchos::null;

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    if (comm->MyPID() == 0)
        printProfile();

    MPI_Finalize();
    return out;
}
