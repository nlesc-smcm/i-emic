#include "TestDefinitions.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{

    RCP<Teuchos::ParameterList> oceanParams;
    RCP<Ocean> ocean;  
    RCP<Epetra_Comm>  comm; 
}

//------------------------------------------------------------------
class IEMIC : public testing::Environment
{
public:
    // constructor
    IEMIC()
        {}

    // destructor
    ~IEMIC()
        {}
};


//------------------------------------------------------------------
TEST(Ocean, Initialization)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        oceanParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
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

    v.PutScalar(1.0);
    ocean->applyMassMat(v, out);

    std::ofstream file;
    file.open("massmat");
    file << out;
    file.close();       

    int numMyElements = out.Map().NumMyElements();
    double rosb = ocean->getPar("Rossby-Number");
    
    for (int i = 0; i < numMyElements; i+=_NUN_)
    {
        if (std::abs(out[i])>0) // UU
            EXPECT_EQ(out[i], rosb);
        
        if (std::abs(out[i+1])>0) // VV
            EXPECT_EQ(out[i+1], rosb);
        
        EXPECT_EQ(out[i+2], 0.0); // WW
        EXPECT_EQ(out[i+3], 0.0); // PP
        EXPECT_EQ(out[i+4], 1.0); // TT
        if (std::abs(out[i+5])>0) // SS
            EXPECT_EQ(out[i+5], 1.0); 
    }

    // check integral condition
    int sres = oceanParams->get("Restoring Salinity Profile", 1);

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
            ocean->setPar(0.1);
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
    int nmax = 2e3;

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
        ocean->setPar(0.0);
        ocean->getState('V')->PutScalar(0.0);
        
        // Create continuation params
        RCP<Teuchos::ParameterList> continuationParams =
            rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params.xml",
                                    continuationParams.ptr());

        // Create contination
        Continuation<RCP<Ocean>, RCP<Teuchos::ParameterList> >
            continuation(ocean, continuationParams);

        // Run continuation
        continuation.run();

        Teuchos::RCP<Epetra_CrsMatrix> mat = ocean->getJacobian();
        DUMPMATLAB("ocean_jac", *mat);
        
        Teuchos::RCP<Epetra_Vector> diagB = ocean->getDiagB();
        EXPECT_NE(Utils::norm(diagB), 0.0);
        DUMP_VECTOR("ocean_B", *diagB);                        
        
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
    ::testing::AddGlobalTestEnvironment(new IEMIC);

    // -------------------------------------------------------
    // TESTING
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    // Get rid of possibly parallel objects for a clean ending.
    ocean = Teuchos::null;

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
