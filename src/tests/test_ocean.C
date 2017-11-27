#include "TestDefinitions.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{

    RCP<Ocean>  ocean;
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
        RCP<Teuchos::ParameterList> oceanParams =
            rcp(new Teuchos::ParameterList);
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
            Teuchos::RCP<Epetra_CrsMatrix> mat = ocean->getJacobian();
            DUMPMATLAB("ocean_jac", *mat);
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
