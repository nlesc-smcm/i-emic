#include "TestDefinitions.H"
#include "NumericalJacobian.H"

#include <limits>

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    RCP<TRIOS::Domain> domain;
    int n, m, l, dof, aux, periodic;
    double xmin,xmax,ymin,ymax;
}

//------------------------------------------------------------------
TEST(Domain, Aux0)
{
    bool failed = false;

    n = 6; m = 6; l = 1; dof = 2;

    xmin = 286 * PI_ / 180;
    xmax = 350 * PI_ / 180;
    ymin =  10 * PI_ / 180;
    ymax =  74 * PI_ / 180;

    periodic = 0;

    aux = 0;

    int dim = n * m * l * dof + aux;

    try
    {
        domain = Teuchos::rcp(new TRIOS::Domain(n, m, l, dof,
                                                xmin, xmax, ymin, ymax,
                                                periodic, 1.0, comm, aux));
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    Teuchos::RCP<Epetra_Map> colmap = domain->GetColMap();

    EXPECT_EQ(colmap->NumGlobalElements(), dim);
    
    std::cout << colmap->NumMyElements() << " "
              << colmap->NumGlobalElements() << " "
              << dim << std::endl;

    domain->Decomp2D();
}

//------------------------------------------------------------------
TEST(Domain, Aux2)
{
    bool failed = false;

    n = 6; m = 6; l = 1; dof = 2;

    xmin = 286 * PI_ / 180;
    xmax = 350 * PI_ / 180;
    ymin =  10 * PI_ / 180;
    ymax =  74 * PI_ / 180;

    periodic = 0;

    aux = 2;

    int dim = n * m * l * dof + aux;

    try
    {
        domain = Teuchos::rcp(new TRIOS::Domain(n, m, l, dof,
                                                xmin, xmax, ymin, ymax,
                                                periodic, 1.0, comm, aux));
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    
    EXPECT_EQ(failed, false);

    //////////////////////////////////////////////

    Teuchos::RCP<Epetra_Map> colmap = domain->GetColMap();
        
    EXPECT_EQ(colmap->NumGlobalElements(), dim);
    
    std::cout << colmap->NumMyElements() << " "
              << colmap->NumGlobalElements() << " "
              << dim << std::endl;
    
    //////////////////////////////////////////////

    domain->Decomp2D();

    Teuchos::RCP<Epetra_Map> standardMap = domain->GetStandardMap();
    Teuchos::RCP<Epetra_Map> assemblyMap = domain->GetAssemblyMap();

    Epetra_Vector vec(*standardMap);
    Epetra_Vector localvec(*assemblyMap);
    
    int numMyStandardElements = standardMap->NumMyElements();
    int numMyAssemblyElements = assemblyMap->NumMyElements();
    
    for (int i = 0; i != numMyStandardElements; ++i)
        vec[i] = i;

    for (int i = 0; i != numMyAssemblyElements; ++i)
        localvec[i] = i;

    EXPECT_NE( numMyAssemblyElements, numMyStandardElements );
    EXPECT_EQ( vec.Map().GID(numMyStandardElements-1), n*m*l*dof + aux - 1 );

    // next: test importing and exporting between assembly and standard mapped vectors
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    initializeEnvironment(argc, argv);
    if (outFile == Teuchos::null)
        throw std::runtime_error("ERROR: Specify output streams");

    ::testing::InitGoogleTest(&argc, argv);

    // -------------------------------------------------------
    // TESTINGv
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;
    
    MPI_Finalize();
    return out;
}
