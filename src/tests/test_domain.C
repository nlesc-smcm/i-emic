#include "TestDefinitions.H"
#include "NumericalJacobian.H"

#include <limits>

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    RCP<TRIOS::Domain> domain;
    RCP<Epetra_Map> standardMap;
    RCP<Epetra_Map> assemblyMap;
    RCP<Epetra_Map> stdSurfMap;
    RCP<Epetra_Map> asmSurfMap;
    RCP<Epetra_Import> as2std;
    RCP<Epetra_Import> as2std_surf;
    RCP<Epetra_Vector> vec;
    RCP<Epetra_Vector> localvec;
    
    
    int n, m, l, dof, aux, periodic;
    double xmin,xmax,ymin,ymax;
}

//------------------------------------------------------------------
TEST(Domain, SimpleInit)
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

    domain->Decomp2D();
}

//------------------------------------------------------------------
TEST(Domain, AuxInit)
{
    // Now we test some auxiliary unknowns
    aux = 2;
    
    bool failed = false;

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


    int dim = n * m * l * dof + aux;
        
    EXPECT_EQ(colmap->NumGlobalElements(), dim);
    
    ////////////////////////////////////////////
    // Create maps
    ////////////////////////////////////////////

    domain->Decomp2D();

    standardMap = domain->GetStandardMap();
    assemblyMap = domain->GetAssemblyMap();

    vec      = Teuchos::rcp(new Epetra_Vector(*standardMap));
    localvec = Teuchos::rcp(new Epetra_Vector(*assemblyMap));
    
    int numMyStandardElements = standardMap->NumMyElements();
    int numMyAssemblyElements = assemblyMap->NumMyElements();
    
    for (int i = 0; i != numMyStandardElements; ++i)
        (*vec)[i] = 100 + (*vec).Map().GID(i);

    int last = FIND_ROW2(dof, n, m, l, n-1, m-1, l-1, dof);
    EXPECT_EQ( vec->Map().LID(last + aux) , numMyStandardElements -1 );

    for (int i = 1; i <= aux; ++i)
        (*vec)[numMyStandardElements - aux - 1  + i] = 10000 + i - 1;
    
    
    for (int i = 0; i != numMyAssemblyElements; ++i)
        (*localvec)[i] = 1000 + localvec->Map().GID(i);

    
    EXPECT_EQ( vec->Map().GID(numMyStandardElements-1), n*m*l*dof + aux - 1 );

    ////////////////////////////////////////////
    // Create surface maps
    ////////////////////////////////////////////
    
    stdSurfMap = domain->CreateStandardMap(1, true);
    asmSurfMap = domain->CreateAssemblyMap(1, true);

    int numMyStdSurfElements = stdSurfMap->NumMyElements();
    int numMyAsmSurfElements = asmSurfMap->NumMyElements();

    if (comm->NumProc() > 1)
    {
        EXPECT_NE( numMyAssemblyElements, numMyStandardElements );
        EXPECT_NE( numMyAsmSurfElements,  numMyStdSurfElements  );
    }
    else
    {
        EXPECT_EQ( numMyAssemblyElements, numMyStandardElements ); 
        EXPECT_EQ( numMyAsmSurfElements,  numMyStdSurfElements  );
    }

    EXPECT_EQ( stdSurfMap->NumGlobalElements(), n*m );
}


//------------------------------------------------------------------
TEST(Domain, Importers)
{
    bool failed = false;
    try
    {
        as2std =
            Teuchos::rcp(new Epetra_Import(*assemblyMap, *standardMap));
        as2std_surf =
            Teuchos::rcp(new Epetra_Import(*asmSurfMap,  *stdSurfMap));

        // Import non-overlapping vec into overlapping localvec
        CHECK_ZERO( localvec->Import(*vec, *as2std, Insert) );
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    // get the local index for the global final unknown
    int last        = FIND_ROW2(dof, n, m, l, n-1, m-1, l-1, dof) + aux;
    int locvec_last = localvec->Map().LID(last);
    int vec_last    = vec->Map().LID(last);
    
    EXPECT_EQ( (*localvec)[locvec_last], 10000 + aux - 1 );

    int numMyStandardElements = standardMap->NumMyElements();
    int numMyAssemblyElements = assemblyMap->NumMyElements();

    // Refill...
    for (int i = 0; i != numMyStandardElements; ++i)
        (*vec)[i] = 100 + (*vec).Map().GID(i);
    
    for (int i = 0; i != numMyAssemblyElements; ++i)
        (*localvec)[i] = 1000 + localvec->Map().GID(i);
    
    for (int i = 1; i <= aux; ++i)
        (*localvec)[locvec_last - aux + i] = 10000 + i - 1;
    
    try
    {
        // Export overlapping localvec into non-overlapping vec
        CHECK_ZERO( vec->Export( *localvec, *as2std, Zero ) );
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    EXPECT_EQ( (*vec)[vec_last], 10000 + aux - 1 );
}

//------------------------------------------------------------------
TEST(Domain, Gather)
{
    int last = FIND_ROW2(dof, n, m, l, n-1, m-1, l-1, dof) + aux;
    Teuchos::RCP<Epetra_MultiVector> gvec = Utils::Gather(*vec, 0);
    EXPECT_EQ( gvec->GlobalLength(), last + 1 );     
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
