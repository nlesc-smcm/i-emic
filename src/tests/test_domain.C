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

    RCP<Epetra_CrsMatrix> mat;

    std::shared_ptr<AtmospherePar> atmos;

    RCP<Teuchos::ParameterList> atmosphereParams;    
    
    int n, m, l, dof, aux, periodic;
    double xmin,xmax,ymin,ymax;
}

//------------------------------------------------------------------
TEST(Atmosphere, Initialization)
{
    bool failed = false;

    // Create atmosphere parameters
    atmosphereParams = rcp(new Teuchos::ParameterList);

    updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());

    try
    {
        atmos = std::make_shared<AtmospherePar>(comm, atmosphereParams);
    }
    catch (std::exception const &e)
    {
        INFO("TEST(Atmosphere, Initialization) exception: " << e.what());
        failed = true;
    }
    catch (int error)
    {
        INFO("TEST(Atmosphere, Initialization) exception: error code = " << error);
        failed = true;
    }
    catch (...)
    {
        INFO("TEST(Atmosphere, Initialization) some exception thrown...");
        failed = true;
        throw;
    }
    
    EXPECT_EQ(failed, false);
}


//------------------------------------------------------------------
TEST(Domain, SimpleInit)
{
    bool failed = false;

    dof = ATMOS_NUN_;

    xmin = 286 * PI_ / 180;
    xmax = 350 * PI_ / 180;
    ymin =  10 * PI_ / 180;
    ymax =  74 * PI_ / 180;

    xmin = atmosphereParams->get("Global Bound xmin", 286.0) * PI_ / 180.0;
    xmax = atmosphereParams->get("Global Bound xmax", 350.0) * PI_ / 180.0;
    ymin = atmosphereParams->get("Global Bound ymin", 10.0)  * PI_ / 180.0;
    ymax = atmosphereParams->get("Global Bound ymax", 74.0)  * PI_ / 180.0;

    n        = atmosphereParams->get("Global Grid-Size n", 16);
    m        = atmosphereParams->get("Global Grid-Size m", 16);
    l        = atmosphereParams->get("Global Grid-Size l", 1);
    periodic = atmosphereParams->get("Periodic", false);

    // without auxiliary unknowns
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


    ////////////////////////////////////////////
    // Create maps
    ////////////////////////////////////////////

    domain->Decomp2D();

    standardMap = domain->GetStandardMap();
    assemblyMap = domain->GetAssemblyMap();
    
    EXPECT_EQ(standardMap->UniqueGIDs(), true);
    EXPECT_EQ(assemblyMap->UniqueGIDs(), false);

}

//------------------------------------------------------------------
TEST(Domain, AuxInit)
{
    // Now we test some auxiliary unknowns
    aux = atmosphereParams->get("Auxiliary unknowns", 2);

    if (aux <= 0)
    {
        WARNING(" aux not set, test has no use", __FILE__, __LINE__);
        std::cout << " aux not set, test has no use" << std::endl;
        return;
    }
        
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
    
    std::cout << *standardMap << std::endl;
    std::cout << *assemblyMap << std::endl;
    getchar();

    EXPECT_EQ(standardMap->UniqueGIDs(), true);
    EXPECT_EQ(assemblyMap->UniqueGIDs(), false);

    vec      = Teuchos::rcp(new Epetra_Vector(*standardMap));
    localvec = Teuchos::rcp(new Epetra_Vector(*assemblyMap));
    
    int numMyStandardElements = standardMap->NumMyElements();
    int numMyAssemblyElements = assemblyMap->NumMyElements();
    
    for (int i = 0; i != numMyStandardElements; ++i)
        (*vec)[i] = 100 + (*vec).Map().GID(i);

    int last = FIND_ROW2(dof, n, m, l, n-1, m-1, l-1, dof);
    if (comm->MyPID() == (comm->NumProc() - 1))
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
    if (aux <= 0)
    {
        WARNING(" aux not set, test has no use", __FILE__, __LINE__);
        std::cout << " aux not set, test has no use" << std::endl;
        return;
    }


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

    EXPECT_EQ(vec->GlobalLength(), standardMap->NumGlobalElements());
    EXPECT_EQ(vec->GlobalLength(), domain->GetSolveMap()->NumGlobalElements());
    
}

//------------------------------------------------------------------
TEST(Domain, Gather)
{
    if (aux <= 0)
    {
        WARNING(" aux not set, test has no use", __FILE__, __LINE__);
        std::cout << " aux not set, test has no use" << std::endl;
        return;
    }

    int last = FIND_ROW2(dof, n, m, l, n-1, m-1, l-1, dof) + aux;

    Teuchos::RCP<Epetra_MultiVector> gvec = Utils::Gather(*vec, comm->NumProc() - 1);
    EXPECT_EQ( gvec->GlobalLength(), last + 1 );

    Teuchos::RCP<Epetra_Vector> intCondCoeff = atmos->getIntCoeff();

    Teuchos::RCP<Epetra_MultiVector> gint =
        Utils::Gather(*intCondCoeff, comm->NumProc() - 1);
    
    EXPECT_EQ( gint->GlobalLength(), last + 1 );    
}


//------------------------------------------------------------------
TEST(Domain, MatVec)
{

    bool failed = false;
    try
    {
        atmos->getState('V')->PutScalar(0.01);
        atmos->setPar(0.01);
        atmos->computeJacobian();
        mat = atmos->getJacobian();
    }
    catch (std::exception const &e)
    {
        INFO("TEST(Domain, MatVec) exception: " << e.what());
        failed = true;
    }
    catch (int error)
    {
        INFO("TEST(Domain, MatVec) exception: error code = " << error);
        failed = true;
    }
    catch (...)
    {
        INFO("TEST(Domain, MatVec) some exception thrown...");
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*standardMap));
    Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(*standardMap));
    
    x->PutScalar(1.0);
    
    mat->Apply(*x, *b);

    // std::cout << *mat << std::endl;
    
    atmos->solve(b);
    
    Teuchos::RCP<Epetra_Vector> x2 = atmos->getSolution('C');
    
    //std::cout << *x2 << std::endl;
    //std::cout << Utils::norm(x2) << std::endl;
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
