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

    if (comm->NumProc() > 1)
        EXPECT_EQ(assemblyMap->UniqueGIDs(), false);
    else
        EXPECT_EQ(assemblyMap->UniqueGIDs(), true);

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
    
    EXPECT_EQ(standardMap->UniqueGIDs(), true);

    if (comm->NumProc() > 1)
        EXPECT_EQ(assemblyMap->UniqueGIDs(), false);
    else
        EXPECT_EQ(assemblyMap->UniqueGIDs(), true);

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

    if (comm->MyPID() == (comm->NumProc() - 1))
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

    if (comm->MyPID() == (comm->NumProc() - 1))
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
        atmos->setPar(0.1);
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

    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*standardMap) );
    Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(*standardMap) );
    
    x->PutScalar(1.0);

    double norm1 = Utils::norm(x);
    
    mat->Apply(*x, *b);

    atmos->solve(b);
    
    Teuchos::RCP<Epetra_Vector> x2 = atmos->getSolution('C');
    
    double norm2 = Utils::norm(x2);
    
    // There may be a difference due to domain overlap in the Ifpack
    // preconditioner that is used as a solver. So these norms should
    // be sort of similar.
    EXPECT_NEAR(norm1, norm2, 1e-1);
}

//------------------------------------------------------------------
TEST(Domain, Values)
{
    Teuchos::RCP<Epetra_Vector> x      = Teuchos::rcp(new Epetra_Vector(*standardMap) );
    Teuchos::RCP<Epetra_Vector> b      = Teuchos::rcp(new Epetra_Vector(*standardMap) );
    Teuchos::RCP<Epetra_Vector> localb = Teuchos::rcp(new Epetra_Vector(*assemblyMap) );

    x->Random();

    mat->Apply(*x, *b);

    domain->Solve2Assembly(*b, *localb);

    comm->Barrier();
    
    double P = 0.0;
    int last = n * m * l * dof + aux - 1;
    int lid  = 0;
    if ( standardMap->MyGID(last) )
    {
        lid = standardMap->LID(last);
        P = (*b)[lid];
    }
    
    comm->SumAll(&P, &P, 1);
    int numMyLocalElements = assemblyMap->NumMyElements();

    // The final element should be P
    EXPECT_EQ(P, (*localb)[numMyLocalElements-1]);    
}

//------------------------------------------------------------------
TEST(Domain, AtmosRHS)
{
    // Use idealized values in atmosphere model
    atmos->getState('V')->PutScalar(0.0);

    EXPECT_EQ(atmos->getState('V')->GlobalLength(), n * m * l * dof + aux);
    
    atmos->getSolution('V')->PutScalar(0.0);

    double defaultP = 123.456;
    
    atmos->idealized(defaultP);
    atmos->computeRHS();
    Teuchos::RCP<Epetra_Vector> rhs = atmos->getRHS('C');

    EXPECT_EQ( rhs->Map().SameAs(*standardMap), 1);
    
    double Prhs = 0.0;
    
    int last = n * m * l * dof + aux - 1;
    int lid = -1;
    
    if ( comm->MyPID() == (comm->NumProc() - 1) )
    {
        lid = standardMap->LID(last);
        Prhs = (*rhs)[lid];
    }
    
    comm->SumAll(&Prhs, &Prhs, 1);
    EXPECT_EQ(Prhs, -defaultP);
}

//------------------------------------------------------------------
TEST(Domain, numericalJacobian)
{
    // only do this test for small problems in serial
    int nmax = 2e3;
    
    if ( (comm->NumProc() == 1) &&
         (atmos->getState('V')->GlobalLength() < nmax) )
    {
        bool failed = false;    
        try
        {
            INFO("compute njC");
            
            NumericalJacobian<std::shared_ptr<AtmospherePar>,
                              Teuchos::RCP<Epetra_Vector> > njC;

            njC.setTolerance(1e-12);
            njC.seth(1e2);
            njC.compute(atmos, atmos->getState('V'));

            std::string fnameJnC("JnC");
            
            INFO(" Printing Numerical Jacobian " << fnameJnC);
            
            njC.print(fnameJnC);

            double njCsum = njC.sumValues();

            // create sum of values in mat
            Epetra_Vector ones(*standardMap);
            Epetra_Vector resl(*standardMap);
            ones.PutScalar(1.0);
            
            mat->Apply(ones, resl);
            
            int numMyElements = standardMap->NumMyElements();
            double sum = 0;
            for (int i = 0; i != numMyElements; ++i)
                sum += resl[i];
            
            comm->SumAll(&sum, &sum, 1);
            INFO("mat elements sum: " << sum);
            INFO("njc elements sum: " << njCsum);
            
            EXPECT_NEAR(sum, njCsum, 1e-4);

        }
        catch (...)
        {
            failed = true;
            throw;
        }
        EXPECT_EQ(failed, false);
    }
    
    if (comm->NumProc() != 1)
    {
        INFO("****Numerical Jacobian test cannot run in parallel****");
    }
    
    if (atmos->getState('V')->GlobalLength() > nmax)
    {
        INFO("****Numerical Jacobian test cannot run for this problem size****");
    }
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
    // TESTING
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;
    
    MPI_Finalize();
    return out;
}
