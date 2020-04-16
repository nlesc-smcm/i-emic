#include "TestDefinitions.H"

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "NumericalJacobian.H"
#include "Atmosphere.H"
#include "Ocean.H"

#include "Epetra_Import.h"

// Testing the serial and parallel atmosphere

//------------------------------------------------------------------
namespace
{
    std::shared_ptr<AtmosLocal>          atmosLoc;
    std::shared_ptr<Atmosphere>          atmosPar;
    Teuchos::RCP<Epetra_Comm>            comm;
    Teuchos::RCP<Teuchos::ParameterList> atmosphereParams;
    Utils::MaskStruct                    mask;
}

//------------------------------------------------------------------
TEST(Atmosphere, Initialization)
{
    bool failed = false;
    // Create atmosphere parameters
    atmosphereParams = Teuchos::rcp(new Teuchos::ParameterList);
    updateParametersFromXmlFile("atmosphere_params.xml", atmosphereParams.ptr());

    try
    {
        atmosLoc = std::make_shared<AtmosLocal>(atmosphereParams);
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    failed = false;
    try
    {
        atmosPar = std::make_shared<Atmosphere>(comm, atmosphereParams);
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

    atmosPar->setPar("Combined Forcing", 0.0);
    atmosLoc->setPar("Combined Forcing", 0.0);

    atmosPar->computeRHS();
    atmosLoc->computeRHS();

    // initial rhs should be equal
    double normPar = Utils::norm(atmosPar->getRHS('V'));
    double normSer = Utils::norm(*atmosLoc->getRHS('V'));
    EXPECT_NEAR(normPar, normSer, 1e-7);
    
    // initial rhs should be small
    EXPECT_LT(normPar, 1e-4);}

//------------------------------------------------------------------
TEST(Atmosphere, RHS)
{
    bool failed = false;
    try
    {
        // Check dimensions
        int dim1 = atmosLoc->dim();
        int dim2 = atmosPar->dim();
        EXPECT_EQ(dim1, dim2);

        // Create random state from atmosPar standardMap state
        Teuchos::RCP<Epetra_Vector> state_par = atmosPar->getState('C');

        state_par->SetSeed(rand() % 100);
        INFO("TEST(Atmosphere, RHS): seed = " << state_par->Seed());

        // Randomize and scale standardMap state
        state_par->Random();

        // Gather state_par into std::vector
        std::shared_ptr<std::vector<double> > state =
            getGatheredVector(state_par);

        // Put the full state in the serial atmosphere
        atmosLoc->setState(state);

        // Put the standardMap parallel state into the parallel atmosphere
        atmosPar->setState(state_par);

        // Compute RHS in both models
        atmosLoc->computeRHS();
        atmosPar->computeRHS();

        // Obtain RHS from both models
        Teuchos::RCP<Epetra_Vector> rhsPar = atmosPar->getRHS('V');
        std::shared_ptr<std::vector<double> > rhsSer = atmosLoc->getRHS('V');

        // Check norms parallel and serial rhs
        EXPECT_NEAR(norm(rhsPar), norm(rhsSer), 1e-7);
        INFO("TEST(Atmosphere, RHS): serial norm = " << norm(rhsSer));
        INFO("TEST(Atmosphere, RHS): parall norm = " << norm(rhsPar));
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, State)
{
    bool failed = false;
    try
    {
        // Obtain state from both models
        Teuchos::RCP<Epetra_Vector> statePar = atmosPar->getState('V');
        std::shared_ptr<std::vector<double> > stateSer = atmosLoc->getState('V');

        // Check norms parallel and serial rhs
        EXPECT_NEAR(norm(statePar), norm(stateSer), 1e-7);
    }
    catch(...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, SurfaceTemperature)
{
    bool failed = false;
    try
    {
        // Obtain map for surface valeus
        Teuchos::RCP<Epetra_BlockMap> surfmap = atmosPar->getStandardSurfaceMap();

        // Create random vector from surface map
        Teuchos::RCP<Epetra_Vector> surfvals = Teuchos::rcp(new Epetra_Vector(*surfmap));
        surfvals->SetSeed(rand() % 100);
        INFO("TEST(Atmosphere, SurfaceTemperature): seed = " << surfvals->Seed());

        // Randomize and scale
        surfvals->Random();
        surfvals->Scale(1e1);

        // Create gather map
        Teuchos::RCP<Epetra_BlockMap> gmap =
            Utils::AllGather(*surfmap);

        // Create import strategy
        Teuchos::RCP<Epetra_Import> gimp =
            Teuchos::rcp(new Epetra_Import(*gmap, *surfmap));

        // Create gather vector
        Teuchos::RCP<Epetra_Vector> gathered = Teuchos::rcp(new Epetra_Vector(*gmap));

        // Gather <surfvals> into <gathered>
        gathered->Import(*surfvals, *gimp, Insert);

        // Create vector for serial atmosphere
        int size = surfmap->NumGlobalElements();
        std::shared_ptr<std::vector<double> > fullsurfvals =
            std::make_shared<std::vector<double> >(size);

        // Extract gathered
        gathered->ExtractCopy(&(*fullsurfvals)[0], size);

        // Check norms randomized parallel and serial sst
        EXPECT_NEAR(norm(surfvals), norm(fullsurfvals), 1e-6);
        INFO("TEST SST: serial SST norm = " << norm(surfvals));
        INFO("TEST SST: parall SST norm = " << norm(fullsurfvals));

        // Put the randomized ocean surface values in both models
        atmosLoc->setOceanTemperature(*fullsurfvals);
        atmosPar->setOceanTemperature(surfvals);

        // Compute RHS in both models
        atmosLoc->computeRHS();
        atmosPar->computeRHS();

        // Obtain RHS from both models
        Teuchos::RCP<Epetra_Vector> rhsPar = atmosPar->getRHS('V');
        std::shared_ptr<std::vector<double> > rhsSer = atmosLoc->getRHS('V');

        EXPECT_NEAR(norm(rhsPar), norm(rhsSer), 1e-6);
        INFO("TEST SST: serial RHS norm = " << norm(rhsSer));
        INFO("TEST SST: parall RHS norm = " << norm(rhsPar));

    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}
//------------------------------------------------------------------
TEST(Atmosphere, Pdist)
{
    atmosLoc->computeRHS();
    atmosLoc->computePrecipitation();
    atmosPar->computeRHS();

    // get serial and parallel integral coefficients
    std::shared_ptr<std::vector<double> > serPco = atmosLoc->getPIntCoeff();
    Teuchos::RCP<Epetra_Vector> parPco = atmosPar->getPIntCoeff();
    Teuchos::RCP<Epetra_Vector> ones =
        Teuchos::rcp(new Epetra_Vector(*parPco));
    ones->PutScalar(1.0);
    double totalArea = Utils::dot(parPco, ones);

    std::vector<double> locPdist = *atmosLoc->getPIntCoeff();
    atmosLoc->fillPdist(&locPdist[0]);
    double locPdistInt = Utils::dot(locPdist, *serPco);

    Teuchos::RCP<Epetra_Vector> parPdist = atmosPar->getPdist();
    double parPdistInt = Utils::dot(parPco, parPdist);
    std::cout << parPdistInt << " " << totalArea << std::endl;

    double parPdistNorm = Utils::norm(parPdist);
    double locPdistNorm = Utils::norm(locPdist);
    std::cout << parPdistNorm << std::endl;
    std::cout << locPdistNorm << std::endl;

    for (auto &e: locPdist)
    {
        if (std::abs(e) > 0.0)
        {
            e = e + 1 - locPdistInt / totalArea;
        }
    }
    
    double locPdistNorm2 = Utils::norm(locPdist);
    EXPECT_NEAR(parPdistNorm, locPdistNorm2, 1e-7);
}

//------------------------------------------------------------------
TEST(Atmosphere, EPfields)
{
    atmosLoc->computeRHS();
    atmosPar->computeRHS();

    // get serial integral coefficients and evaporation field
    std::shared_ptr<std::vector<double> > serPco = atmosLoc->getPIntCoeff();
    std::shared_ptr<std::vector<double> > serE   = atmosLoc->interfaceE();

    // get parallel integral coefficients and evaporation field
    Teuchos::RCP<Epetra_Vector> parPco = atmosPar->getPIntCoeff();
    Teuchos::RCP<Epetra_Vector> parE   = atmosPar->interfaceE();

    // compute area
    double serArea, parArea;
    serArea = Utils::sum(*serPco);
    parPco->Norm1(&parArea);
    EXPECT_NEAR(serArea, parArea,1e-7);

    std::cout << " area = " << std::setprecision(12) << serArea << " (serial model) "
              << parArea << " (parallel model) " << std::endl;
    
    // compute dot products / integrals
    double serInt, parInt;
    serInt = Utils::dot(*serPco, *serE) / serArea;
    parInt = Utils::dot(parPco, parE) / parArea;
    EXPECT_NEAR(serInt, parInt, 1e-7);
    EXPECT_NE(serInt, 0.0);
    EXPECT_NE(parInt, 0.0);

    atmosLoc->computePrecipitation();

    std::shared_ptr<std::vector<double> > serP = atmosLoc->interfaceP();
    Teuchos::RCP<Epetra_Vector> parP = atmosPar->interfaceP();

    double serNrm, parNrm;
    serNrm = Utils::norm(*serP);
    parNrm = Utils::norm(parP);
    std::cout << " serial P: " << serNrm << " parallel P: " << parNrm << std::endl;

    for (auto &e: *serP)
        std::cout << e << " \n";
    std::cout << std::endl;

    std::cout << *parP << std::endl;
    
    EXPECT_NEAR(serNrm, parNrm, 1e-10);
    EXPECT_GT(serNrm, 0.0);
}

//------------------------------------------------------------------
TEST(Atmosphere, MassMatrix)
{
    Epetra_Vector v   = *atmosPar->getState('C');
    Epetra_Vector out = *atmosPar->getState('C');

    v.PutScalar(1.0);
    atmosPar->applyMassMat(v, out);

    std::ofstream file;
    file.open("massmat");
    file << out;
    file.close();

    int numMyElements = out.Map().NumMyElements();

    double rhoa    = atmosphereParams->get("atmospheric density",1.25);
    double hdima   = atmosphereParams->get("atmospheric scale height",8400.);
    double cpa     = atmosphereParams->get("heat capacity",1000.);
    double udim    = atmosphereParams->get("horizontal velocity of the ocean", 0.1e+00);
    double r0dim   = atmosphereParams->get("radius of the earth", 6.37e+06);
    double ce      = atmosphereParams->get("Dalton number",1.3e-03);
    double ch      = atmosphereParams->get("exchange coefficient ch",0.94 * ce);
    double uw      = atmosphereParams->get("mean atmospheric surface wind speed",8.5);
//     double qdim    = atmosphereParams->get("humidity scale", 0.01);  // (kg/kg)

    double muoa    =  rhoa * ch * cpa * uw;
    double Ai      =  rhoa * hdima * cpa * udim / (r0dim * muoa);

    for (int i = 0; i < numMyElements; i+=ATMOS_NUN_)
    {
        EXPECT_EQ(out[i], Ai); // TT
        if (std::abs(out[i+1])>0) // QQ
        {
            EXPECT_EQ(out[i+1], 1.0);
        }
    }

    // check integral equations in test_coupled
}

//------------------------------------------------------------------
TEST(Atmosphere, Jacobian)
{
    bool failed = false;
    try
    {
        atmosLoc->computeJacobian();
        atmosPar->computeJacobian();

        // Obtain vector and randomize
        Teuchos::RCP<Epetra_Vector> x = atmosPar->getState('C');
        x->Random();
        x->Scale(1e1);

        // Gather randomized vector
        std::shared_ptr<std::vector<double> > gx = getGatheredVector(x);

        // Check the gather routine
        EXPECT_NEAR(norm(gx), norm(x), 1e-7);

        std::shared_ptr<std::vector<double> > outSer =
            std::make_shared<std::vector<double> >(gx->size(), 0.0);

        // apply serial Jacobian
        atmosLoc->applyMatrix(*gx, *outSer);

        Teuchos::RCP<Epetra_Vector> outPar = atmosPar->getState('C');
        outPar->PutScalar(0.0);

        // apply parallel Jacobian
        atmosPar->applyMatrix(*x, *outPar);

        // check results
        EXPECT_NEAR(norm(outSer), norm(outPar), 1e-7);

        INFO("TEST(Atmosphere, apply Jacobian): ||serial out|| = " << norm(outSer));
        INFO("TEST(Atmosphere, apply Jacobian): ||parall out|| = " << norm(outPar));

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

//---------------------------------------------------------------------
// We now switch on auxiliary unknowns and only test the parallel
// atmosphere
TEST(Atmosphere, Reinitialization)
{
    atmosphereParams->set("Auxiliary unknowns", 1);
    // create new object
    atmosPar = std::make_shared<Atmosphere>(comm, atmosphereParams);
}

//------------------------------------------------------------------
// FIXME We still need the ocean class to given us the maskstruct...
// this should definitely be factorized
TEST(Atmosphere, SetMasks)
{
    Teuchos::RCP<Teuchos::ParameterList> oceanParams =
        Utils::obtainParams("ocean_params.xml", "Ocean parameters");

    oceanParams->sublist("Belos Solver") =
        *Utils::obtainParams("solver_params.xml", "Solver parameters");
    Ocean ocean(comm, oceanParams);

    Utils::MaskStruct mask = ocean.getLandMask("mask_natl8");
    atmosPar->setLandMask(mask);

    Teuchos::RCP<Epetra_Vector> msi =
        Teuchos::rcp(new Epetra_Vector(*atmosPar->getDomain()->GetStandardSurfaceMap()));

    for (int i = 0; i < msi->MyLength() / 5; ++i)
        (*msi)[i] = 1;

    atmosPar->setSeaIceMask(msi);
}

//------------------------------------------------------------------
TEST(Atmosphere, numericalJacobian)
{
    atmosPar->getState('V')->Random();
    atmosPar->getState('V')->Scale(2.0);
    atmosPar->setPar("Combined Forcing", 0.1);

    atmosPar->computeRHS();
    atmosPar->computeJacobian();

    Teuchos::RCP<Epetra_CrsMatrix> atmosJac  = atmosPar->getJacobian();
    DUMPMATLAB("atmosJac", *atmosJac);

    // only do this test for small problems in serial
    int nmax = 2e3;

    if ( (comm->NumProc() == 1) &&
         (atmosPar->getState('V')->GlobalLength() < nmax) )
    {
        bool failed = false;
        try
        {
            INFO("compute njC");
            NumericalJacobian<std::shared_ptr<Atmosphere>,
                              Teuchos::RCP<Epetra_Vector> > numJac;
            
            numJac.setTolerance(1e-12);
            numJac.seth(1e-6);
            numJac.compute(atmosPar, atmosPar->getState('V'));
            numJac.print("atmosNumJac");

            NumericalJacobian<std::shared_ptr<Atmosphere>,
                              Teuchos::RCP<Epetra_Vector> >::CCS ccs;
            numJac.fillCCS(ccs);

            EXPECT_NE(ccs.beg.back(), 0);

            Teuchos::RCP<Epetra_Vector> x = atmosPar->getState('C');

            testEntries(atmosPar, ccs, x);

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
        std::cout << ("****Numerical Jacobian test cannot run in parallel****\n") ;
        INFO("****Numerical Jacobian test cannot run in parallel****");
    }

    if (atmosPar->getState('V')->GlobalLength() > nmax)
    {
        std::cout << ("****Numerical Jacobian test cannot run for this problem size****\n");
        INFO("****Numerical Jacobian test cannot run for this problem size****");
    }
}


//------------------------------------------------------------------
TEST(Atmosphere, Newton)
{
    Teuchos::RCP<Epetra_Vector> state = atmosPar->getState('V');
    state->PutScalar(0.0);
    atmosPar->setPar("Combined Forcing", 0.4);

    Teuchos::RCP<Epetra_Vector> b = atmosPar->getRHS('V');
    Teuchos::RCP<Epetra_Vector> x = atmosPar->getSolution('V');

    int maxit = 10;
    int niter = 0;
    for (; niter != maxit; ++niter)
    {
        INFO(" old ||b|| = " << Utils::norm(b) );
        atmosPar->computeRHS();
        INFO(" new ||b|| = " << Utils::norm(b) );
        atmosPar->computeJacobian();

        Teuchos::RCP<Epetra_Vector> r = atmosPar->getSolution('C');

        b->Scale(-1.0);
        r->PutScalar(0.0);

        INFO(" old ||x|| = " << Utils::norm(x) );
        atmosPar->solve(b);
        INFO(" new ||x|| = " << Utils::norm(x) );

        atmosPar->applyMatrix(*x, *r);
        r->Update(1.0, *b, -1.0);
        r->Scale(1./Utils::norm(b));

        EXPECT_NEAR(Utils::norm(r), 0, 1e-1);

        INFO("     ||r|| = " << Utils::norm(r) << std::endl);

        state->Update(1.0, *x, 1.0);
        b->Scale(-1.0);
        if (Utils::norm(b) < 1e-7)
            break;
    }

    EXPECT_LT(niter,maxit);
    if (niter < maxit)
    {
        INFO("Newton converged in " << niter << " iterations.");
    }
    else
    {
        WARNING("Newton did not converge", __FILE__, __LINE__);
    }

    Teuchos::RCP<Epetra_CrsMatrix> jac = atmosPar->getJacobian();
    Teuchos::RCP<Epetra_Vector>      B = atmosPar->getMassMat();

    DUMPMATLAB("atmos_jac", *jac);
    DUMP_VECTOR("atmos_B", *B);

    EXPECT_NEAR(Utils::norm(b), 0, 1e-7);
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
    atmosLoc  = std::shared_ptr<AtmosLocal>();
    atmosPar  = std::shared_ptr<Atmosphere>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
