#include "TestDefinitions.H"

#include "NumericalJacobian.H"
#include "Ocean.H"
#include "Atmosphere.H"
#include "SeaIce.H"
#include "CoupledModel.H"

#include "THCMdefs.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Epetra_Comm>      comm;

    std::shared_ptr<Ocean>         ocean;
    std::shared_ptr<Atmosphere>    atmos;
    std::shared_ptr<SeaIce>        seaice;
    std::shared_ptr<CoupledModel>  coupledModel;
    std::vector<Teuchos::RCP<Teuchos::ParameterList> > params;
    enum Ident { OCEAN, ATMOS, SEAICE, COUPLED, CONT};
}

//------------------------------------------------------------------
TEST(ParameterLists, Initialization)
{
    bool failed = false;
    try
    {
        std::vector<std::string> files = {"ocean_params.xml",
                                          "atmosphere_params.xml",
                                          "seaice_params.xml",
                                          "coupledmodel_params.xml",
                                          "continuation_params.xml"};

        std::vector<std::string> names = {"Ocean parameters",
                                          "Atmosphere parameters",
                                          "Sea ice parameters",
                                          "CoupledModel parameters",
                                          "Continuation parameters"};

        for (int i = 0; i != (int) files.size(); ++i)
            params.push_back(Utils::obtainParams(files[i], names[i]));

        params[OCEAN]->sublist("Belos Solver") =
            *Utils::obtainParams("solver_params.xml", "Solver parameters");

        INFO('\n' << "Overwriting:");

        // Allow dominant parameterlists. Note that this trick uses a
        // 'flattened' hierarchy. The Continuation and CoupledModel
        // parameterlists are allowed to overwrite settings.
        Utils::overwriteParameters(params[OCEAN],  params[COUPLED]);
        Utils::overwriteParameters(params[ATMOS],  params[COUPLED]);
        Utils::overwriteParameters(params[SEAICE], params[COUPLED]);

        Utils::overwriteParameters(params[OCEAN],  params[CONT]);
        Utils::overwriteParameters(params[ATMOS],  params[CONT]);
        Utils::overwriteParameters(params[SEAICE], params[CONT]);

        Utils::overwriteParameters(params[COUPLED], params[CONT]);
        INFO('\n');
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed,false);
}

//------------------------------------------------------------------
TEST(Ocean, Initialization)
{
    bool failed = false;
    try
    {
        // Create parallel Ocean
        ocean = std::make_shared<Ocean>(comm, params[OCEAN]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Atmosphere, Initialization)
{
    bool failed = false;
    try
    {
        // Create atmosphere
        atmos = std::make_shared<Atmosphere>(comm, params[ATMOS]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(SeaIce, Initialization)
{
    bool failed = false;
    try
    {
        // Create atmosphere
        seaice = std::make_shared<SeaIce>(comm, params[SEAICE]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, Initialization)
{
    bool failed = false;
    try
    {
        // Create coupledmodel
        coupledModel = std::make_shared<CoupledModel>(ocean,
                                                      atmos,
                                                      seaice,
                                                      params[COUPLED]);
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, inspectState)
{
    bool failed = false;
    try
    {
        std::shared_ptr<Combined_MultiVec> state = coupledModel->getState('V');
        state->Random();
        size_t stateSize = state->Size();

        std::vector<int> lengthsGlobal(stateSize);
        std::vector<int> lengthsLocal(stateSize);
        int sumGlobal = 0, sumLocal = 0;
        for (size_t i = 0; i != stateSize; ++i)
        {
            lengthsGlobal[i] = (*state)(i)->GlobalLength();
            lengthsLocal[i] = (*state)(i)->MyLength();
            sumGlobal += lengthsGlobal[i];
            sumLocal  += lengthsLocal[i];
        }

        EXPECT_EQ(sumGlobal, state->GlobalLength());
        EXPECT_EQ(sumLocal,  state->MyLength());

        std::vector<double> norms(stateSize);
        double sumNorms = 0.0;
        for (size_t i = 0; i != stateSize; ++i)
        {
            norms[i] = Utils::norm((*state)(i));
            sumNorms += pow(norms[i], 2);
        }

        EXPECT_NEAR(Utils::norm(state), sqrt(sumNorms), 1e-7);

        std::shared_ptr<Combined_MultiVec> x = coupledModel->getState('C');
        x->PutScalar(1.0);

        // Check whether the indexing works correctly
        for (int i = 0; i != x->MyLength(); ++i)
        {
            EXPECT_EQ((*x)[i], 1);
        }

    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, MassMatrix)
{
    Combined_MultiVec v   = *coupledModel->getState('C');
    Combined_MultiVec out = *coupledModel->getState('C');

    v.PutScalar(1.0);
    out.Random();

    coupledModel->applyMassMat(v, out);

    EXPECT_GT(Utils::norm(out), 0.0);

    Teuchos::RCP<Epetra_MultiVector> oceanB = out(0);
    Teuchos::RCP<Epetra_MultiVector> atmosB = out(1);

    int n = ocean->getNdim();
    int m = ocean->getMdim();

    std::ofstream file;
    file.open("massmat");
    file << out;
    file.close();

    int rowintconQ = atmos->getRowIntCon();
    int rowintconS = ocean->getRowIntCon();

    // ocean integral condition
    if (oceanB->Map().MyGID(rowintconS))
    {
        int lid = oceanB->Map().LID(rowintconS);
        EXPECT_EQ(0.0, (*oceanB)[0][lid]);
    }

    // atmos integral condition
    if (atmosB->Map().MyGID(rowintconQ))
    {
        int lid = atmosB->Map().LID(rowintconQ);
        EXPECT_EQ((*atmosB)[0][lid], 0.0);
    }

    // auxiliary integral equations
    int atmosLast = ATMOS_NUN_ * m * n - 1;
    if (atmosB->GlobalLength() > ATMOS_NUN_ * m * n)
    {
        int aux = params[ATMOS]->get("Auxiliary unknowns", 0);
        EXPECT_NE(aux, 0);
        for (int i = 1; i <= aux; ++i)
        {
            if (atmosB->Map().MyGID(atmosLast+i))
            {
                int lid = atmosB->Map().LID(atmosLast+i);
                EXPECT_EQ((*atmosB)[0][lid], 0.0);
            }
        }
    }
}

//------------------------------------------------------------------
TEST(CoupledModel, computeJacobian)
{
    bool failed = false;
    try
    {
        // Set a small parameter
        coupledModel->setPar("Combined Forcing", 0.1);

        // put values in state
        coupledModel->getState('V')->Random();
        coupledModel->getState('V')->Scale(1e-4);

        // To check derivatives with sea ice mask components, we make
        // sure that the sea ice mask has active and inactive parts:
        Teuchos::RCP<Epetra_Vector> seaice_state = seaice->getState('V');
        Teuchos::RCP<Epetra_Vector> seaice_mask  = seaice->interfaceM();
        seaice_mask->PutScalar(0.0);
        for (int i = 0; i < seaice_mask->MyLength()/4; ++i)
            (*seaice_mask)[i] = 1.0;

        Epetra_Import importM(seaice_mask->Map(), seaice_state->Map());
        seaice_state->Export(*seaice_mask, importM, Zero);

        coupledModel->computeJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> atmosJac  = atmos->getJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> oceanJac  = ocean->getJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> seaiceJac = seaice->getJacobian();

        Utils::print(atmosJac,  "atmosJac");
        Utils::print(oceanJac,  "oceanJac");
        Utils::print(seaiceJac, "seaiceJac");

        Utils::print(&atmosJac->ColMap(), "atmosJacColMap");
        Utils::print(&atmosJac->DomainMap(), "atmosJacDomainMap");

        Utils::print(&seaiceJac->ColMap(), "seaiceJacColMap");
        Utils::print(&seaiceJac->DomainMap(), "seaiceJacDomainMap");

        coupledModel->dumpBlocks();
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(CoupledModel, numericalJacobian)
{
    // only do this test for small problems in serial
    int nmax = 2e3;

    if ( (comm->NumProc() == 1) &&
         (coupledModel->getState('V')->GlobalLength() < nmax) )
    {
        bool failed = false;
        try
        {
            INFO("compute njC");

            NumericalJacobian<std::shared_ptr<CoupledModel>,
                              std::shared_ptr<Combined_MultiVec> > njC;

            njC.setTolerance(1e-9);
            njC.seth(1e-4);
            njC.compute(coupledModel, coupledModel->getState('V'));

            std::string fnameJnC("JnC");

            INFO(" Printing Numerical Jacobian " << fnameJnC);

            njC.print(fnameJnC);

            // test individual elements
            NumericalJacobian<std::shared_ptr<CoupledModel>,
                              std::shared_ptr<Combined_MultiVec> >::CCS ccs;
            njC.fillCCS(ccs);

            EXPECT_NE(ccs.beg.back(), 0);

            std::shared_ptr<Combined_MultiVec> x = coupledModel->getState('C');

            testEntries(coupledModel, ccs, x);

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

    if (coupledModel->getState('V')->GlobalLength() > nmax)
    {
        std::cout << ("****Numerical Jacobian test cannot run for this problem size****\n");
        INFO("****Numerical Jacobian test cannot run for this problem size****");
    }

}

//------------------------------------------------------------------
// We need this information from THCM
extern "C" _SUBROUTINE_(getdeps)(double*, double*, double*,
                                 double*, double*, double*,
                                 double *);

TEST(CoupledModel, applyMatrix)
{
    bool failed = false;

    try
    {
        // reset model
        coupledModel->getState('V')->PutScalar(0.0);
        coupledModel->setPar("Combined Forcing", 0.0);
        coupledModel->computeRHS(); // synchronize

        // get mask
        Utils::MaskStruct mask = ocean->getLandMask();
        atmos->setLandMask(mask);

        // set values
        coupledModel->getState('V')->Random();
        double parValue = 0.1;
        coupledModel->setPar("Combined Forcing", parValue);

        std::shared_ptr<Combined_MultiVec> x = coupledModel->getState('C');
        std::shared_ptr<Combined_MultiVec> y = coupledModel->getState('C');

        double normIn = Utils::norm(y);

        coupledModel->applyMatrix(*x, *y);

        double normOut = Utils::norm(y);

        EXPECT_NE(normIn, normOut);

        // Test the coupling blocks in the matrix separately

        CouplingBlock<std::shared_ptr<Ocean>,
                      std::shared_ptr<Atmosphere> > C12(ocean, atmos);

        CouplingBlock<std::shared_ptr<Atmosphere>,
                      std::shared_ptr<Ocean> > C21(atmos, ocean);

        Teuchos::RCP<Epetra_MultiVector> oceanVec = (*x)(0);
        Teuchos::RCP<Epetra_MultiVector> atmosVec = (*x)(1);

        int n = ocean->getNdim();
        int m = ocean->getMdim();
        int l = ocean->getLdim();
        double value[3] = {1.234, 2.12, -23.5};

        // Test atmos -> ocean coupling
        int ii = n-2;
        int jj = m-2;

        // Get shortwave radiative heat flux dependence
        Teuchos::RCP<Epetra_MultiVector> suno =
            Utils::AllGather(*ocean->getSunO());

        double So = (*(*suno)(0))[jj*n + ii];

        Atmosphere::CommPars atmosPars;
        atmos->getCommPars(atmosPars);

        for (int v = 0; v != 3; ++v)
        {
            atmosVec->PutScalar(value[v]);
            oceanVec->PutScalar(0.0);

            C12.applyMatrix(*atmosVec, *oceanVec);

            // Get ocean parameters
            double Ooa, Os, nus, eta, lvsc, qdim, pQSnd;
            FNAME(getdeps)(&Ooa, &Os, &nus, &eta, &lvsc, &qdim, &pQSnd);

            // Test center surface element (temperature)
            int surfbT = FIND_ROW2(_NUN_, n, m, l, ii, jj, l-1, TT);

            INFO( " surface element TT " << surfbT );

            // If this point is on land the test does not make any sense
            bool onLand = ((*mask.global_surface)[jj*n+ii] > 0) ? true : false;

            if (onLand)
            {
                ERROR(" ** applyMatrix is testing a land point, not useful... **",
                      __FILE__, __LINE__);
            }
            else
            {
                int lid;
                double surfval;
                if (oceanVec->Map().MyGID(surfbT) && !onLand)
                {
                    lid = oceanVec->Map().LID(surfbT);
                    surfval = (*oceanVec)[0][lid];
                    EXPECT_NEAR(-(Ooa + lvsc * eta * qdim - parValue * So * atmosPars.da ) * value[v],
                                surfval , 1e-7);
                }

                // Test first surface element (salinity)
                int surfbS = FIND_ROW2(_NUN_, n, m, l, ii, jj, l-1, SS);
                INFO( " surface element SS " << surfbS );

                if (oceanVec->Map().MyGID(surfbS))
                {
                    lid = oceanVec->Map().LID(surfbS);
                    surfval = (*oceanVec)[0][lid];
                    EXPECT_NEAR( (nus + nus ) * value[v], surfval , 1e-7);
                }
            }
        }

        // qdim, nuq, eta, dqso, dqdt, Eo0;
        Atmosphere::CommPars pars;
        atmos->getCommPars(pars);

        // Test ocean -> atmos coupling
        for (int v = 0; v != 3; ++v)
        {
            // Put values in temperature domain points. Atmosphere
            // only depends on sst.
            oceanVec->PutScalar(0.0);
            for (int i = TT-1; i < oceanVec->MyLength(); i += _NUN_)
            {
                (*oceanVec)[0][i] = value[v];
            }

            // Initialize atmosphere vector
            CHECK_ZERO(atmosVec->PutScalar(0.0));

            // Perform matvec with coupling block
            C21.applyMatrix(*oceanVec, *atmosVec);

            // assume single atmosphere layer
            l = 1;

            // test atmos temperature point
            int surfbT = FIND_ROW_ATMOS0(ATMOS_NUN_, n, m, l, ii, jj, l-1, ATMOS_TT_);

            INFO( " surface element TT " << surfbT );

            // If this point is on land the test does not make any sense
            bool onLand = ((*mask.global_surface)[jj*n+ii] > 0) ? true : false;

            if (onLand)
            {
                ERROR(" ** applyMatrix is testing a land point, not useful... **",
                      __FILE__, __LINE__);
            }
            else
            {
                int lid;
                double surfval;
                if (atmosVec->Map().MyGID(surfbT))
                {
                    lid = atmosVec->Map().LID(surfbT);
                    surfval = (*atmosVec)[0][lid];
                    EXPECT_NEAR( value[v], surfval, 1e-7);
                }

                // Test atmos humidity point
                int surfbQ = FIND_ROW_ATMOS0(ATMOS_NUN_, n, m, l, ii, jj, l-1, ATMOS_QQ_);

                if (atmosVec->Map().MyGID(surfbQ))
                {
                    lid = atmosVec->Map().LID(surfbQ);
                    surfval = (*atmosVec)[0][lid];
                    EXPECT_NEAR( pars.dqdt * value[v], surfval, 1e-7);
                }

                // Test final element in range (results from sst integral)

                // First check if we have auxiliary unknowns:
                int aux;
                aux = params[ATMOS]->get("Auxiliary unknowns", 0);

                if (atmosVec->GlobalLength() > ATMOS_NUN_ * m * n * l)
                {
                    EXPECT_GT(aux, 0);
                }
                else
                {
                    EXPECT_EQ(aux, 0);
                }

                // If we have auxiliary unknowns, these govern precipitation
                if (aux > 0)
                {
                    Teuchos::RCP<Epetra_Vector> precipintco = atmos->getPIntCoeff();
                    double totalArea;
                    precipintco->Norm1(&totalArea);

                    std::cout << std::endl;
                    std::cout << " totalArea: " << totalArea << std::endl;

                    Teuchos::RCP<Epetra_Vector> ones =
                        Teuchos::rcp(new Epetra_Vector(*precipintco));
                    ones->PutScalar(1.0);

                    double sstInt = Utils::dot(precipintco, ones);
                    std::cout << " sstInt: " << sstInt << std::endl;

                    double intval = sstInt * (1.0 / totalArea) *
                        (pars.tdim / pars.qdim) * pars.dqso;

                    std::cout << " intval: " << intval << std::endl;

                    int last = FIND_ROW_ATMOS0(ATMOS_NUN_, n, m, l, n-1 , m-1, l-1, ATMOS_NUN_);
                    std::cout << " last  : " << last << std::endl;

                    if (atmosVec->Map().MyGID(last + 1))
                    {
                        lid = atmosVec->Map().LID(last + 1);
                        surfval = (*atmosVec)[0][lid];
                        std::cout << "  surfval        : " << surfval << std::endl;
                        std::cout << "  value[v]       : " << value[v] << std::endl;
                        std::cout << "  intval*value[v]: " << intval*value[v] << std::endl;

                        EXPECT_NEAR( intval * value[v], surfval, 1e-7);
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}
//------------------------------------------------------------------
TEST(CoupledModel, Precipitation)
{
    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    stateV->Random();

    // copy construct
    std::shared_ptr<Combined_MultiVec> b =
        std::make_shared<Combined_MultiVec>(*stateV);

    coupledModel->computeJacobian();
    coupledModel->applyMatrix(*stateV, *b);

    Teuchos::RCP<Epetra_MultiVector> atmb = (*b)(1);
    Teuchos::RCP<Epetra_Vector> P = atmos->interfaceP();

    int numMyElements = P->Map().NumMyElements();
    double Pval = 0.0;
    for (int i = 0; i != numMyElements; ++i)
    {
        if (std::abs((*P)[i]) > 0)
        {
            Pval = (*P)[i];
            break;
        }
    }

    std::cout << "Pval = " << Pval << std::endl;

    int n = ocean->getNdim();
    int m = ocean->getMdim();
    // assume single atmosphere layer
    int l = 1;

    int last = FIND_ROW_ATMOS0(ATMOS_NUN_, n, m, l, n-1 , m-1, l-1, ATMOS_NUN_);
    int lid;
    double val;
    if (atmb->Map().MyGID(last+1))
    {
        lid = atmb->Map().LID(last+1);
        val = (*atmb)[0][lid];
        std::cout << "b[last+1] = " << val << std::endl;

    }

    // FIXME TODO
    // EXPECT_EQ(
}


//------------------------------------------------------------------
TEST(CoupledModel, View)
{
    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    std::shared_ptr<Combined_MultiVec> rhsV =
        coupledModel->getRHS('V');


    stateV->PutScalar(0.1);
    coupledModel->computeRHS();

    double norm1 = Utils::norm(rhsV);

    stateV->PutScalar(1.0);
    coupledModel->computeRHS();

    double norm2 = Utils::norm(rhsV);

    EXPECT_NE(norm1, norm2);
}


//------------------------------------------------------------------
TEST(CoupledModel, Synchronization)
{
    bool failed = false;

    std::shared_ptr<Combined_MultiVec> stateV =
        coupledModel->getState('V');

    try
    {
        (*stateV)(0)->PutScalar(1.234);
        (*stateV)(1)->PutScalar(2.345);

        // Set a small parameter
        coupledModel->setPar("Combined Forcing", 0.01);

        // At RHS computation the coupledModel synchronizes the states
        coupledModel->computeRHS();
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    // Obtain atmosphere temperature existing in ocean
    Teuchos::RCP<Epetra_Vector> oceanAtmosT  = ocean->getLocalAtmosT();

    double maxValue, minValue, meanValue;
    oceanAtmosT->MaxValue(&maxValue);
    EXPECT_EQ(maxValue, 2.345);

    oceanAtmosT->MinValue(&minValue);
    EXPECT_EQ(minValue, 2.345);

    if ( oceanAtmosT->Map().UniqueGIDs() )
    {
        oceanAtmosT->MeanValue(&meanValue);
        EXPECT_NEAR(meanValue, 2.345, 1e-7);
    }
    else
        INFO(" oceanAtmosT (overl.) GID's not unique, which is fine in parallel.");

    Teuchos::RCP<Epetra_Vector> atmosOceanT = atmos->getLocalSST();

    Utils::print(atmosOceanT, "atmosOceanT" + std::to_string(comm->MyPID()));

    atmosOceanT->MaxValue(&maxValue);
    EXPECT_EQ(maxValue, 1.234);
    atmosOceanT->MinValue(&minValue);
    EXPECT_EQ(minValue, 1.234);

    if (atmosOceanT->Map().UniqueGIDs())
    {
        atmosOceanT->MeanValue(&meanValue);
        EXPECT_NEAR(meanValue, 1.234, 1e-7);
    }
    else
        INFO(" atmosOceanT (overl.) GID's not unique, which is fine in parallel.");

    // Randomize combined state
    stateV->Random();

    // Make sure that the sea ice mask is inactive
    Teuchos::RCP<Epetra_Vector> seaice_state = seaice->getState('V');
    Teuchos::RCP<Epetra_Vector> seaice_mask  = seaice->interfaceM();
    seaice_mask->PutScalar(0.0);
    Epetra_Import importM(seaice_mask->Map(), seaice_state->Map());
    seaice_state->Export(*seaice_mask, importM, Zero);

    try
    {
        // At RHS computation the coupledModel synchronizes the states
        coupledModel->computeRHS();
    }
    catch (...)
    {
        failed = true;
        throw;
    }
    EXPECT_EQ(failed, false);

    // Evaporation is calculated simultaneously in Ocean and in
    // Atmosphere during the RHS computation above. Here we check
    // whether they return the same fully dimensional norm.
    Teuchos::RCP<Epetra_Vector> atmosE = atmos->interfaceE();
    Teuchos::RCP<Epetra_Vector> oceanE = ocean->interfaceE();

    double nrmAtmosE = Utils::norm(atmosE);
    double nrmOceanE = Utils::norm(oceanE);
    EXPECT_NEAR(nrmAtmosE, nrmOceanE, 1e-7);

    // Precipitation should have been calculated in the atmosphere model
    Teuchos::RCP<Epetra_Vector> oceanAtmosP = ocean->getLocalAtmosP();

#ifdef GNU
    Utils::print(oceanAtmosT,  "oceanAtmosT" +
                 std::to_string(oceanAtmosT->Map().Comm().MyPID()) + ".txt");
    Utils::print(oceanAtmosP,  "oceanAtmosP" +
                 std::to_string(oceanAtmosP->Map().Comm().MyPID()) + ".txt");
#endif

    // There should be something in there
    oceanAtmosP->MaxValue(&maxValue);
    oceanAtmosP->MinValue(&minValue);
    EXPECT_GT(std::max(std::abs(maxValue), std::abs(minValue)), 0.0);
}

//------------------------------------------------------------------
TEST(CoupledModel, Synchronization2)
{
    // reset model
    coupledModel->getState('V')->PutScalar(0.0);
    coupledModel->getSolution('V')->PutScalar(0.0);
    coupledModel->setPar("Combined Forcing", 0.0);

    // perturb model
    coupledModel->setPar("Combined Forcing", 0.01);
    coupledModel->computeRHS();
    coupledModel->computeJacobian();

    // solve with perturbed matrix and rhs
    std::shared_ptr<Combined_MultiVec> b = coupledModel->getRHS('C');
    coupledModel->solve(b);
    std::shared_ptr<Combined_MultiVec> x = coupledModel->getSolution('V');

    // compute norms of solution
    int numModels = b->Size();
    std::vector<double> norms(numModels);
    for (int i = 0; i != numModels; ++i)
    {
        norms[i] = Utils::norm((*x)(i));
        std::cout << norms[i] << " ";
    }
    std::cout << std::endl;

    // perturb last vector in combined_multivec
    (*b)(b->Size()-1)->Scale(1.1);
    coupledModel->solve(b);

    // check whether the solve synchronizes ok
    double norm;
    for (int i = 0; i != numModels; ++i)
    {
        norm = Utils::norm((*x)(i));
        std::cout << norm << " ";
        EXPECT_NE(norms[i], norm);
    }
    std::cout << std::endl;
}

//------------------------------------------------------------------
// Test hashing functions of Combined_MultiVec and Utils
TEST(CoupledModel, Hashing)
{
    std::shared_ptr<Combined_MultiVec> x  = coupledModel->getState('C');
    x->Random();
    std::size_t hash1 = x->hash(); INFO(" hash1 = " << hash1);
    std::size_t hash2 = x->hash(); INFO(" hash2 = " << hash2);
    EXPECT_EQ(hash1, hash2);

    x->Scale(1.0001);

    std::size_t hash3 = x->hash(); INFO(" hash3 = " << hash3);

    EXPECT_NE(hash2, hash3);

    x->PutScalar(1.00001);
    std::size_t hash4 = x->hash(); INFO(" hash4 = " << hash4);

    EXPECT_NE(hash3, hash4);
}

//------------------------------------------------------------------
TEST(CoupledModel, Solve)
{
    bool failed = false;
    try
    {
        // perturb model
        double par = coupledModel->getPar("Combined Forcing");
        coupledModel->setPar("Combined Forcing", par + 0.01);
        coupledModel->computeRHS();
        coupledModel->computeJacobian();

        std::shared_ptr<Combined_MultiVec> x = coupledModel->getSolution('C');
        std::shared_ptr<Combined_MultiVec> b = coupledModel->getRHS('C');
        std::shared_ptr<Combined_MultiVec> c = coupledModel->getRHS('C');

        coupledModel->solve(b);
        std::shared_ptr<Combined_MultiVec> sol = coupledModel->getSolution('C');

        EXPECT_EQ(b->hash(), c->hash());
        coupledModel->applyMatrix(*sol, *c);
        EXPECT_NE(b->hash(), c->hash());
        EXPECT_NEAR(Utils::norm(b), Utils::norm(c), 1e-3);
    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
// Here we are testing the implementation of the integral condition.
// The result of a matrix vector product of the Jacobian with the
// state vector should equal the value of the rhs in the integral
// equation row <rowintcon>.
TEST(CoupledModel, IntegralCondition)
{
    // Randomize coupledModel state
    coupledModel->getState('V')->Random();

    std::shared_ptr<Combined_MultiVec> x = coupledModel->getState('C');
    std::shared_ptr<Combined_MultiVec> b = coupledModel->getState('C');

    // Compute rhs
    coupledModel->computeRHS();
    std::shared_ptr<Combined_MultiVec> F = coupledModel->getRHS('C');

    // Compute Jacobian and matrix vector product
    b->PutScalar(0.0);
    coupledModel->computeJacobian();
    coupledModel->applyMatrix(*x, *b);

    Teuchos::RCP<Epetra_MultiVector> oceanB = (*b)(0);
    Teuchos::RCP<Epetra_MultiVector> oceanF = (*F)(0);

    int rowintcon = ocean->getRowIntCon();
    int lid = -1;

    double valueFloc = 0.0;
    double valueBloc = 0.0;
    double valueF    = 0.0;
    double valueB    = 0.0;

    if (oceanF->Map().MyGID(rowintcon))
    {
        lid = oceanF->Map().LID(rowintcon);
        valueFloc = (*oceanF)[0][lid];
    }

    if (oceanB->Map().MyGID(rowintcon))
    {
        lid = oceanB->Map().LID(rowintcon);
        valueBloc = (*oceanB)[0][lid];
    }

    comm->SumAll(&valueFloc, &valueF, 1);
    comm->SumAll(&valueBloc, &valueB, 1);

    INFO("\n *-* Integral condition           *-* ");
    INFO("     RHS[oceanLast] = " << valueF);
    INFO("   (J*x)[oceanLast] = " << valueB);
    INFO(" *-*                              *-* \n");

    EXPECT_NEAR(valueF, valueB, 1e-12);
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
    ocean        = std::shared_ptr<Ocean>();
    atmos        = std::shared_ptr<Atmosphere>();
    seaice       = std::shared_ptr<SeaIce>();
    coupledModel = std::shared_ptr<CoupledModel>();

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
