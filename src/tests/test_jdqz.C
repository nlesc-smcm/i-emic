#include "TestDefinitions.H"

#include "ComplexVector.H"
#include "JDQZInterface.H"
#include "jdqz.hpp"

//------------------------------------------------------------------
using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
	std::shared_ptr<Ocean>         ocean;
	std::shared_ptr<Atmosphere>    atmos;
    std::shared_ptr<SeaIce>        seaice;
	std::shared_ptr<CoupledModel>  coupledModel;
    RCP<Epetra_Comm>               comm;

}

//------------------------------------------------------------------
class IEMIC : public testing::Environment
{
public:
	// constructor
	IEMIC()
		{
            // Create parameter object for Ocean
            RCP<Teuchos::ParameterList> oceanParams =
                Utils::obtainParams("ocean_params.xml", "Ocean parameters");

            oceanParams->sublist("Belos Solver") =
                *Utils::obtainParams("solver_params.xml", "Solver parameters");

            // Create parameter object for Atmosphere
            RCP<Teuchos::ParameterList> atmosphereParams =
                Utils::obtainParams("atmosphere_params.xml", "Atmosphere parameters");

            // Create parameter object for SeaIce
            RCP<Teuchos::ParameterList> seaIceParams =
                Utils::obtainParams("seaice_params.xml", "Sea ice parameters");

            // Create parameter object for CoupledModel
            RCP<Teuchos::ParameterList> coupledmodelParams =
                Utils::obtainParams("coupledmodel_params.xml", "CoupledModel parameters");

            INFO('\n' << "Overwriting:");
            // The Continuation and CoupledModel parameterlists overwrite settings
            Utils::overwriteParameters(oceanParams,        coupledmodelParams);
            Utils::overwriteParameters(atmosphereParams,   coupledmodelParams);
            Utils::overwriteParameters(seaIceParams,       coupledmodelParams);

            // Create models
 			ocean  = std::make_shared<Ocean>(comm, oceanParams);
			atmos  = std::make_shared<Atmosphere>(comm, atmosphereParams);
            seaice = std::make_shared<SeaIce>(comm, seaIceParams);

			coupledModel =
				std::make_shared<CoupledModel>(ocean,
                                               atmos,
                                               seaice,
                                               coupledmodelParams);
		}

	// destructor
	~IEMIC()
		{}
};

//------------------------------------------------------------------
TEST(JDQZ, CoupledContinuation)
{
        // Create Continuation
        RCP<Teuchos::ParameterList> continuationParams = rcp(new Teuchos::ParameterList);
        updateParametersFromXmlFile("continuation_params.xml", continuationParams.ptr());

        Continuation<std::shared_ptr<CoupledModel>>
            continuation(coupledModel, continuationParams);

        // We perform a continuation with the coupled model
        int status = continuation.run();
        EXPECT_EQ(status, 0);

        // Dump the matrices for checking with another code
        Teuchos::RCP<Epetra_CrsMatrix> oceanJac = ocean->getJacobian();
        Teuchos::RCP<Epetra_CrsMatrix> atmosJac = atmos->getJacobian();

        Teuchos::RCP<Epetra_Vector> oceanB = ocean->getMassMat();
        Teuchos::RCP<Epetra_Vector> atmosB = atmos->getMassMat();

        DUMPMATLAB("ocean_jac", *oceanJac);
        DUMP_VECTOR("ocean_B", *oceanB);
        DUMPMATLAB("atmos_jac", *atmosJac);
        DUMP_VECTOR("atmos_B", *atmosB);
}

//------------------------------------------------------------------
TEST(JDQZ, AtmosphereEigenvalues)
{
    bool failed = false;
    try
    {

        // Now we are going to calculate the eigenvalues of the Atmosphere Jacobian
        INFO("Creating ComplexVector...");
        Teuchos::RCP<Epetra_Vector> x(atmos->getSolution('C'));
        Teuchos::RCP<Epetra_Vector> y(atmos->getSolution('C'));
        x->PutScalar(0.0); y->PutScalar(0.0);

        // JDQZ needs complex arithmetic, so we create a ComplexVector
        ComplexVector<Epetra_Vector> z(*x, *y);
        ComplexVector<Epetra_Vector> residue(*x, *y);
        ComplexVector<Epetra_Vector> tmp(*x, *y);

        INFO("Building JDQZInterface...");
        JDQZInterface<std::shared_ptr<Atmosphere>,
                      ComplexVector<Epetra_Vector > >	matrix(atmos, z);

        INFO("Building JDQZ...");
        JDQZ<JDQZInterface<std::shared_ptr<Atmosphere>,
                           ComplexVector<Epetra_Vector > > > jdqz(matrix, z);

        INFO("Setting parameters...");
        std::map<std::string, double> list;
        list["Shift (real part)"]         = 0.0;
        list["Number of eigenvalues"]     = 5;
        list["Max size search space"]     = 35;
        list["Min size search space"]     = 10;
        list["Max JD iterations"]         = 500;
        list["Tracking parameter"]        = 1e-8;
        list["Criterion for Ritz values"] = 0;
        list["Linear solver"]             = 1;
        list["GMRES search space"]        = 20;
        list["Verbosity"]                 = 5;
        MyParameterList params(list);

        jdqz.setParameters(params);
        jdqz.printParameters();

        INFO("Starting JDQZ solve...");
        jdqz.solve();

        std::vector<ComplexVector<Epetra_Vector> > eivec = jdqz.getEigenVectors();
        std::vector<std::complex<double> >         alpha = jdqz.getAlpha();
        std::vector<std::complex<double> >         beta  = jdqz.getBeta();

        EXPECT_GT(jdqz.kmax(), 0);
        for (int j = 0; j != jdqz.kmax(); ++j)
        {
            residue.zero();
            tmp.zero();
            matrix.AMUL(eivec[j], residue);
            residue.scale(beta[j]);
            matrix.BMUL(eivec[j], tmp);
            residue.axpy(-alpha[j], tmp);
            std::cout << "alpha: " << std::setw(30) << alpha[j]
                      << " beta: " << std::setw(15) << beta[j];
            std::cout << " alpha) / beta: " << std::setw(30)
                      << alpha[j] / beta[j];
            std::cout << " residue: " << residue.norm() << std::endl;
            EXPECT_LT(residue.norm(), 1e-7);
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
TEST(JDQZ, OceanEigenvalues)
{
    bool failed = false;
    try
    {

        // Here we are going to calculate the eigenvalues of the Ocean Jacobian
        INFO("Creating ComplexVector...");
        Teuchos::RCP<Epetra_Vector> x(ocean->getSolution('C'));
        Teuchos::RCP<Epetra_Vector> y(ocean->getSolution('C'));
        x->PutScalar(0.0); y->PutScalar(0.0);

        // JDQZ needs complex arithmetic, so we create a ComplexVector
        ComplexVector<Epetra_Vector> z(*x, *y);
        ComplexVector<Epetra_Vector> residue(*x, *y);
        ComplexVector<Epetra_Vector> tmp(*x, *y);

        INFO("Building JDQZInterface...");
        JDQZInterface<std::shared_ptr<Ocean>,
                      ComplexVector<Epetra_Vector > >	matrix(ocean, z);

        INFO("Building JDQZ...");
        JDQZ<JDQZInterface<std::shared_ptr<Ocean>,
                           ComplexVector<Epetra_Vector > > > jdqz(matrix, z);

        INFO("Setting parameters...");
        std::map<std::string, double> list;
        list["Shift (real part)"]         = 0.0;
        list["Number of eigenvalues"]     = 5;
        list["Max size search space"]     = 35;
        list["Min size search space"]     = 10;
        list["Max JD iterations"]         = 500;
        list["Tracking parameter"]        = 1e-8;
        list["Criterion for Ritz values"] = 0;
        list["Linear solver"]             = 1;
        list["GMRES search space"]        = 20;
        list["Verbosity"]                 = 5;
        MyParameterList params(list);

        jdqz.setParameters(params);
        jdqz.printParameters();

        INFO("Starting JDQZ solve...");
        jdqz.solve();

        std::vector<ComplexVector<Epetra_Vector> > eivec = jdqz.getEigenVectors();
        std::vector<std::complex<double> >         alpha = jdqz.getAlpha();
        std::vector<std::complex<double> >         beta  = jdqz.getBeta();

        EXPECT_GT(jdqz.kmax(), 0);
        for (int j = 0; j != jdqz.kmax(); ++j)
        {
            residue.zero();
            tmp.zero();
            matrix.AMUL(eivec[j], residue);
            residue.scale(beta[j]);
            matrix.BMUL(eivec[j], tmp);
            residue.axpy(-alpha[j], tmp);
            std::cout << "alpha: " << std::setw(30) << alpha[j]
                      << " beta: " << std::setw(15) << beta[j];
            std::cout << " alpha) / beta: " << std::setw(12)
                      << alpha[j] / beta[j];
            std::cout << " residue: " << residue.norm() << std::endl;
            EXPECT_LT(residue.norm(), 1e-7);
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

	// Get rid of possibly parallel objects:
	ocean        = std::shared_ptr<Ocean>();
	atmos        = std::shared_ptr<Atmosphere>();
	coupledModel = std::shared_ptr<CoupledModel>();

	comm->Barrier();
	std::cout << "TEST exit code proc #" << comm->MyPID()
			  << " " << out << std::endl;

	MPI_Finalize();
	return out;
}
