#include "TestDefinitions.H"

#include "AMS.H"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
Teuchos::RCP<Epetra_Comm>  comm; 
}

class TestModel
{
    Teuchos::RCP<Epetra_Map> map_;
    Teuchos::RCP<Epetra_Vector> rhs_;
    Teuchos::RCP<Epetra_Vector> state_;
    Teuchos::RCP<Epetra_Vector> diagB_;

    Teuchos::RCP<Epetra_CrsMatrix> frc_;
public:
    using Vector = Epetra_Vector;
    using VectorPtr = Teuchos::RCP<Vector>;

    TestModel(Teuchos::RCP<Epetra_Map> map)
        :
        map_(map)
        {
            rhs_ = Teuchos::rcp(new Epetra_Vector(*map_));
            state_ = Teuchos::rcp(new Epetra_Vector(*map_));

            std::vector<double> values(2, 1);
            diagB_ = Teuchos::rcp(new Epetra_Vector(Copy, *map_, &values[0]));
        }

    Teuchos::RCP<Epetra_Comm> Comm()
        {
            return comm;
        }

    void computeRHS()
        {
            std::vector<double> values(2);
            values[0] = (*state_)[0] - (*state_)[0] * (*state_)[0] * (*state_)[0];
            values[1] = -2 * (*state_)[1];
            rhs_ = Teuchos::rcp(new Epetra_Vector(Copy, *map_, &values[0]));
        }

    void computeJacobian() {}

    void computeForcing()
        {
            frc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *map_, 2));
            double val = 1.0;
            for (int idx = 0; idx < 2; idx++)
                frc_->InsertGlobalValues(idx, 1, &val, &idx);
            frc_->FillComplete();
        }

    Teuchos::RCP<Epetra_CrsMatrix> getForcing() {return frc_;}

    Teuchos::RCP<Epetra_Vector> getVector(char mode, RCP<Epetra_Vector> vec)
        {
            if (mode == 'C') // copy
            {
                RCP<Epetra_Vector> copy = rcp(new Epetra_Vector(*vec));
                return copy;
            }
            else if (mode == 'V') // view
            {
                return vec;
            }
            else
            {
                WARNING("Invalid mode", __FILE__, __LINE__);
                return Teuchos::null;
            }
        }

    Teuchos::RCP<Epetra_Vector> getState(char mode = 'C')
        {
            return getVector(mode, state_);
        }

    Teuchos::RCP<Epetra_Vector> getRHS(char mode = 'C')
        {
            return getVector(mode, rhs_);
        }

    Teuchos::RCP<Epetra_Vector> getMassMat(char mode = 'C')
        {
            return getVector(mode, diagB_);
        }

    Teuchos::RCP<Epetra_Vector> getSolution(char mode = 'C') { return Teuchos::null; }
    void setShift(double shift) {}
    void applyMatrix(Epetra_MultiVector const &v, Epetra_MultiVector &out) {}
    void solve(Teuchos::RCP<Epetra_MultiVector> rhs = Teuchos::null) {}
};

AMS<Teuchos::RCP<TestModel> > createDoubleWell(Teuchos::RCP<Teuchos::ParameterList> params)
{
    Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(2, 0, *comm));
    Teuchos::RCP<TestModel> model = Teuchos::rcp(new TestModel(map));

    std::vector<double> values(2);

    values[0] = -1;
    values[1] = 0;
    Teuchos::RCP<Epetra_Vector> sol1 = Teuchos::rcp(new Epetra_Vector(Copy, *map, &values[0]));

    values[0] = 1;
    values[1] = 0;
    Teuchos::RCP<Epetra_Vector> sol2 = Teuchos::rcp(new Epetra_Vector(Copy, *map, &values[0]));

    values[0] = 0;
    values[1] = 0;
    Teuchos::RCP<Epetra_Vector> sol3 = Teuchos::rcp(new Epetra_Vector(Copy, *map, &values[0]));

    AMS<Teuchos::RCP<TestModel> > ams(model, params, sol1, sol2, sol3);

    return ams;
}

//------------------------------------------------------------------
TEST(AMS, Restart)
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    params->set("theta", 0.0);
    params->set("sigma", 1.0);
    params->set("dof", 1);
    params->set("var", 0);
    params->set("noise seed", 1);
    params->set("ams seed", 2);

    params->set("time step", 0.01);
    params->set("A distance", 0.05);
    params->set("B distance", 0.05);
    params->set("C distance", 0.1);
    params->set("number of experiments", 100);
    params->set("maximum iterations", 100);
    params->set("write file", "out_data.h5");

    testing::internal::CaptureStdout();

    auto ams = createDoubleWell(params);
    ams.run();

    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_NE(output.find("Initialization"), std::string::npos);
    EXPECT_NE(output.find("AMS: 100"), std::string::npos);
    EXPECT_EQ(output.find("AMS: 101"), std::string::npos);

    params->set("read file", "out_data.h5");
    params->set("maximum iterations", 10000);

    testing::internal::CaptureStdout();

    auto ams2 = createDoubleWell(params);
    ams2.run();

    std::string output2 = testing::internal::GetCapturedStdout();

    EXPECT_EQ(output2.find("Initialization"), std::string::npos);
    EXPECT_EQ(output2.find("AMS: 100"), std::string::npos);
    EXPECT_NE(output2.find("AMS: 200"), std::string::npos);
    EXPECT_EQ(output2.find("AMS: 500"), std::string::npos);
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

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
