#include "TestDefinitions.H"

#include "AMS.H"

#include <stdio.h>

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

template<typename T>
void set_parameter(Teuchos::RCP<Teuchos::ParameterList> &params,
                   std::string const &name, T value)
{
    if (!params->isParameter(name))
        params->set(name, value);
}

void restart_test(Teuchos::RCP<Teuchos::ParameterList> params)
{
    set_parameter(params, "theta", 0.0);
    set_parameter(params, "sigma", 1.0);
    set_parameter(params, "dof", 1);
    set_parameter(params, "var", 0);
    set_parameter(params, "noise seed", 1);
    set_parameter(params, "ams seed", 2);
    set_parameter(params, "method", "TAMS");

    set_parameter(params, "time step", 0.01);
    set_parameter(params, "maximum time", 2.0);
    set_parameter(params, "B distance", 0.05);
    set_parameter(params, "number of experiments", 100);
    set_parameter(params, "maximum iterations", 100);
    set_parameter(params, "write file", "out_data.h5");

    int maxit = params->get("maximum iterations", -1);
    int write_time_steps = params->get("write time steps", -1);

    remove("out_data.h5");

    testing::internal::CaptureStdout();

    auto ams = createDoubleWell(params);
    ams.run();

    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_NE(output.find("Initialization"), std::string::npos);
    if (maxit == 100 && write_time_steps < 0)
    {
        EXPECT_NE(output.find("TAMS: 100"), std::string::npos);
        EXPECT_EQ(output.find("TAMS: 101"), std::string::npos);
    }

    if (maxit == 100 && write_time_steps > 0)
    {
        EXPECT_NE(output.find("TAMS: 80"), std::string::npos);
    }

    params->set("read file", "out_data.h5");
    params->set("write file", "");
    params->set("maximum iterations", 10000);

    testing::internal::CaptureStdout();

    auto ams2 = createDoubleWell(params);
    ams2.run();

    std::string output2 = testing::internal::GetCapturedStdout();

    EXPECT_EQ(output2.find("Initialization"), std::string::npos);
    if (maxit == 100 && write_time_steps < 0)
    {
        EXPECT_EQ(output2.find("TAMS: 100"), std::string::npos);
    }
    else if (maxit < 100)
    {
        EXPECT_NE(output2.find("TAMS: 100"), std::string::npos);
    }
    EXPECT_NE(output2.find("TAMS: 101"), std::string::npos);
    EXPECT_EQ(output2.find("TAMS: 500"), std::string::npos);
}

//------------------------------------------------------------------
TEST(AMS, AMSRestart)
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    params->set("write final state", true);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, AMSRestart2)
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    params->set("write steps", 50);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, AMSRestart3)
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    params->set("write time steps", 1000);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, AMSRestart4)
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    params->set("write steps", 50);
    params->set("maximum iterations", 0);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, TAMSRestart)
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("write final state", true);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, TAMSRestart2)
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("write steps", 50);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, TAMSRestart3)
{
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("write time steps", 1000);
    params->set("write final state", false);

    restart_test(params);
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
