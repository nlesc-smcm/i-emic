#include "TransientFactory.H"

#include "TestDefinitions.H"

#include <stdio.h>

#include "Trilinos_version.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

//------------------------------------------------------------------
namespace // local unnamed namespace (similar to static in C)
{
Teuchos::RCP<Epetra_Comm> comm;
Teuchos::RCP<Epetra_Map> map;
Teuchos::RCP<std::ostringstream> out_stream;
}

class TestModel
{
public:
    using ConstVectorPtr = Teuchos::RCP<const Epetra_Vector>;
protected:
    Teuchos::RCP<Epetra_Map> map_;
    Teuchos::RCP<Epetra_Vector> rhs_;
    Teuchos::RCP<Epetra_Vector> sol_;
    Teuchos::RCP<Epetra_Vector> state_;
    Teuchos::RCP<Epetra_Vector> diagB_;

    Teuchos::RCP<Epetra_CrsMatrix> frc_;
    Teuchos::RCP<Epetra_CrsMatrix> jac_;
public:
    using Vector = Epetra_Vector;
    using VectorPtr = Teuchos::RCP<Vector>;

    TestModel(Teuchos::RCP<Epetra_Map> map)
        :
        map_(map)
        {
            rhs_ = Teuchos::rcp(new Epetra_Vector(*map_));
            sol_ = Teuchos::rcp(new Epetra_Vector(*map_));
            state_ = Teuchos::rcp(new Epetra_Vector(*map_));

            std::vector<double> values(2, 1);
            diagB_ = Teuchos::rcp(new Epetra_Vector(Copy, *map_, &values[0]));
        }

    Teuchos::RCP<Epetra_Comm> Comm() const
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
                CHECK_ZERO(frc_->InsertGlobalValues(idx, 1, &val, &idx));
            CHECK_ZERO(frc_->FillComplete());
        }

    void computeMassMat() {}

    Teuchos::RCP<Epetra_CrsMatrix> getForcing() {return frc_;}

    Teuchos::RCP<Epetra_Vector> getVector(char mode, Teuchos::RCP<Epetra_Vector> vec)
        {
            if (mode == 'C') // copy
            {
                Teuchos::RCP<Epetra_Vector> copy = Teuchos::rcp(new Epetra_Vector(*vec));
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

    Teuchos::RCP<Epetra_Vector> getSolution(char mode = 'C')
        {
            return getVector(mode, sol_);
        }

    Teuchos::RCP<Epetra_CrsMatrix> getJacobian()
        {
            return jac_;
        }

    void solve(Teuchos::RCP<const Epetra_Vector> rhs = Teuchos::null)
        {
            *sol_ = *rhs;
        }

    void applyMatrix(Epetra_MultiVector const &v, Epetra_MultiVector &out) {}
    void applyMassMat(Epetra_MultiVector const &v, Epetra_MultiVector &out) { out = v; }
    void preProcess() {}
};

Teuchos::RCP<Transient<Teuchos::RCP<const Epetra_Vector> > >
createDoubleWell(
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Epetra_MultiVector> V = Teuchos::null)
{
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

    if (V != Teuchos::null)
        return TransientFactory(model, params, sol1, sol2, sol3, V);
    else
        return TransientFactory(model, params, sol1, sol2, sol3);
}

template<typename T>
void set_parameter(Teuchos::RCP<Teuchos::ParameterList> &params,
                   std::string const &name, T value)
{
    if (!params->isParameter(name))
        params->set(name, value);
}

void set_default_parameters(Teuchos::RCP<Teuchos::ParameterList> &params)
{
    set_parameter(params, "theta", 0.0);
    set_parameter(params, "sigma", 1.0);
    set_parameter(params, "dof", 1);
    set_parameter(params, "var", 0);
    set_parameter(params, "noise seed", 5);
    set_parameter(params, "ams seed", 2);

    set_parameter(params, "time step", 0.01);
    set_parameter(params, "maximum time", 2.0);
    set_parameter(params, "B distance", 0.05);
    set_parameter(params, "number of experiments", 20);
    set_parameter(params, "maximum iterations", 10);
}

void restart_test(Teuchos::RCP<Teuchos::ParameterList> params)
{
    const int default_maxit = 10;
    set_default_parameters(params);

    set_parameter(params, "write file", "out_data.h5");
    params->set("noise seed", 42);

    int maxit = params->get("maximum iterations", -1);
    int write_time_steps = params->get("write time steps", -1);

    remove("out_data.h5");

    // testing::internal::CaptureStdout();
    out_stream->str("");

    // Test that the method can start
    auto ams = createDoubleWell(params);
    ams->run();

    // std::string output = testing::internal::GetCapturedStdout();
    std::string output = out_stream->str();

    // std::cout << output;

    EXPECT_NE(output.find("Initialization"), std::string::npos);
    if (maxit > 0)
    {
        EXPECT_NE(output.find("AMS: " + std::to_string(maxit)),
                  std::string::npos);
    }
    EXPECT_EQ(output.find("AMS: " + std::to_string(maxit + 1)),
              std::string::npos);

    // Test that the file is read and maximum iterations are used correctly
    params->set("read file", "out_data.h5");
    params->set("write file", "");
    params->set("maximum iterations", default_maxit);

    out_stream->str("");

    auto ams2 = createDoubleWell(params);
    ams2->run();

    std::string output2 = out_stream->str();

    // std::cout << output2;

    EXPECT_EQ(output2.find("Initialization"), std::string::npos);
    EXPECT_EQ(output2.find("AMS: " + std::to_string(default_maxit + 1)),
              std::string::npos);

    // Test that the file is read and the method actually converges
    params->set("read file", "out_data.h5");
    params->set("write file", "");
    params->set("maximum iterations", 10000);

    out_stream->str("");

    ams2 = createDoubleWell(params);
    ams2->run();

    output2 = out_stream->str();

    // std::cout << output2;

    EXPECT_EQ(output2.find("Initialization"), std::string::npos);
    if (write_time_steps < 0)
    {
        EXPECT_EQ(output2.find("AMS: " + std::to_string(maxit)),
                  std::string::npos);
    }
    else if (maxit < default_maxit)
    {
        EXPECT_NE(output2.find("AMS: " + std::to_string(default_maxit)),
                  std::string::npos);
    }
    EXPECT_NE(output2.find("AMS: " + std::to_string(maxit + 1)),
              std::string::npos);
    EXPECT_EQ(output2.find("AMS: " + std::to_string(default_maxit * 8)),
              std::string::npos);

    if (ams2->get_mfpt() > 0)
    {
        EXPECT_GT(ams2->get_mfpt(), 6);
        EXPECT_LT(ams2->get_mfpt(), 20);
    }
    else
        EXPECT_NEAR(ams2->get_probability(), 0.157, 5e-1);
}

#if TRILINOS_MAJOR_MINOR_VERSION > 121300

//------------------------------------------------------------------
TEST(AMS, AMSRestart)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "AMS");
    params->set("write final state", true);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, AMSRestart2)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "AMS");
    params->set("write steps", 10);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, AMSRestart3)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "AMS");
    params->set("write time steps", 500);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, AMSRestart4)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "AMS");
    params->set("write steps", 10);
    params->set("maximum iterations", 0);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, AMSRestart5)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "AMS");
    params->set("write final state", true);
    params->set("maximum iterations", 30);

    restart_test(params);
}

#endif //TRILINOS_MAJOR_MINOR_VERSION

//------------------------------------------------------------------
TEST(AMS, AMSConvergence)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "AMS");
    params->set("maximum iterations", 1000);
    params->set("number of experiments", 200);
    set_default_parameters(params);

    auto ams = createDoubleWell(params);
    ams->run();

    EXPECT_NEAR(ams->get_mfpt(), 7, 1);
}

//------------------------------------------------------------------
TEST(AMS, ProjectedAMSConvergence)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "AMS");
    params->set("maximum iterations", 10000);
    params->set("number of experiments", 200);
    set_default_parameters(params);

    std::vector<double> values(2);
    values[0] = 1;
    values[1] = 0;

    Teuchos::RCP<Epetra_Vector> V = Teuchos::rcp(new Epetra_Vector(Copy, *map, &values[0]));

    auto ams = createDoubleWell(params, V);
    ams->run();

    EXPECT_NEAR(ams->get_mfpt(), 6.3, 1);
}

#if TRILINOS_MAJOR_MINOR_VERSION > 121300

//------------------------------------------------------------------
TEST(AMS, TAMSRestart)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("write final state", true);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, TAMSRestart2)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("write steps", 10);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, TAMSRestart3)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("write time steps", 1000);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, TAMSRestart4)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("write steps", 10);
    params->set("maximum iterations", 0);
    params->set("write final state", false);

    restart_test(params);
}

//------------------------------------------------------------------
TEST(AMS, TAMSRestart5)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("write final state", true);
    params->set("maximum iterations", 30);

    restart_test(params);
}

#endif //TRILINOS_MAJOR_MINOR_VERSION

//------------------------------------------------------------------
TEST(AMS, TAMSConvergence)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("maximum iterations", 10000);
    params->set("number of experiments", 200);
    set_default_parameters(params);

    auto tams = createDoubleWell(params);
    tams->run();

    EXPECT_NEAR(tams->get_probability(), 0.157, 1e-2);
}

//------------------------------------------------------------------
TEST(AMS, ProjectedTAMSConvergence)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "TAMS");
    params->set("maximum iterations", 10000);
    params->set("number of experiments", 200);
    set_default_parameters(params);

    std::vector<double> values(2);
    values[0] = 1;
    values[1] = 0;

    Teuchos::RCP<Epetra_Vector> V = Teuchos::rcp(new Epetra_Vector(Copy, *map, &values[0]));

    auto tams = createDoubleWell(params, V);
    tams->run();

    EXPECT_NEAR(tams->get_probability(), 0.215, 1e-2);
}

//------------------------------------------------------------------
TEST(AMS, MCConvergence)
{
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("method", "Naive");
    params->set("maximum iterations", 10000);
    params->set("number of experiments", 200);
    set_default_parameters(params);

    auto mc = createDoubleWell(params);
    mc->run();

    EXPECT_NEAR(mc->get_probability(), 0.157, 1e-2);
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    out_stream = Teuchos::rcp(new std::ostringstream);
    comm = initializeEnvironment(argc, argv, out_stream);
    map = Teuchos::rcp(new Epetra_Map(2, 0, *comm));

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
