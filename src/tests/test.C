#include <Epetra_config.h>

#  ifdef HAVE_MPI
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif // HAVE_MPI

#include <EpetraExt_HDF5.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <ios>
#include <iomanip>
#include <memory>
#include <vector>
#include <array>
#include <stack>
#include <string>

#include <cassert>

// Everything I-EMIC
#include "SuperVector.H"
#include "CoupledModel.H"
#include "CouplingBlock.H"
#include "Atmosphere.H"
#include "ThetaStepper.H"
#include "Continuation.H"
#include "GlobalDefinitions.H"
#include "THCMdefs.H"
#include "IDRSolver.H"
#include "GMRESSolver.H"
#include "NumericalJacobian.H"

#include "gtest/gtest.h" // google test

using Teuchos::RCP;
using Teuchos::rcp;

//------------------------------------------------------------------
// A few globals (see GlobalDefinitions.H)
//------------------------------------------------------------------
RCP<std::ostream> outFile;      // output file
ProfileType       profile;      // profile
std::stack<Timer> timerStack;   // timing stack
RCP<Epetra_Comm>  comm;         // communicator object

//------------------------------------------------------------------
RCP<std::ostream> outputFiles()
{
	Teuchos::RCP<std::ostream> outFile;
	if (comm->MyPID() < 1)
	{
		std::ostringstream infofile;     // setting up a filename

		infofile    << "info_"    << comm->MyPID()   << ".txt";

		std::cout << "info for CPU" << comm->MyPID() << " is written to "
				  << infofile.str().c_str() << std::endl;

		outFile = Teuchos::rcp(new std::ofstream(infofile.str().c_str()));
	}
	else
	{
		outFile = Teuchos::rcp(new Teuchos::oblackholestream());
	}
	return outFile;
}

//------------------------------------------------------------------
void initializeEnvironment(int argc, char **argv)
{
#ifdef HAVE_MPI           // Initialize communicator
	MPI_Init(&argc, &argv);
	comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	comm = rcp(new Epetra_SerialComm());
#endif
	outFile = outputFiles(); 	// Initialize output files
}

//------------------------------------------------------------------
namespace // local (unnamed) namespace (similar to static in C)
{	
	RCP<Ocean>                    ocean;
	std::shared_ptr<Atmosphere>   atmos;
	std::shared_ptr<CoupledModel> coupledModel;
}

//------------------------------------------------------------------
class IEMIC : public testing::Environment
{
public:
	// constructor
	IEMIC()
		{
			// Create parallel Ocean 
			RCP<Teuchos::ParameterList> oceanParams =
				rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("ocean_test.xml", oceanParams.ptr());
			ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
	
			// Create Atmosphere object
			RCP<Teuchos::ParameterList> atmosphereParams =
				rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("atmos_test.xml",
										atmosphereParams.ptr());
			atmos = std::make_shared<Atmosphere>(atmosphereParams);
		
			// Create CoupledModel
			RCP<Teuchos::ParameterList> coupledmodelParams =
				rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("coupledmodel_test.xml",
										coupledmodelParams.ptr());
			coupledModel =
				std::make_shared<CoupledModel>(ocean, atmos, coupledmodelParams);
		}

	// destructor
	~IEMIC()
		{}
};

//------------------------------------------------------------------
TEST(IEMIC, AtmosphereView)
{
	SuperVector stateView = *atmos->getState('V');
	SuperVector stateCopy = *atmos->getState('C');

	double normView = stateView.norm();
	double normCopy = stateCopy.norm();

	EXPECT_EQ(normView, normCopy);
}

//------------------------------------------------------------------
TEST(IEMIC, AtmosphereCopy)
{
	SuperVector stateView = *atmos->getState('V');
	SuperVector stateCopy = *atmos->getState('C');

	double normView = stateView.norm();
	double normCopy = stateCopy.norm();

	EXPECT_EQ(normView, normCopy);
}

//------------------------------------------------------------------
TEST(IEMIC, CouplingBlock)
{
	std::vector<double> C12values;
	std::vector<int>    C12rows;
	std::vector<double> C21values;
	std::vector<int>    C21rows;

	ocean->getAtmosBlock(C12values, C12rows);
	atmos->getOceanBlock(C21values, C21rows);

	SuperVector vec0(ocean->getState('C')->getOceanVector(),
					 atmos->getState('C')->getAtmosVector());

	SuperVector vec1(vec0);
	SuperVector vec2(vec0);

	vec0.putScalar(0.0);
	vec1.putScalar(3.14);	

	CouplingBlock C12("AO", C12values, C12rows, C21rows);
	CouplingBlock C21("OA", C21values, C21rows, C12rows);

	C12.applyMatrix(vec1, vec0);
	EXPECT_EQ(3.14, vec0.norm());
}

//------------------------------------------------------------------
TEST(IEMIC, OceanScaling)
{
	Teuchos::RCP<Epetra_Vector> rowScaling = ocean->getRowScaling();
	Teuchos::RCP<Epetra_Vector> colScaling = ocean->getColScaling();
	std::ofstream rowScalFile, colScalFile;
	rowScalFile.open("rowscaling.txt");	
	colScalFile.open("colscaling.txt");
	Teuchos::RCP<Epetra_MultiVector> rowScalingGath = Utils::Gather(*rowScaling, 0);
	Teuchos::RCP<Epetra_MultiVector> colScalingGath = Utils::Gather(*colScaling, 0);
	std::vector<double> rowScalingVec(rowScaling->GlobalLength(), 0.0);
	std::vector<double> colScalingVec(colScaling->GlobalLength(), 0.0);
	if (rowScaling->Map().Comm().MyPID() == 0)
	{
		(*rowScalingGath)(0)->ExtractCopy(&rowScalingVec[0]);
		(*colScalingGath)(0)->ExtractCopy(&colScalingVec[0]);

		INFO("    ++++ printing rowscaling ++++");
		for (auto &r : rowScalingVec)
			rowScalFile << std::setprecision(12) << r << '\n';

		INFO("    ++++ printing colscaling ++++");
		for (auto &c : rowScalingVec)
			colScalFile << std::setprecision(12) << c << '\n';
	}

	EXPECT_EQ(rowScaling->GlobalLength(), colScaling->GlobalLength());
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
	// Initialize the environment:
	initializeEnvironment(argc, argv);
	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");

	::testing::InitGoogleTest(&argc, argv);
	::testing::AddGlobalTestEnvironment(new IEMIC);
	// -------------------------------------------------------
	// TESTING 
	int out = RUN_ALL_TESTS();
	// -------------------------------------------------------
	
	// Get rid of possibly parallel objects:
	ocean        = Teuchos::null;
	atmos        = nullptr;
	coupledModel = nullptr;
	
	comm->Barrier();
	std::cout << "TEST exit code proc #" << comm->MyPID()
			  << " " << out << std::endl;
	MPI_Finalize();
	return out;
}
