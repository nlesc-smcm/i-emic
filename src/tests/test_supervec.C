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

#include "SuperVector.H"
#include "CoupledModel.H"
#include "Atmosphere.H"
#include "GlobalDefinitions.H"
#include "THCMdefs.H"

#include "gtest/gtest.h" // google test

//------------------------------------------------------------------
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
			updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
			ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
	
			// Create Atmosphere object
			RCP<Teuchos::ParameterList> atmosphereParams =
				rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("atmosphere_params.xml",
										atmosphereParams.ptr());
			atmos = std::make_shared<Atmosphere>(atmosphereParams);
		
			// Create CoupledModel
			RCP<Teuchos::ParameterList> coupledmodelParams =
				rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("coupledmodel_params.xml",
										coupledmodelParams.ptr());
			coupledModel =
				std::make_shared<CoupledModel>(ocean, atmos, coupledmodelParams);
		}

	// destructor
	~IEMIC()
		{}
};

//------------------------------------------------------------------
TEST(SuperVector, ComplexDot)
{
	SuperVector x1 = *coupledModel->getSolution();
	SuperVector x2 = *coupledModel->getSolution();
	SuperVector y1 = *coupledModel->getSolution();
	SuperVector y2 = *coupledModel->getSolution();

	x1.putScalar(1.0);
	x2.putScalar(2.0);
	y1.putScalar(3.0);
	y2.putScalar(4.0);

	ComplexSuperVector z1(x1,y1);
	ComplexSuperVector z2(x2,y2);

	std::complex<double> result = z1.dot(z2);
	EXPECT_EQ(result.real(),43904);
	EXPECT_EQ(result.imag(),-6272);
}

//------------------------------------------------------------------
TEST(SuperVector, Wrapper)
{
	std::vector<double> vec1 = {1,2,3,4,5,6,7,8,9,10};
	std::vector<double> vec2 = {2,2,3,2,5,2,2,2,9,2};
	std::shared_ptr<std::vector<double> > spvec1 =
		std::make_shared<std::vector<double> >(vec1);
	std::shared_ptr<std::vector<double> > spvec2 =
		std::make_shared<std::vector<double> >(vec2);
	
 	SuperVector wrvec1(spvec1);
	SuperVector wrvec2(spvec2);

	EXPECT_EQ(wrvec1.length(), vec1.size());
	EXPECT_EQ(wrvec2.length(), vec2.size());

	EXPECT_NEAR(wrvec1.norm(), 19.621416870348583, 1e-10);
	EXPECT_NEAR(wrvec2.norm(), 11.958260743101398, 1e-10);

	EXPECT_EQ(wrvec1.dot(wrvec2), 191);
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
	return 0;
}
