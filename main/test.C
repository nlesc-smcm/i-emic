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
// This is will function as a fixture for the test routines
class IEMIC : public testing::Test
{

public:
	// constructor
	IEMIC()
		{
			// Create parallel Ocean 
			RCP<Teuchos::ParameterList> oceanParams =
				rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("ocean_params.xml", oceanParams.ptr());
			ocean_ = Teuchos::rcp(new Ocean(comm, oceanParams));
	
			// Create Atmosphere object
			RCP<Teuchos::ParameterList> atmosphereParams =
				rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("atmosphere_params.xml",
										atmosphereParams.ptr());
			atmos_ = std::make_shared<Atmosphere>(atmosphereParams);
		
			// Create CoupledModel
			RCP<Teuchos::ParameterList> coupledmodelParams =
				rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("coupledmodel_params.xml",
										coupledmodelParams.ptr());
			coupledModel_ =
				std::make_shared<CoupledModel>(ocean_, atmos_, coupledmodelParams);

		}
	
	// destructor
	~IEMIC()
		{}
	
	RCP<Ocean>                    ocean_;
	std::shared_ptr<Atmosphere>   atmos_;
	std::shared_ptr<CoupledModel> coupledModel_;
};

//------------------------------------------------------------------
TEST_F(IEMIC, AtmosphereView)
{
	SuperVector stateView = *atmos_->getState('V');
	SuperVector stateCopy = *atmos_->getState('C');

	double normView = stateView.norm();
	double normCopy = stateCopy.norm();

	EXPECT_EQ(normView, normCopy);
}

//------------------------------------------------------------------
TEST_F(IEMIC, AtmosphereCopy)
{
	SuperVector stateView = *atmos_->getState('V');
	SuperVector stateCopy = *atmos_->getState('C');

	double normView = stateView.norm();
	double normCopy = stateCopy.norm();

	EXPECT_EQ(normView, normCopy);
}

//------------------------------------------------------------------
TEST_F(IEMIC, CouplingBlock)
{
	std::vector<double> C12values;
	std::vector<int>    C12rows;
	std::vector<double> C21values;
	std::vector<int>    C21rows;

	ocean_->getAtmosBlock(C12values, C12rows);
	atmos_->getOceanBlock(C21values, C21rows);

	SuperVector vec0(ocean_->getState('C')->getOceanVector(),
					 atmos_->getState('C')->getAtmosVector());

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
int main(int argc, char **argv)
{
	// Initialize the environment:
	initializeEnvironment(argc, argv);

	if (outFile == Teuchos::null)
		throw std::runtime_error("ERROR: Specify output streams");

	// -------------------------------------------------------
	// TEST_FING 
	::testing::InitGoogleTest(&argc, argv);
	int out = RUN_ALL_TESTS();
	// -------------------------------------------------------

	MPI_Finalize();
	return out;
}

//------------------------------------------------------------------
void testCouplingBlock()
{
	

}

