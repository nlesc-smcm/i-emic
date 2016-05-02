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
#include "Continuation.H"
#include "Atmosphere.H"
#include "GlobalDefinitions.H"
#include "THCMdefs.H"
#include "JDQZInterface.H"
#include "MyParameterList.H"

#include "jdqz.H"

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
	Continuation<std::shared_ptr<CoupledModel>,
				 RCP<Teuchos::ParameterList> >	continuation;
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
			
			// Create Continuation
			RCP<Teuchos::ParameterList> continuationParams = rcp(new Teuchos::ParameterList);
			updateParametersFromXmlFile("continuation_params.xml", continuationParams.ptr());
			
			continuation = 
				Continuation<std::shared_ptr<CoupledModel>, RCP<Teuchos::ParameterList> >
				(coupledModel, continuationParams);
		}

	// destructor
	~IEMIC()
		{}
};

TEST(JDQZ, General)
{
	continuation.run();

	coupledModel->printJacobian("testjac");
	
	SuperVector x = *coupledModel->getSolution();
	SuperVector y = *coupledModel->getSolution();
	x.zero(); y.zero();
	ComplexSuperVector z(x, y);
	ComplexSuperVector residue(x,y);
	ComplexSuperVector tmp(x,y);

	JDQZInterface<CoupledModel> matrix(*coupledModel);	

	JDQZ<JDQZInterface<CoupledModel> > jdqz(matrix, z);
	
	std::map<std::string, double> list;	
	list["Shift (real part)"]         = 0.0;
	list["Number of eigenvalues"]     = 3;
	list["Max size search space"]     = 20;
	list["Min size search space"]     = 10;
	list["Max JD iterations"]         = 200;
	list["Tracking parameter"]        = 1e-8;
	list["Criterion for Ritz values"] = 0;
	list["Linear solver"]             = 1;
	list["GMRES search space"]        = 20;		
	MyParameterList params(list);
	
	jdqz.setParameters(params);
	jdqz.printParameters();
	jdqz.solve();

	std::vector<ComplexSuperVector>    eivec = jdqz.getEigenVectors();
	std::vector<std::complex<double> > alpha = jdqz.getAlpha();
	std::vector<std::complex<double> > beta  = jdqz.getBeta();


	for (int j = 0; j != jdqz.kmax(); ++j)
	{
		residue.zero();
		tmp.zero();
		matrix.AMUL(eivec[j], residue);
		residue.scale(beta[j]);
		matrix.BMUL(eivec[j], tmp);
		residue.axpy(-alpha[j], tmp);
		std::cout << "alpha: " <<  std::setw(20) << alpha[j]
				  << " beta: " << std::setw(20) << beta[j]
				  << " " << residue.norm() << std::endl;
		
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
