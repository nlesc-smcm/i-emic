// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>

#  ifdef HAVE_MPI
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif // HAVE_MPI

#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Epetra_Vector.h>

#include "GlobalDefinitions.H"
#include "Ocean.H"
#include "Newton.H"
#include "Singleton.H"

using Teuchos::RCP;
using Teuchos::rcp;

RCP<std::ostream> outFile;

int main(int argc, char **argv)
{
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	RCP<Epetra_MpiComm> Comm =
		rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	RCP<Epetra_SerialComm> Comm =
		rcp(new Epetra_SerialComm());
#endif
	// Specify output streams
	outFile = rcp(&std::cout, false);
	// Initialize ocean model (THCM)
	RCP<OceanTheta> ocean = rcp(new OceanTheta(Comm));
	Newton<RCP<OceanTheta>, RCP<Epetra_Vector> > newton(ocean);

	ocean->parkModel();
	ocean->computeJacobian();
	ocean->computeRHS();
	ocean->solve();
	
    //--------------------------------------------------------
	// Finalize MPI
	//--------------------------------------------------------
	MPI_Finalize();	
}
