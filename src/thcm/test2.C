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
#include <Teuchos_Version.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#include <sstream>
#include "GlobalDefinitions.H"


using Teuchos::RCP;
using Teuchos::rcp;

RCP<std::ostream> outFile;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
	Epetra_SerialComm Comm;
#endif
	// -------------------------------------------------------------------
	// Setup stream for INFO, DEBUG, ERROR and WARNING Macros
	// -------------------------------------------------------------------
	outFile = rcp(&std::cout, false);

	const int numDomainElements = 16;
	const int numRangeElements = 32;
	const int indexBase = 0;
	Epetra_Map map(numDomainElements, indexBase, Comm);
	Epetra_Map map2(numRangeElements, indexBase, Comm);
	RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(Copy, map, 8));
	RCP<Epetra_CrsMatrix> B = rcp(new Epetra_CrsMatrix(Copy, map, map2, 8));
	int numMyElmnts = map.NumMyElements();
	DEBVAR(numMyElmnts);
	DEBVAR(map);
	int *myGlobalElements = map.MyGlobalElements();
	int rowidx;
	double valuesA[3];
	int   indicesA[3];
	double valuesB[3];
	int   indicesB[3];
	for (int i = 0; i != numMyElmnts; ++i)
	{
		rowidx = myGlobalElements[i];
		valuesB[0]  = 121;
		valuesB[1]  = 32;
		indicesB[0] = rowidx;
		indicesB[1] = rowidx + numMyElmnts;
		B->InsertGlobalValues(rowidx, 2, valuesB, indicesB);

		if (rowidx == 0)
		{
			valuesA[0]  = 6;
			valuesA[1]  = 3;
			indicesA[0] = rowidx;
			indicesA[1] = rowidx + 1;
			A->InsertGlobalValues(rowidx, 2, valuesA, indicesA);
		}
		else if (rowidx == numDomainElements - 1)
		{
			valuesA[0]  = 3;
			valuesA[1]  = 6;
			indicesA[0] = rowidx;
			indicesA[1] = rowidx - 1;
			A->InsertGlobalValues(rowidx, 2, valuesA, indicesA);
		}
		else
		{
			valuesA[0]  = 3;
			valuesA[1]  = 6;
			valuesA[1]  = 3;
			indicesA[0] = rowidx - 1;
			indicesA[1] = rowidx;
			indicesA[2] = rowidx + 1;
			A->InsertGlobalValues(rowidx, 2, valuesA, indicesA);
		}
	}
	A->FillComplete(map, map);
	B->FillComplete(map, map2);
	DEBVAR(*A);
	DEBVAR(*B);
	
	//------------------------------------------------------------------
	// Finalize MPI
	//------------------------------------------------------------------
	MPI_Finalize();
	return 0;
}
