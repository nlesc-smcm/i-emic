#include <Epetra_config.h>

#ifdef HAVE_MPI
#  include <mpi.h>
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif // HAVE_MPI

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Version.h>

#include <sstream>
#include <stdexcept>

int main (int argc, char *argv[])
{
	using std::cout;
	using std::endl;
	
#ifdef HAVE_MPI
	MPI_Init (&argc, &argv);
	Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
	Epetra_SerialComm comm;
#endif // HAVE_MPI
	
	const int myRank = comm.MyPID ();
	const int numProcs = comm.NumProc ();
	
	if (myRank == 0)
	{
		// Print out the Epetra software version.
		cout << Epetra_Version () << endl << endl
			 << "Total number of processes: " << numProcs << endl;
	}


	// The number of rows and columns in the matrix.
	const int numGlobalElements = 50;

	// Construct a Map that puts approximately the same number of
	// equations on each processor.
	const int indexBase = 0;
	Epetra_Map map (numGlobalElements, indexBase, comm);

	// Get the list of global indices that this process owns.  In this
	// example, this is unnecessary, because we know that we created a
	// contiguous Map (see above).  (Thus, we really only need the min
	// and max global index on this process.)  However, in general, we
	// don't know what global indices the Map owns, so if we plan to add
	// entries into the sparse matrix using global indices, we have to
	// get the list of global indices this process owns.
	const int numMyElements = map.NumMyElements ();

	int* myGlobalElements = NULL;

	myGlobalElements = map.MyGlobalElements ();

	// In general, tests like this really should synchronize across all
	// processes.  However, the likely cause for this case is a
	// misconfiguration of Epetra, so we expect it to happen on all
	// processes, if it happens at all.
	if (numMyElements > 0 && myGlobalElements == NULL) {
		throw std::logic_error ("Failed to get the list of global indices");
	}

	if (myRank == 0) {
		cout << endl << "Creating the sparse matrix" << endl;
	}

	// Create a Epetra sparse matrix whose rows have distribution given
	// by the Map.  The max number of entries per row is 3.
	Epetra_CrsMatrix A (Copy, map, 3);

	// Local error code for use below.
	int lclerr = 0;

	// Fill the sparse matrix, one row at a time.  InsertGlobalValues
	// adds entries to the sparse matrix, using global column indices.
	// It changes both the graph structure and the values.
	double tempVals[3];
	int tempGblInds[3];
	for (int i = 0; i < numMyElements; ++i) {
		// A(0, 0:1) = [2, -1]
		if (myGlobalElements[i] == 0) {
			tempVals[0] = 2.0;
			tempVals[1] = -1.0;
			tempGblInds[0] = myGlobalElements[i];
			tempGblInds[1] = myGlobalElements[i] + 1;
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues (myGlobalElements[i], 2, tempVals, tempGblInds);
			}
			if (lclerr != 0) {
				break;
			}
		}
		// A(N-1, N-2:N-1) = [-1, 2]
		else if (myGlobalElements[i] == numGlobalElements - 1) {
			tempVals[0] = -1.0;
			tempVals[1] = 2.0;
			tempGblInds[0] = myGlobalElements[i] - 1;
			tempGblInds[1] = myGlobalElements[i];
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues (myGlobalElements[i], 2, tempVals, tempGblInds);
			}
			if (lclerr != 0) {
				break;
			}
		}
		// A(i, i-1:i+1) = [-1, 2, -1]
		else {
			tempVals[0] = -1.0;
			tempVals[1] = 2.0;
			tempVals[2] = -1.0;
			tempGblInds[0] = myGlobalElements[i] - 1;
			tempGblInds[1] = myGlobalElements[i];
			tempGblInds[2] = myGlobalElements[i] + 1;
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues (myGlobalElements[i], 3, tempVals, tempGblInds);
			}
			if (lclerr != 0) {
				break;
			}
		}
	}

	// If any process failed to insert at least one entry, throw.
	int gblerr = 0;
	(void) comm.MaxAll (&lclerr, &gblerr, 1);
	if (gblerr != 0) {
		throw std::runtime_error ("Some process failed to insert an entry.");
	}

	// Tell the sparse matrix that we are done adding entries to it.
	gblerr = A.FillComplete ();
	if (gblerr != 0) {
		std::ostringstream os;
		os << "A.FillComplete() failed with error code " << gblerr << ".";
		throw std::runtime_error (os.str ());
	}

	// Number of iterations
	const int niters = 500;
	// Desired (absolute) residual tolerance
	const double tolerance = 1.0e-2;


	//
	// Now we're going to change values in the sparse matrix and run the
	// power method again.
	//

	//
	// Increase diagonal dominance
	//
	if (myRank == 0) {
		cout << endl << "Increasing magnitude of A(0,0), solving again" << endl;
	}

	if (A.RowMap ().MyGID (0)) {
		// Get a copy of the row with with global index 0.  Modify the
		// diagonal entry of that row.  Submit the modified values to the
		// matrix.
		const int gidOfFirstRow = 0;
		// Since 0 is a GID in the row Map on the calling process,
		// lidOfFirstRow is a valid LID of that GID in the row Map.
		const int lidOfFirstRow = A.RowMap ().LID (gidOfFirstRow);
		int numEntriesInRow = A.NumMyEntries (lidOfFirstRow);
		double* rowvals = new double [numEntriesInRow];
		int* rowinds = new int [numEntriesInRow];

		// Get a copy of the entries and column indices of global row 0.
		// Get global column indices, so that we can figure out which
		// entry corresponds to the diagonal entry in this row.  (The row
		// Map and column Map of the matrix may differ, so their local
		// indices may not be the same.)
		//
		// Note that it's legal (though we don't exercise it in this
		// example) for the row Map of the sparse matrix not to be one to
		// one.  This means that more than one process might own entries
		// in the first row.  In general, multiple processes might own the
		// (0,0) entry, so that the global A(0,0) value is really the sum
		// of all processes' values for that entry.  However, scaling the
		// entry by a constant factor distributes across that sum, so it's
		// OK to do so.
		if (lclerr == 0) {
			lclerr = A.ExtractGlobalRowCopy (gidOfFirstRow,
											 numEntriesInRow, numEntriesInRow,
											 rowvals, rowinds);
		}
		if (lclerr == 0) { // no error
			for (int i = 0; i < numEntriesInRow; ++i) {
				if (rowinds[i] == gidOfFirstRow) {
					// We have found the diagonal entry; modify it.
					rowvals[i] *= 10.0;
				}
			}
			// "Replace global values" means modify the values, but not the
			// structure of the sparse matrix.  If the specified columns
			// aren't already populated in this row on this process, then this
			// method fails and returns nonzero.  Since we have already called
			// FillComplete() on this matrix, we may not modify its graph
			// structure any more.
			if (lclerr == 0) {
				lclerr = A.ReplaceGlobalValues (gidOfFirstRow, numEntriesInRow,
												rowvals, rowinds);
			}
		}

		if (rowvals != NULL) {
			delete [] rowvals;
		}
		if (rowinds != NULL) {
			delete [] rowinds;
		}
	}

	// If the owning process(es) of global row 0 failed to replace the
	// one entry, throw.
	gblerr = 0;
	(void) comm.MaxAll (&lclerr, &gblerr, 1);
	if (gblerr != 0) {
		throw std::runtime_error ("One of the owning process(es) of global "
								  "row 0 failed to replace an entry.");
	}
	

#ifdef HAVE_MPI
	(void) MPI_Finalize ();
#endif // HAVE_MPI

	return 0;
}
