#include <Epetra_config.h>

#  ifdef HAVE_MPI
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif // HAVE_MPI

#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Teuchos_RCP.hpp>

#include "gtest/gtest.h" // google test

namespace
{
    Teuchos::RCP<Epetra_Comm>  comm;
    Teuchos::RCP<Epetra_Map>    map;
    Teuchos::RCP<Epetra_Map> submap;
    Teuchos::RCP<Epetra_Import> imp;
}

//------------------------------------------------------------------
TEST(Map, Construction)
{
    // Simple map construction. 
    int N    = 100;
    int pid  = comm->MyPID();
    int size = comm->NumProc();
    int myN  = N / size;

    int lastpid = size - 1;
    if (pid == lastpid)
    {
        myN += N % size;
    }

    int sum;
    comm->SumAll(&myN, &sum, 1);
    EXPECT_EQ(sum, N);

    int *myGlobalElements = new int[myN];
    for (int i = 0; i != myN; ++i)
        myGlobalElements[i] = pid * (N / size) + i;

    map = Teuchos::rcp(new Epetra_Map(N, myN, myGlobalElements, 0, *comm));

    int k = 0;
    for (int i = 0; i < myN; i+=4)
        myGlobalElements[k++] = N - 1;

    submap = Teuchos::rcp(new Epetra_Map(-1, k, myGlobalElements, 0, *comm));

    imp = Teuchos::rcp(new Epetra_Import(*submap, *map));
    
    delete [] myGlobalElements;
}


//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
    comm = Teuchos::rcp(new Epetra_SerialComm() );
#endif

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
