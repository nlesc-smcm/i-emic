#include "TestDefinitions.H"

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"

namespace
{
    Teuchos::RCP<Epetra_Comm>  comm;
    Teuchos::RCP<Epetra_Map>    map;
    Teuchos::RCP<Epetra_Map> submap;
    Teuchos::RCP<Epetra_Import> imp;

    int N;
}

//------------------------------------------------------------------
TEST(Map, Construction)
{
    bool failed = false;
    try
    {
        // domain size
        N = 100;

        // Simple domain decomposition
        int pid  = comm->MyPID();
        int size = comm->NumProc();
        int myN  = N / size;

        int lastpid = size - 1;
        if (pid == lastpid)
        {
            myN += N % size;
        }

        // check domain decomposition
        int sum;
        comm->SumAll(&myN, &sum, 1);
        EXPECT_EQ(sum, N);

        // create map
        int *myGlobalElements = new int[myN];
        for (int i = 0; i != myN; ++i)
            myGlobalElements[i] = pid * (N / size) + i;

        map = Teuchos::rcp(new Epetra_Map(N, myN, myGlobalElements, 0, *comm));

        // create submap
        int k = 0;
        myGlobalElements[k++] = N - 1;

        submap = Teuchos::rcp(new Epetra_Map(-1, k, myGlobalElements, 0, *comm));

        // create importer
        imp = Teuchos::rcp(new Epetra_Import(*submap, *map));

        delete [] myGlobalElements;
    }
    catch (...)
    {
        failed = true;
    }
    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(Map, Usage)
{
    // create and randomize vector
    Epetra_Vector x(*map);
    Epetra_Vector y(*submap);
    x.Random();

    // set last element
    int lastgid = N-1;
    int lid = map->LID(lastgid);
    double val = 3.141592;
    if (lid >= 0)
        x[lid] = val;

    y.Import(x, *imp, Insert);

    for (int i = 0; i != submap->NumMyElements(); ++i)
    {
        EXPECT_EQ(y[i], val);
    }
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
