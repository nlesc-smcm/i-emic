#include "TestDefinitions.H"

//------------------------------------------------------------------
double norm(Teuchos::RCP<Epetra_Vector> vec)
{
    double nrm;
    vec->Norm2(&nrm);
    return nrm;
}

//------------------------------------------------------------------
double norm(std::shared_ptr<std::vector<double> > vec)
{
    int dim = (int) vec->size();
    int incX = 1;
    int incY = 1;
    double dot;
    dot = ddot_(&dim, &(*vec)[0], &incX, &(*vec)[0], &incY);
    return sqrt(dot);
}

//------------------------------------------------------------------
std::shared_ptr<std::vector<double> >
getGatheredVector(Teuchos::RCP<Epetra_Vector> vec)
{
    // get size
    int dim = (int) vec->GlobalLength();

    // Create gather scheme: map and importer
    Teuchos::RCP<Epetra_BlockMap> gmap =
        Utils::AllGather(vec->Map());
    Teuchos::RCP<Epetra_Import> imp =
        Teuchos::rcp(new Epetra_Import(*gmap, vec->Map()));

    // Create vector to hold the gathered state
    Teuchos::RCP<Epetra_Vector> gathered =
        Teuchos::rcp(new Epetra_Vector(*gmap));

    // Gather the state into <gathered>
    gathered->Import(*vec, *imp, Insert);

    // Create full state for serial atmosphere
    std::shared_ptr< std::vector<double> > gvec =
        std::make_shared<std::vector<double> >(dim, 0.0);

    // Extract gvec from gathered
    gathered->ExtractCopy(&(*gvec)[0], dim);

    return gvec;
}
