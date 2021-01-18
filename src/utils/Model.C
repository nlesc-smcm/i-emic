#include "Model.H"

#include "TRIOS_Domain.H"

#include "EpetraExt_Exception.h"
#include "EpetraExt_HDF5.h"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"

// Implementations
//=============================================================================
int Model::loadStateFromFile(std::string const &filename)
{
    INFO("_________________________________________________________");
    if (loadState_)
    {
        INFO("Loading state and parameters from " << filename);
    }
    else
    {
        INFO("Performing only model specific import operations from " << filename);
    }

    // Check whether file exists
    std::ifstream file(filename);
    if (!file)
    {
        WARNING("Can't open " << filename
                << ", continue with trivial state", __FILE__, __LINE__);

        // create trivial state
        state_->PutScalar(0.0);
        return 1;
    }
    else file.close();

    // Create HDF5 object
    EpetraExt::HDF5 HDF5(*comm_);
    Epetra_MultiVector *readState;

    // Open file
    HDF5.Open(filename);

    if (loadState_)
    {
        // Check contents
        if (!HDF5.IsContained("State"))
        {
            ERROR("The group <State> is not contained in hdf5 " << filename,
                  __FILE__, __LINE__);
        }

        // Read the state. To be able to restart with different
        // numbers of procs we do not include the Map in the hdf5. The
        // state as read here will have a linear map and we import it
        // into the current domain decomposition.
        HDF5.Read("State", readState);

        if ( readState->GlobalLength() != getDomain()->GetSolveMap()->NumGlobalElements() )
        {
            WARNING("Loading state from differ #procs", __FILE__, __LINE__);
        }

        // Create importer
        // target map: domain StandardMap
        // source map: state with linear map as read by HDF5.Read
        Teuchos::RCP<Epetra_Import> lin2solve =
            Teuchos::rcp(new Epetra_Import(*(getDomain()->GetSolveMap()),
                                           readState->Map() ));

        // Import state from HDF5 into state_ datamember
        CHECK_ZERO(state_->Import(*((*readState)(0)), *lin2solve, Insert));

        delete readState;

        INFO(" state: ||x|| = " << Utils::norm(state_));

        // Interface between HDF5 and the parameters,
        // put all the <npar> parameters back in atmos.
        std::string parName;
        double parValue;

        // Check contents
        if (!HDF5.IsContained("Parameters"))
        {
            ERROR("The group <Parameters> is not contained in hdf5 " << filename,
                  __FILE__, __LINE__);
        }

        for (int par = 0; par < npar(); ++par)
        {
            parName  = int2par(par);

            // Read continuation parameter and set them in model
            try
            {
                HDF5.Read("Parameters", parName.c_str(), parValue);
            }
            catch (EpetraExt::Exception &e)
            {
                e.Print();
                continue;
            }

            setPar(parName, parValue);
            INFO("   " << parName << " = " << parValue);
        }
    }

    additionalImports(HDF5, filename);

    INFO("_________________________________________________________");
    return 0;
}

//=============================================================================
int Model::saveStateToFile(std::string const &filename)
{
    INFO("_________________________________________________________");
    INFO("Create backup of " << outputFile_);
    copyState(".bak");

    INFO("Writing state and parameters to " << filename);
    INFO("   state: ||x|| = " << Utils::norm(state_));

    // Write state, map and continuation parameter
    EpetraExt::HDF5 HDF5(*comm_);
    HDF5.Create(filename);
    HDF5.Write("State", *state_);

    // Interface between HDF5 and the parameters,
    // store all the <npar> parameters in an HDF5 file.
    std::string parName;
    double parValue;
    for (int par = 0; par < npar(); ++par)
    {
        parName  = int2par(par);
        parValue = getPar(parName);
        INFO("   " << parName << " = " << parValue);
        HDF5.Write("Parameters", parName.c_str(), parValue);
    }

    // Write grid information available in domain object
    HDF5.Write("Grid", "n",   getDomain()->GlobalN());
    HDF5.Write("Grid", "m",   getDomain()->GlobalM());
    HDF5.Write("Grid", "l",   getDomain()->GlobalL());
    HDF5.Write("Grid", "nun", getDomain()->Dof());
    HDF5.Write("Grid", "aux", getDomain()->Aux());

    HDF5.Write("Grid", "xmin", getDomain()->Xmin());
    HDF5.Write("Grid", "xmax", getDomain()->Xmax());
    HDF5.Write("Grid", "ymin", getDomain()->Ymin());
    HDF5.Write("Grid", "ymax", getDomain()->Ymax());
    HDF5.Write("Grid", "hdim", getDomain()->Hdim());


    // Write grid arrays
    const TRIOS::Grid& grid = getDomain()->GetGlobalGrid();
    HDF5.Write("Grid", "x", H5T_NATIVE_DOUBLE, grid.x_.size(), grid.x_.get());
    HDF5.Write("Grid", "y", H5T_NATIVE_DOUBLE, grid.y_.size(), grid.y_.get());
    HDF5.Write("Grid", "z", H5T_NATIVE_DOUBLE, grid.z_.size(), grid.z_.get());
    HDF5.Write("Grid", "xu", H5T_NATIVE_DOUBLE, grid.xu_.size(), grid.xu_.get());
    HDF5.Write("Grid", "yv", H5T_NATIVE_DOUBLE, grid.yv_.size(), grid.yv_.get());
    HDF5.Write("Grid", "zw", H5T_NATIVE_DOUBLE, grid.zw_.size(), grid.zw_.get());

    additionalExports(HDF5, filename);
    comm_->Barrier();

    INFO("_________________________________________________________");
    return 0;
}

//=============================================================================
int Model::copyState(std::string const &append)
{
    if (comm_->MyPID() == 0)
    {
        if (saveState_)
        {
            std::stringstream ss;
            ss << outputFile_ << append;
            INFO("copying " << outputFile_ << " to " << ss.str());
            std::ifstream src(outputFile_.c_str(), std::ios::binary);
            std::ofstream dst(ss.str(), std::ios::binary);
            dst << src.rdbuf();
        }
        else
        {
            WARNING("No use in copying a state when saveState = false",
                    __FILE__, __LINE__);
        }
    }
    return 0;
}

//=============================================================================
void Model::gid2coord(int const &gid, int &mdl,
                             int &i, int &j, int &k, int &xx)
{
    mdl = modelIdent();

    int N   = getDomain()->GlobalN();
    int M   = getDomain()->GlobalM();
    int L   = getDomain()->GlobalL();
    int dof = getDomain()->Dof();
    int dim = N*M*L*dof;

    int aux = gid - dim;
    if (aux >= 0) // this is an auxiliary unknown
    {
        xx = dof + aux;
        i  = 0;
        j  = 0;
        k  = 0;
    }
    else
    {
        int tmp = gid;

        xx  = tmp % dof;
        tmp = (tmp - xx) / dof;
        i   = tmp % N;
        tmp = (tmp - i) / N;
        j   = tmp % M;
        k   = (tmp - j) / M;
    }
}

void Model::initializeState()
{
    state_->PutScalar(0.0);
}
