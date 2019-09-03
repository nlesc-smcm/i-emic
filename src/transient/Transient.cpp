#include "Transient.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Utils.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "EpetraExt_HDF5.h"

#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>

#include "Trilinos_version.h"

std::string mem2string(long long mem)
{
    double value = mem;
    std::string unit = "B";
    if (std::abs(value) > 1.0e3) {value *= 1.0e-3; unit = "kB";}
    if (std::abs(value) > 1.0e3) {value *= 1.0e-3; unit = "MB";}
    if (std::abs(value) > 1.0e3) {value *= 1.0e-3; unit = "GB";}
    if (std::abs(value) > 1.0e3) {value *= 1.0e-3; unit = "TB";}

    std::ostringstream ss;
    ss << std::fixed;
    ss.precision(2);
    ss << value << " " << unit;
    return ss.str();
}

// This read/write mechanism need Trilinos pull request #3381, which
// is present in Trilinos 12.14
#if TRILINOS_MAJOR_MINOR_VERSION > 121300

template<>
void Transient<Teuchos::RCP<const Epetra_Vector> >::read(
    std::string const &name,
    std::vector<AMSExperiment<Teuchos::RCP<const Epetra_Vector> > > &experiments) const
{
    if (name == "")
        return;

    std::cout << "Reading state from " << name << std::endl;
    int experiments_size = experiments.size();
    if (experiments_size == 0)
        return;

    Epetra_Comm const &comm = experiments[0].x0->Comm();
    EpetraExt::HDF5 HDF5(comm);
    HDF5.Open(name);

    int num_exp = -1;
    HDF5.Read("data", "num exp", num_exp);
    if (num_exp_ != num_exp)
    {
        std::cerr << "Error: Number of experiments is " << num_exp_
                  << " instead of " << num_exp << std::endl;
        return;
    }

    int num_init_exp = -1;
    HDF5.Read("data", "num init exp", num_init_exp);
    if (experiments_size != num_exp && experiments_size < num_init_exp)
    {
        std::cerr << "Error: Number of experiments is " << experiments_size
                  << " instead of " << num_exp << " or " << num_init_exp << std::endl;
        return;
    }

    HDF5.Read("data", "its", its_);

    int ell_size = -1;
    HDF5.Read("data", "ell size", ell_size);
    if (ell_size > 0)
    {
       ell_.resize(ell_size);
       HDF5.Read("data", "ell", H5T_NATIVE_INT, ell_size, &ell_[0]);
    }

    Teuchos::RCP<Epetra_Import> import = Teuchos::null;
    Epetra_MultiVector *tmp;
    Epetra_BlockMap const &map = experiments[0].x0->Map();

    for (int i = 0; i < experiments_size; i++)
    {
        int size = -1;
        HDF5.Read("experiments/" + Teuchos::toString(i),
                  "size", size);
        HDF5.Read("experiments/" + Teuchos::toString(i),
                  "max distance", experiments[i].max_distance);
        HDF5.Read("experiments/" + Teuchos::toString(i),
                  "time", experiments[i].time);
        HDF5.Read("experiments/" + Teuchos::toString(i),
                  "initial time", experiments[i].initial_time);
        HDF5.Read("experiments/" + Teuchos::toString(i),
                  "return time", experiments[i].return_time);
        int converged = 0;
        HDF5.Read("experiments/" + Teuchos::toString(i),
                  "converged", converged);
        experiments[i].converged = converged;
        int initialized = 0;
        HDF5.Read("experiments/" + Teuchos::toString(i),
                  "initialized", initialized);
        experiments[i].initialized = initialized;

        if (size > 0)
        {
            if (import == Teuchos::null)
            {
                HDF5.Read("experiments/" + Teuchos::toString(i) +
                          "/xlist/0", tmp);
                import = Teuchos::rcp(new Epetra_Import(map, tmp->Map()));
                delete tmp;
            }

            experiments[i].xlist.resize(size);
            experiments[i].dlist.resize(size);
            experiments[i].tlist.resize(size);

            for (int j = 0; j < size; j++)
            {
                HDF5.Read("experiments/" + Teuchos::toString(i) +
                          "/xlist/" + Teuchos::toString(j), tmp);
                Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(map));
                x->Import(*tmp, *import, Insert);
                experiments[i].xlist[j] = x;
                delete tmp;
            }

            HDF5.Read("experiments/" + Teuchos::toString(i), "dlist",
                      H5T_NATIVE_DOUBLE, size, &experiments[i].dlist[0]);

            HDF5.Read("experiments/" + Teuchos::toString(i), "tlist",
                      H5T_NATIVE_DOUBLE, size, &experiments[i].tlist[0]);
        }
    }
}

template<>
void Transient<Teuchos::RCP<const Epetra_Vector> >::write(
    std::string const &name,
    std::vector<AMSExperiment<Teuchos::RCP<const Epetra_Vector> > > const &experiments) const
{
    if (name == "")
        return;

    std::cout << "Writing current state to " << name << std::endl;
    if (experiments.size() == 0 || experiments[0].xlist.size() == 0)
        return;

    Epetra_Comm const &comm = experiments[0].x0->Comm();
    int lock_file = -1;
    if (comm.MyPID() == 0)
    {
        lock_file = open("lock", O_RDWR | O_CREAT, 0666);
        flock(lock_file, LOCK_EX);
    }
    comm.Barrier();

    EpetraExt::HDF5 HDF5(comm);
    HDF5.Create(name);
    HDF5.CreateGroup("experiments");
    for (int i = 0; i < (int)experiments.size(); i++)
    {
        HDF5.CreateGroup("experiments/" + Teuchos::toString(i));
        int size = experiments[i].dlist.size();
        HDF5.Write("experiments/" + Teuchos::toString(i),
                   "size", size);
        HDF5.Write("experiments/" + Teuchos::toString(i),
                   "max distance", experiments[i].max_distance);
        HDF5.Write("experiments/" + Teuchos::toString(i),
                   "time", experiments[i].time);
        HDF5.Write("experiments/" + Teuchos::toString(i),
                   "initial time", experiments[i].initial_time);
        HDF5.Write("experiments/" + Teuchos::toString(i),
                   "return time", experiments[i].return_time);
        HDF5.Write("experiments/" + Teuchos::toString(i),
                   "converged", experiments[i].converged);
        HDF5.Write("experiments/" + Teuchos::toString(i),
                   "initialized", experiments[i].initialized);

        if (size > 0)
        {
            HDF5.CreateGroup("experiments/" + Teuchos::toString(i) + "/xlist");
            for (int j = 0; j < (int)experiments[i].xlist.size(); j++)
                HDF5.Write("experiments/" + Teuchos::toString(i) +
                           "/xlist/" + Teuchos::toString(j), *experiments[i].xlist[j]);

            HDF5.Write("experiments/" + Teuchos::toString(i), "dlist",
                       H5T_NATIVE_DOUBLE, experiments[i].dlist.size(),
                       &experiments[i].dlist[0]);

            HDF5.Write("experiments/" + Teuchos::toString(i), "tlist",
                       H5T_NATIVE_DOUBLE, experiments[i].tlist.size(),
                       &experiments[i].tlist[0]);
        }
    }

    HDF5.Write("data", "its", its_);
    HDF5.Write("data", "num exp", num_exp_);
    HDF5.Write("data", "num init exp", num_init_exp_);

    int ell_size = ell_.size();
    HDF5.Write("data", "ell size", ell_size);
    if (ell_size > 0)
        HDF5.Write("data", "ell", H5T_NATIVE_INT, ell_.size(), &ell_[0]);

    if (lock_file >= 0)
        close(lock_file);
}
#endif //TRILINOS_MAJOR_MINOR_VERSION
