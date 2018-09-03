#include "TimeStepper.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Utils.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "EpetraExt_HDF5.h"

template<>
void TimeStepper<Teuchos::RCP<Epetra_Vector> >::read(
    std::string const &name,
    std::vector<AMSExperiment<Teuchos::RCP<Epetra_Vector> > >  &experiments) const
{
    std::cout << "Reading state from " << name << std::endl;
    int experiments_size = experiments.size();
    if (experiments_size == 0)
        return;

    EpetraExt::HDF5 HDF5(experiments[0].x0->Comm());
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
                experiments[i].xlist[j] = Teuchos::rcp(new Epetra_Vector(map));
                experiments[i].xlist[j]->Import(*tmp, *import, Insert);
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
void TimeStepper<Teuchos::RCP<Epetra_Vector> >::write(
    std::string const &name,
    std::vector<AMSExperiment<Teuchos::RCP<Epetra_Vector> > > const &experiments) const
{
    std::cout << "Writing current state to " << name << std::endl;
    if (experiments.size() == 0 || experiments[0].xlist.size() == 0)
        return;

    EpetraExt::HDF5 HDF5(experiments[0].xlist[0]->Comm());
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
}
