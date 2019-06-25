#ifndef STOCHASTICBASE_H
#define STOCHASTICBASE_H

#include <random>

#include "Epetra_Comm.h"

class StochasticBase
{
protected:
    //! RNG
    unsigned int noise_seed_;
    Teuchos::RCP<std::mt19937_64> engine_;

public:
    //-------------------------------------------------------
    //! constructor
    template<typename ParameterList>
    StochasticBase(Epetra_Comm const &comm, ParameterList params)
        :
        noise_seed_(params->get("noise seed", 0))
        {
            std::random_device::result_type seed = noise_seed_;
            if (seed == 0)
            {
                std::random_device rd;
                seed = rd();
            }
            std::seed_seq seeder{seed};
            engine_ = Teuchos::rcp(new std::mt19937_64(seeder));

            write_seed(comm, seed, "noise seed");
        }

    virtual ~StochasticBase() {}

    static void write_seed(Epetra_Comm const &comm, unsigned int seed, std::string const &label)
        {
            unsigned int *seeds = new unsigned int[comm.NumProc()];

            int *seedptr = reinterpret_cast<int *>(&seed);
            int *seedsptr = reinterpret_cast<int *>(seeds);

            CHECK_ZERO(comm.GatherAll(seedptr, seedsptr, 1));
            for (int i = 0; i < comm.NumProc(); i++)
                INFO(label << ": " << seeds[i]);
            delete[] seeds;
        }
};

#endif
