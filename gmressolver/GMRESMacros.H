#ifndef GMRESMACROS_H
#define GMRESMACROS_H

// PARALLEL PRINT MACRO
#ifndef PRINT
# ifdef HAVE_MPI
#  define PRINT(msg, verbosity)						\
	{												\
		if (verbosity > 0)							\
		{											\
			int rank;								\
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);	\
			if (rank == 0)							\
				std::cout << msg << std::endl;		\
		}											\
	}
# else
#  define PRINT(msg) std::cout << msg << std::endl;
# endif
#endif

#endif
