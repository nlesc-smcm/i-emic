/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

#include <fstream>
#include "Teuchos_RCP.hpp"
Teuchos::RCP<std::ofstream> info;
#ifdef DEBUGGING
Teuchos::RCP<std::ofstream> debug;
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

void Error(std::string msg, std::string file, int line)
{
	std::cerr << "ERROR: "   << msg  << std::endl;
	std::cerr <<         "(" << file << ", line " << line << ")\n";
	(*info)   << "ERROR: "   << msg  << std::endl;
	(*info)   <<         "(" << file << ", line " << line << ")\n";
#ifdef HAVE_MPI
	MPI_Abort(MPI_COMM_WORLD, -1);
#endif
	exit(-1);
}
