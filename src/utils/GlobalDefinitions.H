/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
/**********************************************************************
 * Modified/Extended by T.E. Mulder, Utrecht University 2014/15/16    *
 * contact: t.e.mulder@uu.nl                                           *
 **********************************************************************/

#ifndef GLOBALDEFINITIONS_H
#define GLOBALDEFINITIONS_H

#include <Teuchos_RCP.hpp>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>
#include <mpi.h>

#include <iostream>
#include <cmath>
#include <map>
#include <array>
#include <string>
#include <stack>
#include <ctime>  // std::clock()

// These outstreams need to be defined in the main routine.
extern Teuchos::RCP<std::ostream> outFile;
extern Teuchos::RCP<std::ostream> cdataFile;

//------------------------------------------------------------------
// Setup profile:

// Number of entries in profile
#ifndef PROFILE_ENTRIES
# define PROFILE_ENTRIES 4
#endif

typedef std::map<std::string, std::array<double, PROFILE_ENTRIES> > ProfileType;

// This profile container needs to be defined in the main routine.
extern ProfileType profile;

// Forward declaration of Timer class
class Timer;

// We define a global stack for Timer objects, so that we can nest timings
extern std::stack<Timer> timerStack;

//=========================================================================
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

// Precision in scientific output (even)
#ifndef _PRECISION_
# define _PRECISION_ 6
#endif

// Set field width in data outputs (even)
#ifndef _FIELDWIDTH_
# define _FIELDWIDTH_ 14
#endif 

#ifndef INFO
#  define INFO(s) (*outFile) << s << std::endl; 
#endif

#ifndef WRITECDATA
#  define WRITECDATA(s) (*cdataFile) << s << std::endl; 
#endif

//! Fatal error handling
#ifndef ERROR
#  define ERROR(x,y,z)                                                  \
	{                                                                   \
		(*outFile) << "**ERROR**: " << x << " "                         \
				   << y << " " << z << std::endl;                       \
                                                                        \
        std::cerr << "Fatal Error: "<< x << std::endl;                  \
        std::cerr << "(in " << y << ", line "                           \
                  << z << ")" << std::endl;                             \
                                                                        \
        int initialized = 0;                                            \
        MPI_Initialized(&initialized);                                  \
        if (initialized)                                                \
            MPI_Abort(MPI_COMM_WORLD, -1);                              \
        std::exit(-1);                                                  \
}
#endif

#ifndef WARNING
#  define WARNING(x,y,z)							\
	{												\
		(*outFile) << "**WARNING**: " << x << " "	\
				   << y << " " << z << std::endl;	\
	}
#endif

#ifndef DEBUG
#  ifdef DEBUGGING_OLD
#    define DEBUG(s) (*outFile) << s << std::endl;
// # stringifies variable s: 
#    define DEBVAR(s) (*outFile) << #s << ": " << s << std::endl;
#  else
#    define DEBUG(s)
#    define DEBVAR(s)
#  endif
#endif

// For new debug statements
#ifndef DEB
#  ifdef DEBUGGING_NEW
#    define DEB(s) (*outFile) << s << std::endl;
// # stringifies variable s: 
#    define DEBV(s) (*outFile) << #s << ": " << s << std::endl;
#  else
#    define DEB(s)
#    define DEBV(s)
#  endif
#endif

// Don't trust these dumps!! I think calling print on such objects is safer.
#ifndef DUMP 
#  define DUMP(fname,mat) EpetraExt::RowMatrixToMatrixMarketFile(fname, mat);
#endif

#ifndef DUMP_VECTOR
#  define DUMP_VECTOR(fname,mat) EpetraExt::VectorToMatlabFile(fname, mat);
#endif

#ifndef DUMPMATLAB
#  define DUMPMATLAB(fname,mat) EpetraExt::RowMatrixToMatlabFile(fname, mat);
#endif

#ifndef MAX
#  define MAX(x,y) (x > y ? x : y)
#endif

#ifndef MIN
#  define MIN(x,y) (x < y ? x : y)
#endif

#ifndef SGN
#  define SGN(x) (x < 0 ? -1 : 1)
#endif

// Timer macro using a global stack, allows for nesting
#ifndef TIMER_START
#  define TIMER_START(msg)											\
	{																\
		Timer timer;												\
		timer.ResetStartTime();										\
		if (profile.find(msg) == profile.end())   \
            profile[msg] = std::array<double, PROFILE_ENTRIES>(); \
		timerStack.push(timer);										\
	} 
#endif

#ifndef TIMER_STOP
#  define TIMER_STOP(msg)										\
	{															\
		double time = timerStack.top().ElapsedTime();			\
		timerStack.pop();										\
        profile[msg][0] += time;                               \
        profile[msg][1] += 1;                                  \
        profile[msg][2] = profile[msg][0] / profile[msg][1]; \
	}
#endif

// We can also keep track of different things, to distinguish these
// from timings we prepend NOTIME to the message
#ifndef TRACK_ITERATIONS
# define TRACK_ITERATIONS(charmsg, iters)						\
	{															\
		std::string msg(charmsg);								\
		msg.insert(0, "_NOTIME_");								\
        profile[msg] = (profile.count(msg)) ?         \
			profile[msg] :										\
			std::array<double, PROFILE_ENTRIES>();				\
        profile[msg][0] += iters;                              \
        profile[msg][1] += 1;                                  \
        profile[msg][2] = profile[msg][0] / profile[msg][1]; \
	}
#endif

#ifndef TRACK_RESIDUAL
# define TRACK_RESIDUAL(charmsg, residual)						\
	{															\
		std::string msg(charmsg);								\
		msg.insert(0, "_NOTIME_");								\
        profile[msg] = (profile.count(msg)) ?         \
			profile[msg] :										\
			std::array<double, PROFILE_ENTRIES>();				\
        profile[msg][0] += residual;                           \
        profile[msg][1] += 1;                                  \
        profile[msg][2] = profile[msg][0] / profile[msg][1]; \
	}
#endif

//------------------------------------------------------------------
// Timer class, compatible with macros here in GlobalDefinitions.H
class Timer
{
	double startTime_;
public:
	Timer()
		:
		startTime_(0.0)
		{}
	
	double wallTime()
		{
			int mpiInit;
			MPI_Initialized(&mpiInit);
			if (mpiInit)
				return MPI_Wtime();
			else
				return (double) std::clock() / CLOCKS_PER_SEC;
		}
	
	void ResetStartTime()
		{ startTime_ = wallTime();}
	
	double ElapsedTime()
		{ return (double) (wallTime() - startTime_); }
};

//=========================================================================

// macros for checking Epetra calls
#ifndef CHECK_ZERO
#define CHECK_ZERO(funcall) {                                           \
        int ierr = funcall;                                             \
        if (ierr) { ERROR("Trilinos Error " << ierr << " returned from call " \
                          << #funcall, __FILE__, __LINE__);             \
        }                                                               \
    }                                                                   \
 
#define CHECK_TRUE(funcall) {bool berr = funcall;                       \
		if (!berr) {ERROR("Trilinos call " << #funcall << " returned false", __FILE__, __LINE__);}} \


#define CHECK_NONNEG(funcall) {int ierr = funcall;                      \
        if (ierr<0) { ERROR("Trilinos Error " << ierr << " returned from call " \
                            << #funcall, __FILE__, __LINE__); }}       \
        
#endif

#ifndef MOD
//! our own 'modulo' function, which behaves like mod in matlab.
//! the C++ built-in '%' operator returns -1%n=-1 and is therefore
//! not very useful for periodic boundaries...
#  define MOD(x,y) (((double)(y)==0.0)? (double)(x): ((double)(x) - floor((double)(x)/((double)(y)))*((double)(y))))
#endif

#endif
