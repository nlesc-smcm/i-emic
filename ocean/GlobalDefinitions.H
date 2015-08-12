/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
/**********************************************************************
 * Modified by T.E. Mulder, Utrecht University 2014/15                *
 * contact: t.e.mulder@uu.nl                                          *
 **********************************************************************/

#ifndef GLOBALDEFINITIONS_H
#define GLOBALDEFINITIONS_H

#include <Teuchos_RCP.hpp>
#include <EpetraExt_RowMatrixOut.h>
#include <mpi.h>

#include <iostream>
#include <cmath>
#include <map>
#include <array>
#include <string>
#include <stack>
#include <ctime>  // std::clock()

// This outstream needs to be defined in the main routine.
extern Teuchos::RCP<std::ostream> outFile;

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

#ifndef INFO
#  define INFO(s) (*outFile) << s << std::endl; 
#endif

#ifndef ERROR
#  define ERROR(x,y,z) (*outFile) << "**ERROR**: " << x << " "	\
	<< y << " " << z << std::endl;
#endif

#ifndef WARNING
#  define WARNING(x,y,z) (*outFile) << "**WARNING**: " << x << "\n"	\
	<< y << " " << z << std::endl;
#endif

#ifndef DEBUG
#  ifdef DEBUGGING
#    define DEBUG(s) (*outFile) << s << std::endl;
// The use of # stringifies variable s: 
#    define DEBVAR(s) (*outFile) << #s << ": " << s << std::endl;
#  else
#    define DEBUG(s)
#    define DEBVAR(s)
#  endif
#endif

#ifndef DUMP
#  define DUMP(x,y) EpetraExt::RowMatrixToMatrixMarketFile(x, y);
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
#  define TIMER_START(msg)					\
	{												\
		Timer timer;								\
		timer.ResetStartTime();						\
		profile[msg] = (profile.count(msg)) ?		\
			profile[msg] :							\
			std::array<double, PROFILE_ENTRIES>();	\
		timerStack.push(timer);						\
	} 
#endif

#ifndef TIMER_STOP
#  define TIMER_STOP(msg)										\
	{															\
		double time = timerStack.top().ElapsedTime();			\
		timerStack.pop();										\
		profile[msg][0] += time;								\
		profile[msg][1] += 1;									\
		profile[msg][2] = profile[msg][0] / profile[msg][1];	\
	}
#endif

#ifndef TRACK_ITERATIONS
# define TRACK_ITERATIONS(msg, iters)							\
	{															\
		profile[msg] = (profile.count(msg)) ?					\
			profile[msg] :										\
			std::array<double, PROFILE_ENTRIES>();				\
		profile[msg][0] += iters;								\
		profile[msg][1] += 1;									\
		profile[msg][2] = profile[msg][0] / profile[msg][1];	\
	}
#endif	

//------------------------------------------------------------------
// Timer class, compatible with macros here in GlobalDefinitions.H
// Not used...
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
#define CHECK_ZERO(funcall) {int ierr = funcall;						\
		if (ierr) {INFO("Trilinos Error " << ierr << " returned from call " \
						<< #funcall);									\
			/*throw std::runtime_error("Trilinos Error...");*/}}
#define CHECK_TRUE(funcall) {bool berr = funcall;						\
		if (!berr) {INFO("Trilinos call " << #funcall << " returned false"); \
			/*throw std::runtime_error("Trilinos Error...");*/}}
#define CHECK_NONNEG(funcall) {int ierr = funcall;						\
		if (ierr<0) {INFO("Trilinos Error " << ierr << " returned from call " \
						  << #funcall);									\
			/*throw std::runtime_error("Trilinos Error...");*/}}
#endif

#ifndef MOD
//! our own 'modulo' function, which behaves like mod in matlab.
//! the C++ built-in '%' operator returns -1%n=-1 and is therefore
//! not very useful for periodic boundaries...
#  define MOD(x,y) (((double)(y)==0.0)? (double)(x): ((double)(x) - floor((double)(x)/((double)(y)))*((double)(y))))
#endif



#endif
