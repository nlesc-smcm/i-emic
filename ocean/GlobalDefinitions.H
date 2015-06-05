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

#include <iostream>
#include <cmath>
#include "Teuchos_RCP.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include <map>
#include <string>

// This outstream needs to be defined in the main routine.
extern Teuchos::RCP<std::ostream> outFile;

extern std::map<std::string, double> profile;

#ifndef INFO
#  define INFO(s) (*outFile) << s << std::endl; 
#endif

#ifndef ERROR
#  define ERROR(x,y,z) (*outFile) << "**ERROR**: " << x << " "\
	                              << y << " " << z << std::endl;
#endif

#ifndef WARNING
#  define WARNING(x,y,z) (*outFile) << "**WARNING**: " << x << " "\
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

#ifndef TIMER_START
#  define TIMER_START(message, timerPtr)		\
	(*outFile) << message << std::endl;			\
	timerPtr->ResetStartTime();					\
	profile[message] = (profile.count(message)) ? profile[message] : 0;
#endif

#ifndef TIMER_END
#  define TIMER_END(message, timerPtr)	{							\
	double time = timerPtr->ElapsedTime();							\
	profile[message] += time;										\
	(*outFile) << message << " done ("								\
			   << time << ", total: "								\
			   << profile[message] << ")" <<  std::endl; }
#endif
	
//==================================================================
//==================================================================

// macros for checking Epetra calls
#ifndef CHECK_ZERO
#define CHECK_ZERO(funcall) {int ierr = funcall;						\
		if (ierr) {INFO("Trilinos Error " << ierr << " returned from call " \
						                          << #funcall);\
			/*throw std::runtime_error("Trilinos Error...");*/}}
#define CHECK_TRUE(funcall) {bool berr = funcall;						\
		if (!berr) {INFO("Trilinos call " << #funcall << " returned false"); \
			/*throw std::runtime_error("Trilinos Error...");*/}}
#define CHECK_NONNEG(funcall) {int ierr = funcall;						\
		if (ierr<0) {INFO("Trilinos Error " << ierr << " returned from call " \
						                    << #funcall);\
			/*throw std::runtime_error("Trilinos Error...");*/}}
#endif

#ifndef MOD
//! our own 'modulo' function, which behaves like mod in matlab.
//! the C++ built-in '%' operator returns -1%n=-1 and is therefore
//! not very useful for periodic boundaries...
#  define MOD(x,y) (((double)(y)==0.0)? (double)(x): ((double)(x) - floor((double)(x)/((double)(y)))*((double)(y))))
#endif

#endif
