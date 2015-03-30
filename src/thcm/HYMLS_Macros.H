#ifndef HYMLS_MACROS_H
#define HYMLS_MACROS_H

#include "Teuchos_StandardCatchMacros.hpp"

#ifdef DEBUGGING
#ifndef TESTING
#define TESTING
#endif
#endif

#ifdef USE_SCOREP
#include "SCOREP_User.h"
#else
#define SCOREP_USER_REGION(s,t)
#define SCOREP_USER_PARAMETER_INT64(s,p)
#endif

#ifndef HYMLS_SMALL_ENTRY
#include <limits>
//#define HYMLS_SMALL_ENTRY std::numeric_limits<double>::epsilon()
#define HYMLS_SMALL_ENTRY 1.0e-14
#endif

// the timing macro's are also used for function tracing, so enable them all:
#ifdef DEBUGGING
#ifndef TIMING_LEVEL
#define TIMING_LEVEL 10
#endif
#ifndef PRINT_TIMING
#define PRINT_TIMING false
#endif
#else
#ifndef TIMING_LEVEL
#define TIMING_LEVEL 2
#endif
#ifndef PRINT_TIMING
#define PRINT_TIMING false
#endif
#endif

#ifdef TESTING
#ifndef FUNCTION_TRACING
#define FUNCTION_TRACING
#endif
#endif

// the TIMING_LEVEL macro defines how many functions are included
// in the timing output at the end of the run.
//
// 0: no timing at all
// 1: timing of major functions like preconditioner setup
// 2: detailed timing of things like applying orthogonal transformations
// 3: time everything, this is mainly for the case where TESTING/DEBUGGING
//    is defined and the timing routines are abused for function tracing. 
//    setting this may infringe the performance substantially.


// the HYMLS_PROF macros can be placed anywhere in the code,
// but only once in each scope ({}-block). The
// timing is stopped at the end of the current scope. So for 
// instance, to time a function foo and a part of it as well,
// do something like this:
// void bar::foo()
//   {
//   HYMLS_PROF("bar","foo");
//   ... do something that is not timed separately ...
//   {
//   HYMLS_PROF2("bar","foo-part")
//      ... do foo-part
//   }
//   ... do some more work that is not timed separately
//  }


#if TIMING_LEVEL>0
/*! @def HYMLS_PROF - start profiling region

 s1 and s2 are concatenated to form the profiler/timer label.
 s1 may be e.g. an object label and s2 a function name.
 */
#define HYMLS_PROF(s1,s2) HYMLS::TimerObject \
Error_You_are_trying_to_start_multiple_timers_in_one_scope \
(((std::string)(s1))+": "+(((std::string)(s2))),PRINT_TIMING); \
SCOREP_USER_REGION((std::string(s1)+std::string(s2)).c_str(),SCOREP_USER_REGION_TYPE_FUNCTION)

/*! @def HYMLS_LPROF(s1,s2): like HYMLS_PROF, but for HYMLS' recursively constructed 
classes. The macro assumes that the calling scope (e.g. the class) has an int variable 
'myLevel_', which it will append to s1.*/ 
#define HYMLS_LPROF(s1,s2) HYMLS_PROF(s1+"_L"+Teuchos::toString(myLevel_),s2) \
SCOREP_USER_PARAMETER_INT64("level",(int64_t)myLevel_);
#else
#define HYMLS_PROF(s1,s2)
#define HYMLS_LPROF(s1,s2)
#endif

#ifdef DEBUGGING
#define SET_CHECKPOINT(s1,s2,s3,s4,s5) \
HYMLS::Tools::SetCheckPoint(((std::string)s1) + ": "+((std::string)s2),s3,s4,s5);
#define BREAK_ON_CHECKPOINT(s1,s2) {std::string msg; std::string file; int line; \
if (HYMLS::Tools::GetCheckPoint(s1+": "+s2,msg,file,line)){ \
HYMLS::Tools::Fatal("aborting at check point '"+msg+"'",file.c_str(),line);}}
#else
#define SET_CHECKPOINT(s1,s2,s3,s4,s5)
#define BREAK_ON_CHECKPOINT(s1,s2)
#endif

#if TIMING_LEVEL>1
#define HYMLS_PROF2(s1,s2) HYMLS_PROF(s1,s2)
#define HYMLS_LPROF2(s1,s2) HYMLS_LPROF(s1,s2)
#else
#define HYMLS_PROF2(s1,s2)
#define HYMLS_LPROF2(s1,s2)
#endif

#if TIMING_LEVEL>2
#define HYMLS_PROF3(s1,s2) HYMLS_PROF(s1,s2)
#define HYMLS_LPROF3(s1,s2) HYMLS_LPROF(s1,s2)
#else
#define HYMLS_PROF3(s1,s2)
#define HYMLS_LPROF3(s1,s2)
#endif

#ifndef NO_MEMORY_TRACING
#define REPORT_SUM_MEM(label,description,doubles,ints,comm) \
{ \
double lval=(double)((doubles)*sizeof(double)+(ints)*sizeof(int)); \
double gval=lval; \
(comm)->SumAll(&lval,&gval,1); \
HYMLS::Tools::ReportMemUsage((std::string)(label)+": "+(std::string)(description), gval); \
}
#define REPORT_MEM(label,description,doubles,ints) \
{ \
double gval=(double)((doubles)*sizeof(double)+(ints)*sizeof(int)); \
HYMLS::Tools::ReportMemUsage((std::string)(label)+": "+(std::string)(description), gval); \
}
#else
#define REPORT_SUM_MEM(label,description,doubles,ints,comm)
#define REPORT_MEM(label,description,doubles,ints)
#endif


#ifdef DEBUGGING
#ifndef DEBUG
#define DEBUG(s) HYMLS::Tools::deb() << s << std::endl << std::flush;
#endif
#ifndef DEBVAR
#define DEBVAR(s) HYMLS::Tools::deb() << #s << " = "<<s << std::endl << std::flush;
#endif
#else
#ifndef DEBUG
#define DEBUG(s) 
#endif
#ifndef DEBVAR
#define DEBVAR(s) 
#endif
#endif

#ifndef CHECK_ZERO
#define CHECK_ZERO(funcall) {int ierr=0; bool status=true; try{\
ierr = funcall; } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,status)\
if (!status) {\
std::string msg="Caught an exception in call "+std::string(#funcall);\
std::cerr << std::flush; \
HYMLS::Tools::Fatal(msg,__FILE__,__LINE__);} \
if (ierr) {\
std::string msg="error code "+Teuchos::toString(ierr)+" returned from call "+std::string(#funcall);\
HYMLS::Tools::Error(msg,__FILE__,__LINE__);}\
}
#endif

#ifndef CHECK_TRUE
#define CHECK_TRUE(funcall) {bool berr = funcall;\
if (!berr) {std::cerr<<"Trilinos call "<<#funcall<<" returned false"<<std::endl;}}
#endif
#ifndef CHECK_NONNEG
#define CHECK_NONNEG(funcall) {int ierr = funcall;\
if (ierr<0) {std::cerr<<"Trilinos Error "<<ierr<<" returned from call "<<#funcall<<std::endl;}}
#endif

//! our own 'modulo' function, which behaves like mod in matlab.
//! the C++ built-in '%' operator returns -1%n=-1 and is therefore
//! not very useful for periodic boundaries...
#ifndef MOD
#define MOD(x,y) (((double)(y)==0.0)? (double)(x): ((double)(x) - floor((double)(x)/((double)(y)))*((double)(y))))
#endif

#ifndef PEC_INIT
#define PEC_INIT HYMLS::ParallelExceptionCatcher PEC
#endif

#ifndef PEC_PROTECT
#define PEC_PROTECT if (PEC.CheckErrors()) continue;
#endif

#ifndef PEC_IFPACK_CHK_ERR
#define PEC_IFPACK_CHK_ERR(x) PEC.IfpackCheckErr(x); PEC_PROTECT
//#define PEC_IFPACK_CHK_ERR(x) IFPACK_CHK_ERR(x)
#endif

#ifndef PEC_IFPACK_CHK_ERR2
#define PEC_IFPACK_CHK_ERR2(x) PEC.IfpackCheckErr(x)
//#define PEC_IFPACK_CHK_ERR2(x) IFPACK_CHK_ERR(x)
#endif

#ifndef PEC_ERROR
#define PEC_ERROR(x___) PEC.Error(x___); PEC_PROTECT
#endif

#ifndef PEC_CATCH
#define PEC_CATCH(x___) try { x___;} catch (HYMLS::Exception &ex___) {PEC.Error(ex___);} PEC_PROTECT
//#define PEC_CATCH(x) x
#endif

#ifndef PEC_CATCH2
#define PEC_CATCH2(x___) try { x___;} catch (HYMLS::Exception &ex___) {PEC.Error(ex___);}
//#define PEC_CATCH2(x) x
#endif

#ifndef PEC_HANDLE_ERRORS
#define PEC_HANDLE_ERRORS {int err___ = PEC.HandleErrors(); if (err___ != 0) return err___;}
#endif

#endif
