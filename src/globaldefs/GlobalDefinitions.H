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
#include <mpi.h>

#include "Teuchos_StandardCatchMacros.hpp"

#include <iostream>
#include <cmath>
#include <map>
#include <array>
#include <string>
#include <stack>

// These outstreams need to be defined in the main routine.
extern Teuchos::RCP<std::ostream> outFile;
extern Teuchos::RCP<std::ostream> cdataFile;
extern Teuchos::RCP<std::ostream> tdataFile;

class Epetra_Comm;

Teuchos::RCP<Epetra_Comm> initializeEnvironment(
    int argc, char **argv,
    Teuchos::RCP<std::ostream> info = Teuchos::null,
    Teuchos::RCP<std::ostream> cdata = Teuchos::null,
    Teuchos::RCP<std::ostream> tdata = Teuchos::null);

//------------------------------------------------------------------
// Profile definition:

// Number of entries in profile
#ifndef PROFILE_ENTRIES
# define PROFILE_ENTRIES 4
#endif

typedef std::map<std::string, std::array<double, PROFILE_ENTRIES> > ProfileType;

//=========================================================================
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

// Precision in scientific output (even)
#ifndef _PRECISION_
# define _PRECISION_ 8
#endif

// Set field width in data outputs (even)
#ifndef _FIELDWIDTH_
# define _FIELDWIDTH_ 16
#endif

#ifndef INFO
#  define INFO(s) (*outFile) << s << std::endl;
#endif

#ifndef WRITECDATA
#  define WRITECDATA(s) (*cdataFile) << s << std::endl;
#endif

#ifndef WRITETDATA
#  define WRITETDATA(s) (*tdataFile) << s << std::endl; 
#endif

// ! Error reporting
#ifndef ERROR
#  define ERROR(x,y,z)                                  \
    {                                                   \
    std::cerr << "\n**ERROR**: "<< x << std::endl;      \
    std::cerr << "(in " << y << ", line "               \
    << z << ")\n" << std::endl;                         \
                                                        \
    (*outFile) << "\n**ERROR**: "<< x << std::endl;     \
    (*outFile) << "(in " << y << ", line "              \
    << z << ")\n" << std::endl;                         \
                                                        \
    int initialized = 0;                                \
    MPI_Initialized(&initialized);                      \
    if (initialized)                                    \
        MPI_Abort(MPI_COMM_WORLD, -1);                  \
    std::exit(-1);                                      \
    }
#endif

#ifndef WARNING
#  define WARNING(x,y,z)                                \
    {                                                   \
        std::cerr << "\n**WARNING**: " << x << " ";     \
        std::cerr << "(in " << y << ", line "           \
                  << z << ")\n" << std::endl;           \
        (*outFile) << "\n**WARNING**: " << x << " ";    \
        (*outFile) << "(in " << y << ", line "          \
                   << z << ")\n" << std::endl;          \
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
#  define DUMP_VECTOR(fname,vec) {                                      \
        Teuchos::RCP<Epetra_MultiVector> mVec = Utils::Gather(vec, 0);  \
        if (mVec->Map().Comm().MyPID() == 0)                            \
        {                                                               \
            std::ofstream file;                                         \
            file.open(fname);                                           \
            for (int i = 0; i != mVec->MyLength(); ++i)                 \
            {                                                           \
                file << std::setprecision(17) << std::scientific        \
                     << (*(*mVec)(0))[i] << "\n";                       \
            }                                                           \
        }                                                               \
    }                                                                   \
 
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

//------------------------------------------------------------------
// Timer class, compatible with macros here in GlobalDefinitions.H
class Timer
{
    double startTime_;
    std::string label_;
public:
    Timer(std::string const &label);
    double wallTime();
    void ResetStartTime();
    double ElapsedTime();
    std::string label();
};

class ScopeTimer
{
    std::string label_;
public:
    ScopeTimer(std::string const &label);
    ~ScopeTimer();
};

void printProfile();

extern "C" {
void timer_start_(char const *msg);
void timer_stop_(char const *msg);

void track_iterations_(char const *charmsg, int iters);
void track_residual_(char const *charmsg, double residual);
}

// Timer macro using a global stack, allows for nesting
#ifndef TIMER_START
#  define TIMER_START(msg) timer_start_(msg);
#endif

#ifndef TIMER_STOP
#  define TIMER_STOP(msg) timer_stop_(msg);
#endif

#ifndef TIMER_SCOPE
#  define TIMER_SCOPE(msg) ScopeTimer Trying_to_start_multiple_timers_in_the_same_scope(msg);
#endif

// We can also keep track of different things, to distinguish these
// from timings we prepend NOTIME to the message
#ifndef TRACK_ITERATIONS
# define TRACK_ITERATIONS(msg, iters) track_iterations_(msg, iters);
#endif

#ifndef TRACK_RESIDUAL
# define TRACK_RESIDUAL(msg, residual) track_residual_(msg, residual);
#endif

//=========================================================================

// macros for checking Epetra calls
#ifndef CHECK_ZERO
#define CHECK_ZERO(funcall) {                                           \
        int ierr = 0;                                                   \
        bool status = true;                                             \
        try {                                                           \
            ierr = funcall;                                             \
        } TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, status);   \
        if (!status) {                                                  \
            ERROR("Caught an exception in call " << #funcall,           \
                  __FILE__, __LINE__);                                  \
        }                                                               \
        if (ierr) {                                                     \
            ERROR("Trilinos Error " << ierr << " returned from call "   \
                  << #funcall, __FILE__, __LINE__);                     \
        }                                                               \
    }                                                                   \

#define CHECK_TRUE(funcall) {bool berr = funcall;                       \
        if (!berr) {ERROR("Trilinos call " << #funcall                  \
                          << " returned false", __FILE__, __LINE__);}}  \


#define CHECK_NONNEG(funcall) {int ierr = funcall;                      \
        if (ierr<0) { ERROR("Trilinos call " << #funcall                \
                            << " returned negative: "                   \
                            << ierr , __FILE__, __LINE__); }}           \

#endif

// macro for checking and replacing maps
#ifndef CHECK_MAP
#define CHECK_MAP(vec, map) {                                           \
        if (!(vec->Map().SameAs(*map)))                                 \
        {                                                               \
            CHECK_ZERO(vec->ReplaceMap(*map));                          \
        }}
#endif

#ifndef MOD

//! our own 'modulo' function, which behaves like mod in matlab.
//! the C++ built-in '%' operator returns -1%n=-1 and is therefore
//! not very useful for periodic boundaries...
#  define MOD(x,y) (((double)(y)==0.0)? (double)(x): ((double)(x) - floor((double)(x)/((double)(y)))*((double)(y))))
#endif

#endif
