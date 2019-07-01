#include "GlobalDefinitions.H"

#include <ctime>  // std::clock()
#include <fstream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Epetra_config.h>

#  ifdef HAVE_MPI
// Epetra's wrapper for MPI_Comm. This header file only exists if
// Epetra was built with MPI enabled.
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif // HAVE_MPI

Teuchos::RCP<std::ostream> outFile;      // output file
Teuchos::RCP<std::ostream> cdataFile;    // cdata file
Teuchos::RCP<std::ostream> tdataFile;    // tdata file

//-----------------------------------------------------------------------------
void outputFiles(Teuchos::RCP<Epetra_Comm> Comm,
                 Teuchos::RCP<std::ostream> info,
                 Teuchos::RCP<std::ostream> cdata,
                 Teuchos::RCP<std::ostream> tdata)
{
    // Setup output files "fname_#.txt" for P==0 && P==1, other processes
    // will get a blackholestream.

    if (info != Teuchos::null)
        outFile = info;
    else if (Comm->MyPID() < 10)
    {
        std::ostringstream infofile;   // setting up a filename

        infofile  << "info_" << Comm->MyPID() << ".txt";

        std::cout << "info for CPU" << Comm->MyPID() << " is written to "
                  << infofile.str().c_str() << std::endl;

        outFile =
            Teuchos::rcp(new std::ofstream(infofile.str().c_str()));
    }
    else
        outFile = Teuchos::rcp(new Teuchos::oblackholestream());

    // Write continuation data on proc 0
    if (cdata != Teuchos::null)
        cdataFile = cdata;
    else if (Comm->MyPID() == 0)
    {
        std::ostringstream cdatafile;  // setting up a filename

        cdatafile << "cdata.txt";

        std::cout << "continuation data is written to " <<  cdatafile.str().c_str()
                  << " by CPU " << Comm->MyPID() << std::endl;

        cdataFile =
            Teuchos::rcp(new std::ofstream(cdatafile.str().c_str()));
    }
    else
        cdataFile = Teuchos::rcp(new Teuchos::oblackholestream());

    // Write timestepping data on proc 0
    if (tdata != Teuchos::null)
        tdataFile = tdata;
    else if (Comm->MyPID() == 0)
    {
        std::ostringstream tdatafile;  // setting up a filename

        tdatafile << "tdata.txt";

        std::cout << "continuation data is written to " << tdatafile.str().c_str()
                  << " by CPU " << Comm->MyPID() << std::endl;

        tdataFile =
            Teuchos::rcp(new std::ofstream(tdatafile.str().c_str()));
    }
    else
        tdataFile = Teuchos::rcp(new Teuchos::oblackholestream());
}

//-----------------------------------------------------------------------------
Teuchos::RCP<Epetra_Comm> initializeEnvironment(
    int argc, char **argv,
    Teuchos::RCP<std::ostream> info,
    Teuchos::RCP<std::ostream> cdata,
    Teuchos::RCP<std::ostream> tdata)
{
    // Setup MPI communicator

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Teuchos::RCP<Epetra_MpiComm> Comm =
        Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
    Teuchos::RCP<Epetra_SerialComm> Comm =
        Teuchos::rcp(new Epetra_SerialComm() );
#endif

    // Specify output files
    outputFiles(Comm, info, cdata, tdata);
    return Comm;
}

//------------------------------------------------------------------
// Setup profile:

// This profile container needs to be defined in the main routine.
ProfileType profile;

// We define a global stack for Timer objects, so that we can nest timings
std::stack<Timer> timerStack;

//------------------------------------------------------------------

void track_iterations_(char const *charmsg, int iters)
{
    std::string msg(charmsg);
    msg.insert(0, "_NOTIME_");
    profile[msg] = (profile.count(msg)) ?
        profile[msg] :
        std::array<double, PROFILE_ENTRIES>();
    profile[msg][0] += iters;
    profile[msg][1] += 1;
    profile[msg][2] = profile[msg][0] / profile[msg][1];
}

void track_residual_(char const *charmsg, double residual)
{
    std::string msg(charmsg);
    msg.insert(0, "_NOTIME_");
    profile[msg] = (profile.count(msg)) ?
        profile[msg] :
        std::array<double, PROFILE_ENTRIES>();
    profile[msg][0] += residual;
    profile[msg][1] += 1;
    profile[msg][2] = profile[msg][0] / profile[msg][1];
}

void timer_start_(char const *msg)
{
    Timer timer(msg);
    timer.ResetStartTime();
    if (profile.find(msg) == profile.end())
        profile[msg] = std::array<double, PROFILE_ENTRIES>();
    timerStack.push(timer);
}

void timer_stop_(char const *msg)
{
    std::string label = timerStack.top().label();
    double time = timerStack.top().ElapsedTime();
    bool sane = true;
    if (msg != label)
    {
        WARNING("Timer msg and label not equal!\n" <<
                "   msg = " << msg << "\n"
                " label = " << label << "\n",
                __FILE__, __LINE__);
        sane = false;
    }
    assert(sane);
    timerStack.pop();
    profile[msg][0] += time;
    profile[msg][1] += 1;
    profile[msg][2] = profile[msg][0] / profile[msg][1];
}

//------------------------------------------------------------------
// Timer class implementation
Timer::Timer(std::string const &label)
    :
    startTime_(0.0),
    label_(label)
{}

double Timer::wallTime()
{
    int mpiInit;
    MPI_Initialized(&mpiInit);
    if (mpiInit)
        return MPI_Wtime();
    else
        return (double) std::clock() / CLOCKS_PER_SEC;
}

void Timer::ResetStartTime()
{
    startTime_ = wallTime();
}

double Timer::ElapsedTime()
{
    return (double) (wallTime() - startTime_);
}

std::string Timer::label()
{
    return label_;
}

ScopeTimer::ScopeTimer(std::string const &label)
    :
    label_(label)
{
    timer_start_(label.c_str());
}

ScopeTimer::~ScopeTimer()
{
    timer_stop_(label_.c_str());
}

//-----------------------------------------------------------------------------
void printProfile()
{
    bool sane = true;
    while (timerStack.empty() == false)
    {
        WARNING("Unequal amount of TIMER_START and TIMER_STOP uses" <<
                " label = " << timerStack.top().label(),
                __FILE__, __LINE__);
        timerStack.pop();
        sane = false;
    }
    assert(sane);

    std::ostringstream profilefile("profile_output");   // setting up a filename
    std::ofstream file(profilefile.str().c_str());      // setup output file

    // Set format flags
    file << std::left;

    // Define line format
#ifndef LINE
# define LINE(s1, s2, s3, s4, s5, s6, s7, s8, s9)               \
    {                                                           \
        int sp = 3;  int it = 5;  int id = 5;                   \
        int db = 12; int st = 45;                               \
        file << std::setw(id) << s1     << std::setw(sp) << s2  \
             << std::setw(st) << s3 << std::setw(sp) << s4      \
             << std::setw(db) << s5     << std::setw(sp) << s6  \
             << std::setw(it) << s7     << std::setw(sp) << s8  \
             << std::setw(db) << s9     << std::endl;           \
    }
#endif

    // Header
    LINE("", "", "", "", "cumul.", "", "calls", "", "average");

    // Display timings of the separate models, summing
    int counter = 0;
    for (auto const &map : profile)
        if (map.first.compare(0,8,"_NOTIME_") != 0)
        {
            counter++;
            std::stringstream s;
            s << " (" << counter << ")";
            LINE(s.str(), "", map.first, ":", map.second[0], "",
                 map.second[1], "", map.second[2]);
        }

    // Newline
    file << std::endl;

    // Display iteration information
    for (auto const &map : profile)
        if (map.first.compare(0,8,"_NOTIME_") == 0 )
        {
            counter++;
            std::stringstream s;
            s << " (" << counter << ")";
            LINE(s.str(), "", map.first.substr(8), ":", map.second[0], "",
                 map.second[1], "", map.second[2]);
        }

}
