#include "GlobalDefinitions.H"

#include <ctime>  // std::clock()
#include <fstream>

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
    double time = timerStack.top().ElapsedTime();
    bool sane = true;
    if (!(msg == timerStack.top().label()))
    {
        WARNING("msg and label not equal!\n" <<
                "   msg = " << msg << "\n"
                " label = " << timerStack.top().label() << "\n",
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
