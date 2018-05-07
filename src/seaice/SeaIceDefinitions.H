#ifndef SEAICE_DEFS
#define SEAICE_DEFS

#ifndef PI_
# define PI_ 3.14159265358979323846
#endif

// Sea ice is calculated without using any neighbours:
#ifndef SEAICE_NP_
# define SEAICE_NP_ 1
#endif

// Number of unknowns
#ifndef SEAICE_NUN_
# define SEAICE_NUN_ 3
#endif

// The unknowns we have in this model 
#ifndef SEAICE_HH_
# define SEAICE_HH_  1 // Atmospheric temperature
# define SEAICE_QQ_  2 // Vapor anomaly
# define SEAICE_MM_  3 // Vapor anomaly
#endif

#endif
