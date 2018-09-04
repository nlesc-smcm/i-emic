#ifndef ATMOSPHERE_DEFS
#define ATMOSPHERE_DEFS

// Pi
#ifndef PI_
# define PI_ 3.14159265358979323846
#endif

//! Number of points in neighbourhood:
//! +----------++-------++----------+
//! | 12 15 18 || 3 6 9 || 21 24 27 |
//! | 11 14 17 || 2 5 8 || 20 23 26 |
//! | 10 13 16 || 1 4 7 || 19 22 25 |
//! |  below   || center||  above   |
//! +----------++-------++----------+
#ifndef ATMOS_NP_
# define ATMOS_NP_ 27
#endif

// Number of ordinary unknowns (excluding auxiliary single-valued
// unknowns).
#ifndef ATMOS_NUN_
# define ATMOS_NUN_ 3
#endif

// The unknowns we have in this model 
#ifndef ATMOS_TT_
# define ATMOS_TT_  1 // Atmospheric temperature
# define ATMOS_QQ_  2 // Vapor anomaly
# define ATMOS_AA_  3 // Albedo

// Auxiliary with a dedicated case in find_row
# define ATMOS_PP_  (ATMOS_NUN_ + 1) // Global precipitation

#endif

// Difference Celsius and Kelvin
#ifndef C2K_
# define C2K_ 273.15
#endif

// Input is expected 1-based
// Output is 1-based
#ifndef FIND_ROW_ATMOS1
# define FIND_ROW_ATMOS1(nun_, n, m, l, i, j, k, XX)				\
	((nun_) * ((k-1)*(n)*(m)+ (n)*(j-1) + (i-1)) + (XX))
#endif

// Input is expected 0-based (except XX and nun)
// Output is 0-based
#ifndef FIND_ROW_ATMOS0
# define FIND_ROW_ATMOS0(nun_, n, m, l, i, j, k, XX)				\
	((nun_) * ((k)*(n)*(m)+ (n)*(j) + (i)) + (XX) - 1)
#endif

#endif
