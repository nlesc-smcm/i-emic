/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006-2010.            *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
/**********************************************************************
 * Modified by T.E. Mulder, Utrecht University 2014/15                *
 * contact: t.e.mulder@uu.nl                                          *
 **********************************************************************/
#ifndef TRIOS_MACROS_H
#define TRIOS_MACROS_H

#include "GlobalDefinitions.H"

// ordering of variables (should be consistent with trilinos_thcm's 

#ifndef FIND_ROW2
# define FIND_ROW2(nun_,n,m,l,i,j,k,XX) ((nun_)*((k)*(n)*(m)+ (n)*(j) + (i))+(XX)-1)
#endif

#ifndef UU
# define UU 1
#endif
#ifndef VV
# define VV 2
#endif
#ifndef WW
# define WW 3
#endif
#ifndef PP
# define PP 4
#endif
#ifndef TT
# define TT 5
#endif
#ifndef SS
# define SS 6
#endif

#endif
