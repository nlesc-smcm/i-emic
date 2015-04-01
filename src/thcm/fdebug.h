#if 0
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
#endif

#ifdef DEBUGGING
#  define _DEBUG_(s) write(*,*)      "(fdb)......", s
#  define _DEBUG2_(s1,s2) write(*,*) " (f#)......", s1, s2
#else
#  define _DEBUG_(s)
#  define _DEBUG2_(s1,s2)
#endif

