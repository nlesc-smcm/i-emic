/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef MY_F2C_H
#define MY_F2C_H

// define how fortran functions are called
#if defined(GNU)
// g++ and gfortran
//   ##: merging,
//   Append _ to function name
#define FNAME(fun) fun ## _ 
//   Prepend method name with __<module>_MOD_ 
#define F90NAME(mod,fun) __ ## mod ## _MOD_ ## fun  

#elif defined(G95)
// g++ and g95
#define FNAME(fun) fun ## _
#define F90NAME(mod,fun) mod ## _MP_ ## fun 
#elif defined(INTEL)||defined(ASTER)||defined(ORBIT)||defined(HOPF)
// icc and ifort v9.0
#define FNAME(fun) fun ## _
#define F90NAME(mod,fun) mod ## _mp_ ## fun ## _
#elif defined(PGI)
// this is for the pgi compilers pgcc and pgf90
#define FNAME(fun) fun ## _
#define F90NAME(mod,fun) mod ## _ ## fun ## _
#elif defined(HUYGENS)
#define FNAME(fun) fun
#define F90NAME(mod,fun) __ ## mod ## _NMOD_ ## fun
#else
#error no compiler identifier present, cannot convert Fortran names!
#endif

// not sure if this has to be different for some compilers
#define F77CHAR(s) &(s)

#define _SUBROUTINE_(fun) void FNAME(fun)
#define _MODULE_SUBROUTINE_(mod,fun) void F90NAME(mod,fun)

// some global definitions for the C++ part of this project
// They should be consistent with usr.F90 (formerly usr.com)

#endif
