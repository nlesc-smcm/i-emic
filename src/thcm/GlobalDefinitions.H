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

#ifndef INFO
#  define INFO(s) std::cout << s << std::endl; 
#endif

#ifndef ERROR
#  define ERROR(x,y,z) std::cerr << x << " " << y << " " << z << std::endl;
#endif

#ifndef DEBUG
#  ifdef DEBUGGING
#    define DEBUG(s) std::cout << s << std::endl;
     // The use of # stringifies variable s:
#    define DEBVAR(s) std::cout << #s << ": " << s << std::endl;
#  else
#    define DEBUG(s)
#    define DEBVAR(s)
#  endif
#endif

// macros for checking Epetra calls
#ifndef CHECK_ZERO
#define CHECK_ZERO(funcall) {int ierr = funcall;\
if (ierr) {INFO("Trilinos Error " << ierr << " returned from call " << #funcall);}}
#define CHECK_TRUE(funcall) {bool berr = funcall;\
if (!berr) {INFO("Trilinos call "<<#funcall<<" returned false");}}
#define CHECK_NONNEG(funcall) {int ierr = funcall;\
if (ierr<0) {INFO("Trilinos Error "<<ierr<<" returned from call "<<#funcall);}}
#endif

#ifndef MOD
//! our own 'modulo' function, which behaves like mod in matlab.
//! the C++ built-in '%' operator returns -1%n=-1 and is therefore
//! not very useful for periodic boundaries...
#  define MOD(x,y) (((double)(y)==0.0)? (double)(x): ((double)(x) - floor((double)(x)/((double)(y)))*((double)(y))))
#endif

#endif
