/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef GLOBDEFS_H
#define GLOBDEFS_H

#include <iostream>
#include <cmath>
#include "Teuchos_RCP.hpp"

#ifdef HAVE_MPI
#  include <mpi.h>
#endif

#include "my_f2c.H"

// some global definitions for the C++ part of this project
// They should be consistent with usr.F90 (formerly usr.com)

#define PI_ 3.14159265358979323846

// ordering of variables (should be consistent with fortran 
// function find_row2 in matetc.f, beware, however, that    
// this one is 0-based!)                                    
#define FIND_ROW2(nun_,n,m,l,i,j,k,XX) ((nun_)*((k)*(n)*(m)+ (n)*(j) + (i))+(XX)-1)

// number of parameters currently used in THCM (excluding our own params 
// "Time" and "Exponent", which shouldn't be passed to THCM directly).   
#define _NPAR_ 30
// we added time (0, doesn't count) and some more parameters.
// these are not passed on to THCM, though.
#define _NPAR_TRILI 6

// number of unknowns
#define _NUN_ 6

// number of grid neighbours
#define _NP_ 27

#ifndef UU
#define UU 1
#define VV 2
#define WW 3
#define PP 4
#define TT 5
#define SS 6
#endif

#ifndef INFO
#define INFO(s) std::cout << s << std::endl; 
#endif

#ifndef DEBUG
#  define DEBUG(s) std::cout << s << std::endl;
// The use of # stringifies variable s:
#  define DEBVAR(s) std::cout << #s << std::endl;
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
