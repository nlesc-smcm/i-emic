/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef THCMDEFS_H
#define THCMDEFS_H

#include <iostream>
#include <cmath>
#include "my_f2c.H"
#include "GlobalDefinitions.H"

// some global definitions for the C++ part of this project
// They should be consistent with usr.F90 (formerly usr.com)
#define PI_ 3.14159265358979323846
// ordering of variables (should be consistent with fortran 
// function find_row2 in matetc.f, beware, however, that    
// this one is 0-based! <except for the XX-part> )                                    
#define FIND_ROW2(nun_,n,m,l,i,j,k,XX)				\
	((nun_)*((k)*(n)*(m)+ (n)*(j) + (i))+(XX)-1)

// number of parameters currently used in THCM (excluding our own params 
// "Time" and "Exponent", which shouldn't be passed to THCM directly).   
#define _NPAR_ 30
// number of unknowns
#define _NUN_ 6
// number of grid neighbours
#define _NP_ 27
//
#ifndef UU
#  define UU 1
#  define VV 2
#  define WW 3
#  define PP 4
#  define TT 5
#  define SS 6
#endif

#endif
