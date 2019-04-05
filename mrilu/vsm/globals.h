#ifndef _GLOBALS_H
#define _GLOBALS_H

#define VERSION "VSM Version 1.6.2"
#define DATE    "May 12, 2009"
#define RUG     "RuG/Math.Dep."

#ifdef __hpux
#ifndef _HPUX_SOURCE
#define _HPUX_SOURCE
#endif
#endif

#ifdef __STDC__
#define __P(protos)     protos          /* ANSI C prototypes */
#define _PROTOTYPES
#else
#define __P(protos)     ()              /* K&R C preprocessor */
#endif

/* to ensure to correct loading of XK_Return and XK_BackSpace in keysymdef.h: */
#define XK_MISCELLANY
#include <X11/keysymdef.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "inttypes.h"


#define MAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define BETWEEN(x,min,max) ( (min) < (x) && (x) < (max) )

#define COLORS 21

typedef  uint32  Color;

void LXinit __P(( Display*, Drawable, int ));
void LXflushAll __P(());
void LXdrawPoint __P(( short, short, Color ));
void LXfillRectangle __P(( short, short, unsigned short, unsigned short, Color ));

/***

  ( 5.6 0.0 0.0 ) nRows = 3
  ( 1.4 7.3 2.3 ) nVals = 6
  ( 0.1 0.0 1.0 )

row:  0   1   4   6-----------
      |   |   |_______        |
      |   |           |       |
      v   v           v       v
col:  0   0   1   2   0   2
      |   |   |   |   |   |
val: 5.6 1.4 7.3 2.3 0.1 1.0

row_i == col[row[i]..row[i+1])
len(val) == len(col) == nVals
row[nRows] == nVals
***/

typedef struct {
  /* structure of the matrix */
  int32   nRows, nRowsOrig, nVals;
  int32   *Col, *Row;
  double  *Val;
  double  AbsMin, AbsMax;
  /* the file it came from */
  struct {
    char *name; /* name of the file containing this sparse matrix */
    int  kind;  /* method: ascii, ascii-single, binary or fortran */
  } file ;
  /* a view on this matrix MinRange .. MaxRange */
  struct {
    double min, max;
  } view;
  /* simplified structure */
  int32  *Color;
} SparseMatrix;

#define METHOD_ASCII 0
#define METHOD_BINARY 1
#define METHOD_FORTRAN 2
#define METHOD_ASCII_SINGLE 3


void ReadSparseMatrix __P(( SparseMatrix*, char*, int ));
void RemoveSparseMatrix __P(( SparseMatrix* ));

void VisualizeSparseMatrix __P(( SparseMatrix*, int ));

#define NEW(n,type) (type*) malloc( (n)*sizeof(type) )

#endif /* _GLOBALS_H */

