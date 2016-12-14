#include "globals.h"

#define MAXGRID 32
#define MAXPOINT 12

/* #define FONT "lucidasanstypewriter-bold-14" */
/*
courB08.pcf.Z -adobe-courier-bold-r-normal--8-80-75-75-m-50-iso8859-1
courB18.pcf.Z -adobe-courier-bold-r-normal--18-180-75-75-m-110-iso8859-1
courB12.pcf.Z -adobe-courier-bold-r-normal--12-120-75-75-m-70-iso8859-1
courB24.pcf.Z -adobe-courier-bold-r-normal--24-240-75-75-m-150-iso8859-1
courB10.pcf.Z -adobe-courier-bold-r-normal--10-100-75-75-m-60-iso8859-1
courB14.pcf.Z -adobe-courier-bold-r-normal--14-140-75-75-m-90-iso8859-1
*/
/*
  #define BIGFONTNAME "-adobe-courier-bold-r-normal--14-140-75-75-m-90-iso8859-1"
  #define SMALLFONTNAME "-adobe-courier-bold-r-normal--12-120-75-75-m-70-iso8859-1"
  #define TINYFONTNAME "-adobe-courier-bold-r-normal--10-100-75-75-m-60-iso8859-1"
*/
#define BIGFONTNAME "9x15"
#define SMALLFONTNAME "8x13"
#define TINYFONTNAME "6x13"

#define STACKSIZE 32

#define UNKNOWN 0
#define GLOBALLOG 1
#define LOCALLOG 2
#define GLOBALSIGN 3

#define QUIT "Quit"
#define GLOBLOG "glob.log"
#define LOCLOG "loc.log"
#define SIGN "sign"
#define PRINT "print"

#define MINCUT ((double) 1.0E-4)
#define MAXCUT ((double) 1.0E4)

#ifdef _PROTOTYPES
void InitColors( SparseMatrix *M, 
		 double MinVal, double MaxVal, int graphtype )
#else
void InitColors( M, MinVal, MaxVal, graphtype )
     SparseMatrix *M;
     double MinVal, MaxVal;
     int graphtype;
#endif
{
  int32  i, color;
  double logmax, logmin, stepsize, dtmp;
  static int    oldgraphtype = UNKNOWN;
  static double oldMinVal, oldMaxVal, oldviewmin, oldviewmax;
  
  if ( graphtype == oldgraphtype )
    if ( MinVal == oldMinVal && MaxVal == oldMaxVal &&
	 M->view.min == oldviewmin && M->view.max == oldviewmax )
      return;
  oldgraphtype = graphtype;
  oldMinVal = MinVal;
  oldMaxVal = MaxVal;
  oldviewmin = M->view.min;
  oldviewmax = M->view.max;
  
  /* Fill M->Color */
  switch (graphtype) {
  case GLOBALLOG:
  case LOCALLOG:
    logmax = log(MaxVal);
    logmin = log(MinVal);
    stepsize = ( logmax - logmin )/COLORS;
    for ( i = 0; i < M->nVals; i++ ) {
      if ( fabs( M->Val[i] ) < M->view.min || fabs( M->Val[i] ) > M->view.max ) 
        M->Color[i] = 0; 
      else { 
        dtmp = ( logmax - log( fabs( M->Val[i] ) ) )/stepsize;
        if ( dtmp < 0.0 ) dtmp = 0.0; else
        if ( dtmp > COLORS ) dtmp = COLORS;
        color = COLORS - (int32)(dtmp);
        if ( color < 1 ) M->Color[i] = 1;
        else if ( color > COLORS ) M->Color[i] = COLORS;
        else M->Color[i] = color;
      }
    }
    break;
  case GLOBALSIGN:
    for ( i = 0; i < M->nVals; i++ ) {
      if ( fabs( M->Val[i] ) < M->view.min || fabs( M->Val[i] ) > M->view.max ) 
        M->Color[i] = 0; else
      if ( M->Val[i] < 0.0 ) M->Color[i] = 1; 
      else if ( M->Val[i] > 0.0 ) M->Color[i] = ( COLORS + 1 )/2;
      else M->Color[i] = 0;
    }
    break;
  }
}

uint32  fg, bg, bd;
int32   row0, col0, rdim, cdim, blocksize;
int xofs, yofs;
int graphtype = GLOBALLOG;
XWindowAttributes xwa;
char posstring[20] = "", oldposstring[20];
char valstring[20] = "", oldvalstring[20];
int printing = 0;  
int baselineskip, leftmargin;
int matrixToggle = False;

XRectangle b[8];

#define SMALLFONT 0
#define BIGFONT 1
#define TINYFONT 2
struct {
  char *name;
  XFontStruct *fontStruct;
} fontlist[3];
int curFont = BIGFONT;

#define FORMAT1 ( curFont == BIGFONT ? "%#10.5f" : "%#10.5f" )
#define FORMAT2 ( curFont == BIGFONT ? "%#10.4e" : "%#10.4e" )
#define FORMAT3 ( curFont == BIGFONT ? "%#14.8f" : "%#14.8f" )
#define FORMAT4 ( curFont == BIGFONT ? "%#14.7e" : "%#14.7e" )

#ifdef _PROTOTYPES
void InputReal( Display *dpy, Window win, GC gc, 
		int x0, int y0, int x1, int y1, 
		double *Val, double Def )
#else
void InputReal(dpy, win, gc, x0, y0, x1, y1, Val, Def)
  Display *dpy;
  Window win;
  GC gc;
  int x0, y0, x1, y1;
  double *Val, Def;
#endif
{ 
  XEvent event;
  char s[20];
  int cp, key, x, y;

  while (XCheckTypedEvent(dpy, ExposureMask|ButtonPressMask|
    ButtonReleaseMask|PointerMotionMask|KeyPressMask, &event));
  XSetForeground(dpy, gc, fg); XSetFunction(dpy, gc, GXcopy);
  if BETWEEN( *Val, MINCUT, MAXCUT ) sprintf(s, FORMAT1, *Val);
  else sprintf(s, FORMAT2, *Val);
  s[10] = 0; 
  cp = strlen(s);
  do {
    XDrawString(dpy, win, gc, x0 + 2, y0+7*baselineskip/6, s, strlen(s));
    x = x0 + XTextWidth(fontlist[curFont].fontStruct, s, cp) + 2;
    y = y0 + y1; 
    XSetFunction(dpy, gc, GXxor); XSetForeground(dpy, gc, fg^bg);
    XDrawLine(dpy, win, gc, x, (y + baselineskip)/2,
      x, (y-baselineskip)/2);
    XSetFunction(dpy, gc, GXcopy); XSetForeground(dpy, gc, fg);
    while (XPeekEvent(dpy, &event), 
      event.type == MotionNotify || event.type == ButtonRelease)
      XNextEvent(dpy, &event);
    if (event.type == KeyPress) {
      XNextEvent(dpy, &event);
      XSetFunction(dpy, gc, GXxor); XSetForeground(dpy, gc, fg^bg);
      XDrawLine(dpy, win, gc, x, (y + baselineskip)/2,
        x, (y-baselineskip)/2);
      XSetFunction(dpy, gc, GXcopy); XSetForeground(dpy, gc, fg);
      key = XLookupKeysym(&event.xkey, event.xkey.state);
      if ( ( key >= '0' && key <= '9' ) || 
	   key == 'E' || key == 'e' || 
	   key == '.' || key == '+' || key == '-' || key == ' ') {
        if (cp < 10) { s[cp++] = key; s[cp] = 0; }
      } else if (key == XK_BackSpace && cp > 0) { 
        XSetFunction(dpy, gc, GXxor); XSetForeground(dpy, gc, fg^bg);
        XDrawString(dpy, win, gc, x0 + 2 + XTextWidth(fontlist[curFont].fontStruct, s, cp-1),
                    y0 + 7*baselineskip/6, &s[cp-1], 1);
        XSetFunction(dpy, gc, GXcopy); XSetForeground(dpy, gc, fg);
        s[--cp] = 0;
      }
    }
  } while (event.type == KeyPress && key != XK_Return );
  if (cp == 0) *Val = Def;
  else sscanf(s, "%lg", Val);
}

#ifdef _PROTOTYPES
void GetLocalMinMax( SparseMatrix *M, double *LocMin, double *LocMax )
#else
void GetLocalMinMax( M, LocMin, LocMax )
     SparseMatrix *M;
     double *LocMin, *LocMax;
#endif
{ 
  double Min, Max;
  int First;
  int32  Row, Col;

  Min = 0.0; Max = 0.0;
  for ( Row = row0, First = True; Row < row0 + rdim; Row++ )
    for ( Col = M->Row[Row]; Col < M->Row[Row+1]; Col++ ) {
      if ( M->Col[Col] < col0 ) continue;
      if ( M->Col[Col] >= col0 + cdim ) break;
      if ( fabs( M->Val[Col] ) >= M->view.min && fabs( M->Val[Col] ) <= M->view.max) {
        if (First) { First = False; Min = Max = fabs( M->Val[Col] ); } else
	  if ( fabs( M->Val[Col] ) < Min ) Min = fabs( M->Val[Col] ); else
	    if ( fabs( M->Val[Col] ) > Max ) Max = fabs( M->Val[Col]);
      }
    }
  *LocMin = Min; *LocMax = Max;
}

#ifdef _PROTOTYPES
void MakeLegenda ( Display *dpy, Window win, GC gc,
		   int x0, int y0, int width, int height, /* dimensions east */
		   int x0b, int y0b, int widthb, int heightb, /* dimensions south */
		   SparseMatrix *M,
		   double MinVal, double MaxVal,
		   FILE *f )
#else
void MakeLegenda (dpy, win, gc, x0, y0, width, height, x0b, y0b, widthb, heightb, 
  M, MinVal, MaxVal, f)
Display *dpy;
Window win;
GC gc;
int x0, x0b;
int y0, y0b;
int width, widthb;
int height, heightb;
double MinVal, MaxVal;
SparseMatrix *M;
FILE *f;
#endif
{
  int i, y, barwidth, rightmargin;
  char s[20];
  double dtmp, logmin, logmax;
  int newFont;

  XSetForeground(dpy, gc, fg);

  newFont = BIGFONT;
  baselineskip = fontlist[newFont].fontStruct->max_bounds.ascent + fontlist[newFont].fontStruct->max_bounds.descent;
  while( ( 9 + COLORS )*baselineskip > height - 11*baselineskip/2 - 7*baselineskip/6 )
    baselineskip--;
  if ( baselineskip < fontlist[newFont].fontStruct->max_bounds.ascent ||
       ( widthb - 6*(widthb/31) )/5 < XTextWidth( fontlist[newFont].fontStruct, GLOBLOG, strlen(GLOBLOG) ) ||
       ( width < XTextWidth( fontlist[newFont].fontStruct, VERSION, strlen(VERSION) ) ) ) {
    newFont = SMALLFONT;
    baselineskip = fontlist[newFont].fontStruct->max_bounds.ascent + fontlist[newFont].fontStruct->max_bounds.descent;
    while( ( 9 + COLORS )*baselineskip > height - 11*baselineskip/2 - 7*baselineskip/6 )
      baselineskip--;
    if ( baselineskip < fontlist[newFont].fontStruct->max_bounds.ascent ||
	 ( widthb - 6*(widthb/31) )/5 < XTextWidth( fontlist[newFont].fontStruct, GLOBLOG, strlen(GLOBLOG) ) ||
	 ( width < XTextWidth( fontlist[newFont].fontStruct, VERSION, strlen(VERSION) ) ) ) {
      newFont = TINYFONT;
      baselineskip = fontlist[newFont].fontStruct->max_bounds.ascent + fontlist[newFont].fontStruct->max_bounds.descent;
      while( ( 9 + COLORS )*baselineskip > height - 11*baselineskip/2 - 7*baselineskip/6 )
	baselineskip--;
    }
  }

  if ( newFont != curFont ) {
    curFont = newFont;
    XSetFont(dpy, gc, fontlist[curFont].fontStruct->fid);
  }


  leftmargin = width/20;
  rightmargin = leftmargin/2;
  barwidth = baselineskip*7/5;

  /* First line: version */
  XDrawString( dpy, win, gc, 
	       x0 + ( width - XTextWidth( fontlist[curFont].fontStruct, VERSION, strlen(VERSION) ) )/2,
	       y0 + baselineskip, VERSION, strlen(VERSION) );
  /* Second line: copyright */
  XDrawString( dpy, win, gc, 
	       x0 + ( width - XTextWidth( fontlist[curFont].fontStruct, RUG, strlen(RUG) ) )/2,
	       y0 + 2*baselineskip, RUG, strlen(RUG) );

  /* Last lines: Max and Min */
  b[4].x = b[5].x = x0 + leftmargin + rightmargin +
    MAX( XTextWidth( fontlist[curFont].fontStruct, "Min", strlen( "Min" ) ),
	 XTextWidth( fontlist[curFont].fontStruct, "Min", strlen( "Min" ) ) );
  b[4].width = b[5].width = width - b[4].x + x0 - rightmargin;
  b[4].y = height - (11*baselineskip/2) - (7*baselineskip/6);
  b[4].height = b[5].height = 5*baselineskip/3;
  b[5].y = b[4].y + 2*baselineskip;

  XSetFunction( dpy, gc, GXcopy ); XSetForeground( dpy, gc, fg );
  XDrawString( dpy, win, gc, 
	       x0 + leftmargin, height - 11*baselineskip/2,
	       "Max", strlen( "Max" ) );
  XDrawRectangles(dpy, win, gc, &b[4], 1 );
  if BETWEEN( M->view.max, MINCUT, MAXCUT ) sprintf(s, FORMAT1, M->view.max);
  else sprintf( s, FORMAT2, M->view.max );
  s[10] = 0;
  XDrawString( dpy, win, gc, 
	       b[4].x+2, height - 11*baselineskip/2,
	       s, strlen(s) );

  XDrawString( dpy, win, gc, 
	       x0 + leftmargin, height - 7*baselineskip/2,
	       "Min", strlen( "Min" ) );
  if BETWEEN( M->view.min, MINCUT, MAXCUT ) sprintf(s, FORMAT1, M->view.min);
  else sprintf( s, FORMAT2, M->view.min );
  s[10] = 0;
  XDrawString( dpy, win, gc, 
	       b[5].x+2, height - 7*baselineskip/2,
	       s, strlen(s) );
  XDrawRectangles(dpy, win, gc, &b[5], 1 );

  /* bottom bar */
  b[0].y = b[1].y = b[2].y = b[3].y = b[6].y = y0b + heightb/6;
  b[0].height = b[1].height = b[2].height = b[3].height = b[6].height = heightb - 2*(heightb/6);
  {
    int offset = (widthb/31);
    int width = (widthb - 6*offset)/5;
    b[0].width = b[1].width = b[2].width = b[3].width = b[6].width = width;
    b[0].x = x0b + offset;
    b[1].x = b[0].x + width + offset;
    b[2].x = b[1].x + width + offset;
    b[3].x = b[2].x + width + offset;
    b[6].x = b[3].x + width + offset;
  }
  { 
    XFontStruct *font = fontlist[curFont].fontStruct;
    XDrawString( dpy, win, gc, 
		 b[0].x + ( b[0].width - XTextWidth(font, QUIT, strlen(QUIT) ) )/2, 
		 b[0].y + ( b[0].height + font->max_bounds.ascent - font->max_bounds.descent )/2, 
		 QUIT, strlen(QUIT) );
    XDrawString( dpy, win, gc, 
		 b[1].x + ( b[1].width - XTextWidth(font, GLOBLOG, strlen(GLOBLOG) ) )/2, 
		 b[1].y + ( b[1].height + font->max_bounds.ascent - font->max_bounds.descent )/2, 
		 GLOBLOG, strlen(GLOBLOG) );
    XDrawString( dpy, win, gc, 
		 b[2].x + ( b[2].width - XTextWidth(font, LOCLOG, strlen(LOCLOG) ) )/2, 
		 b[2].y + ( b[2].height + font->max_bounds.ascent - font->max_bounds.descent )/2, 
		 LOCLOG, strlen(LOCLOG) );
    XDrawString( dpy, win, gc, 
		 b[3].x + ( b[3].width - XTextWidth(font, SIGN, strlen(SIGN) ) )/2, 
		 b[3].y + ( b[3].height + font->max_bounds.ascent - font->max_bounds.descent )/2, 
		 SIGN, strlen(SIGN) );
    XDrawString( dpy, win, gc, 
		 b[6].x + ( b[6].width - XTextWidth(font, PRINT, strlen(PRINT) ) )/2, 
		 b[6].y + ( b[6].height + font->max_bounds.ascent - font->max_bounds.descent )/2, 
		 PRINT, strlen(PRINT) );
  }
  XSetFunction( dpy, gc, GXxor ); 
  XSetForeground( dpy, gc, fg^bg );
  XDrawRectangles( dpy, win, gc, &b[0], 1 );
  if ( graphtype == GLOBALLOG ) XFillRectangles( dpy, win, gc, &b[1], 1 );
  else XDrawRectangles( dpy, win, gc, &b[1], 1 );
  if ( graphtype == LOCALLOG ) XFillRectangles( dpy, win, gc, &b[2], 1 );
  else XDrawRectangles( dpy, win, gc, &b[2], 1 );
  if ( graphtype == GLOBALSIGN ) XFillRectangles( dpy, win, gc, &b[3], 1 );
  else XDrawRectangles( dpy, win, gc, &b[3], 1 );
  XDrawRectangles( dpy, win, gc, &b[6], 1 );
  XSetFunction( dpy, gc, GXcopy ); 
  XSetForeground( dpy, gc, fg );

  if ( matrixToggle ) {
    int boxwidth = 11*baselineskip/2;
    int x, y, w, h;
    b[7].x = x0 + ( width - boxwidth )/2; b[7].y = y0 + 5*baselineskip/2;
    b[7].height = b[7].width = boxwidth;
    XDrawRectangles( dpy, win, gc, &b[7], 1 );
    x = col0*( boxwidth - 2 )/( M->nRows - 1);
    y = row0*( boxwidth - 2 )/( M->nRows - 1);
    w = col0 + cdim == M->nRows ? boxwidth - 1 - x: MAX( 1, cdim*( boxwidth - 2 )/( M->nRows - 1) );
    h = row0 + rdim == M->nRows ? boxwidth - 1 - y: MAX( 1, rdim*( boxwidth - 2 )/( M->nRows - 1) );
    LXfillRectangle( b[7].x + x + 1, b[7].y + y + 1, w, h, 1 );
    if ( w == 1 && h == 1 ) {
      int len0 = 2, len = 10;
      XDrawLine( dpy, win, gc, b[7].x + x + 1, b[7].y + y + 1 + len0, b[7].x + x + 1, b[7].y + y + 1 + len );
      XDrawLine( dpy, win, gc, b[7].x + x + 1, b[7].y + y + 1 - len0, b[7].x + x + 1, b[7].y + y + 1 - len );
      XDrawLine( dpy, win, gc, b[7].x + x + 1 + len0, b[7].y + y + 1, b[7].x + x + 1 + len, b[7].y + y + 1 );
      XDrawLine( dpy, win, gc, b[7].x + x + 1 - len0, b[7].y + y + 1, b[7].x + x + 1 - len, b[7].y + y + 1 );
    }
  } else {
    b[7].x = x0 + leftmargin/2; b[7].y = y0 + 5*baselineskip/2;
    b[7].height = 11*baselineskip/2;
    b[7].width = width - leftmargin/2 - rightmargin;
    XDrawRectangles( dpy, win,gc, &b[7], 1 );
    
    y = y0 + 7*baselineskip/2;
    sprintf( s, "Rows: %ld -", row0 + 1 );
    XDrawString( dpy, win, gc, x0 + leftmargin, y, s, strlen(s) ); y += baselineskip;
    sprintf( s, "      %ld", row0 + rdim );
    XDrawString( dpy, win, gc, x0 + leftmargin, y, s, strlen(s) ); y += baselineskip;
    sprintf( s, "Cols: %ld -", col0 + 1 );
    XDrawString( dpy, win, gc, x0 + leftmargin, y, s, strlen(s) ); y += baselineskip;
    sprintf( s, "      %ld", col0 + cdim );
    XDrawString( dpy, win, gc, x0 + leftmargin, y, s, strlen(s) ); y += baselineskip;
    if ( blocksize > 0 ) sprintf( s, "Magn: %ld:1", blocksize );
    else sprintf( s, "Magn: 1:%ld", -blocksize );
    XDrawString( dpy, win, gc, x0 + leftmargin, y, s, strlen(s) ); y += baselineskip;
  }

  y = y0 + 9*baselineskip;
  switch (graphtype) {
  case GLOBALLOG:
  case LOCALLOG:
    for (i = 0; i < COLORS; i++) {
      LXfillRectangle( x0 + leftmargin, y + i*baselineskip - baselineskip/3, 
		       barwidth, baselineskip, COLORS - i ); 
      if (printing) {
	fprintf(f, "1 3 0 1 -1 0 0 21 0.00000 1 0.000 %ld %ld %d %d %ld %ld %ld %ld\n",
		MAXGRID*(x0)/blocksize + 150, 
		MAXGRID*(y0)/blocksize + i*100,
		1+(COLORS - i)*MAXPOINT/COLORS, 1+(COLORS - i)*MAXPOINT/COLORS,
		MAXGRID*(x0)/blocksize + 150, 
		MAXGRID*(y0)/blocksize + i*100,
		MAXGRID*(x0)/blocksize + 150 +1+(COLORS - i)*MAXPOINT/COLORS,
		MAXGRID*(y0)/blocksize + i*100 +1+(COLORS - i)*MAXPOINT/COLORS);
      }
    }
    XSetForeground( dpy, gc, fg );
    for ( i = 0; i <= COLORS; i++ ) {     
      XDrawLine( dpy, win, gc, 
        x0 + leftmargin + barwidth, y + i*baselineskip - baselineskip/3,
        (x0 + leftmargin + barwidth) + 2*leftmargin, y + i*baselineskip - baselineskip/3 );
      logmin = log(MinVal);
      logmax = log(MaxVal);
      dtmp = exp( logmin + (COLORS - i)*(logmax - logmin)/COLORS );
      if BETWEEN( dtmp, MINCUT, MAXCUT ) sprintf(s, FORMAT1, dtmp);
      else sprintf( s, FORMAT2, dtmp );
      s[10] = 0;
      XDrawString( dpy, win, gc, 
		   x0 + leftmargin + barwidth + 2*leftmargin, y + i*baselineskip, 
		   s, strlen(s));
      if (printing) {
        fprintf(f, "4 0 0 72 0 -1 0 0.00000 4 23 96 %ld %ld %s\001\n",
		MAXGRID*(x0)/blocksize + 250,
		MAXGRID*(y0)/blocksize + i*100 - 26,
		s);
      }
    }
    break;
  case GLOBALSIGN:
    i = 0;
    LXfillRectangle( x0 + leftmargin, y + i*baselineskip - baselineskip/3,
		     barwidth, baselineskip, COLORS/2 + 1);
    i = 1;
    LXfillRectangle( x0 + leftmargin, y + i*baselineskip - baselineskip/3,
		     barwidth, baselineskip, 1);

    XSetForeground(dpy, gc, fg);
    i = 0;
    XDrawString( dpy, win, gc, x0 + leftmargin + barwidth + 2*leftmargin, 
      y + i*baselineskip + baselineskip/2, " > 0", 4);
    i = 1;
    XDrawString( dpy, win, gc, x0 + leftmargin + barwidth + 2*leftmargin, 
      y + i*baselineskip + baselineskip/2, " < 0", 4);

    break;
  }
}

#ifdef _PROTOTYPES
void DoManyToOne ( Display *dpy, Window win, GC gc,
		   SparseMatrix *M,
		   int32 row0, int32 col0, int32 rdim, int32 cdim,
		   int32 blocksize,
		   int xofs, int yofs )
#else
void DoManyToOne( dpy, win, gc, M, row0, col0, rdim, cdim, blocksize, xofs, yofs )
     Display *dpy;
     Window win;
     GC gc;
     SparseMatrix *M;
     int32  row0;
     int32  col0;
     int32  rdim;
     int32  cdim;
     int32  blocksize;
     int xofs;
     int yofs;
#endif
{ 
  int32  ColIndex[1024], ColLast[1024];
  int    i, isequal, first, block;
  int32  rownr1, colmin = 0, maxcolor;
  double maxvalue;
  
  for (rownr1 = row0; rownr1 < row0+rdim; rownr1 += blocksize) {
    for (i = 0; i < blocksize && rownr1+i < row0+rdim; i++) {
      ColIndex[i] = M->Row[rownr1+i];
      ColLast[i] = M->Row[rownr1+i+1];
      while (ColIndex[i] != ColLast[i] && (M->Col)[ColIndex[i]] < col0) ColIndex[i]++;
      while (ColIndex[i] != ColLast[i] && (M->Col)[ColLast[i]-1] >= col0+cdim) ColLast[i]--;
    } 
    while (1) {
      for (i = 0, isequal = True; 
           i < blocksize && rownr1+i < row0+rdim && isequal; i++) {
        isequal = (ColIndex[i] == ColLast[i]);
      }
      if (isequal) break;
      for (i = 0, first = True; i < blocksize && rownr1+i < row0+rdim; i++)
        if (ColIndex[i] != ColLast[i]) {
          if (first) {
            first = False; colmin = (M->Col)[ColIndex[i]];
	  } else if (colmin > (M->Col)[ColIndex[i]]) {
            colmin = (M->Col)[ColIndex[i]];
	  }
	}
      /* colmin is now first new column number in rows */
      block = (colmin - col0)/blocksize;
      switch (graphtype) {
      case GLOBALLOG:
      case LOCALLOG:
	maxcolor = 0;
        for (i = 0, first = True; i < blocksize && rownr1+i < row0+rdim; i++) {
          while (ColIndex[i] != ColLast[i] && 
                 ((M->Col)[ColIndex[i]] - col0)/blocksize == block) {
	    if ( (M->Color)[ColIndex[i]] != 0 ) {
	      if (first) { first = False; maxcolor = (M->Color)[ColIndex[i]]; }
	      else if (maxcolor < (M->Color)[ColIndex[i]]) 
		maxcolor = (M->Color)[ColIndex[i]];
	    }
            ColIndex[i]++;
          }
        }
        LXdrawPoint(xofs + block, yofs + (rownr1 - row0)/blocksize, maxcolor);
        break;
      case GLOBALSIGN:
	maxvalue = 0;
        for (i = 0, first = True; i < blocksize && rownr1+i < row0+rdim; i++) {
          while (ColIndex[i] != ColLast[i] && 
                 ((M->Col)[ColIndex[i]] - col0)/blocksize == block) {
            if (fabs((M->Val)[ColIndex[i]]) >= M->view.min && 
                fabs((M->Val)[ColIndex[i]]) <= M->view.max) {
              if (first) { first = False; maxvalue = (M->Val)[ColIndex[i]]; }
              else if (fabs(maxvalue) < fabs((M->Val)[ColIndex[i]]))
		maxvalue = (M->Val)[ColIndex[i]];
            }
            ColIndex[i]++;
          }
        }
        LXdrawPoint(xofs + block, yofs + (rownr1 - row0)/blocksize, 
		    first ? 0 : (maxvalue < 0.0 ? 1 : COLORS/2+1));
        break;
      }
    }
  }
  LXflushAll();
}

#ifdef _PROTOTYPES
void DoOneToMany( Display *dpy, Window win, GC gc,
		  SparseMatrix *M,
		  int32 row0, int32 col0, int32 rdim, int32 cdim,
		  int32 blocksize,
		  int xofs, int yofs,
		  FILE *f )
#else
void DoOneToMany (dpy, win, gc, M, row0, col0, rdim, cdim, blocksize, xofs, yofs, f)
     Display *dpy;
     Window win;
     GC gc;
     SparseMatrix *M;
     int32  row0;
     int32  col0;
     int32  rdim;
     int32   cdim;
     int32  blocksize;
     int xofs;
     int yofs;
     FILE *f;
#endif
{
  int32  colnr1, rownr1;

  for (rownr1 = row0; rownr1 < row0+rdim; rownr1++) {
    for (colnr1 = (M->Row)[rownr1]; colnr1 < (M->Row)[rownr1+1]; colnr1++) {
      if ((M->Col)[colnr1] < col0) continue;
      if ((M->Col)[colnr1] >= col0+cdim) break;
      if ( blocksize == 1 ) {
        LXdrawPoint(xofs + blocksize*((M->Col)[colnr1]-col0),  
		    yofs + blocksize*(rownr1-row0), (M->Color)[colnr1]);
	if (printing) {  
	  fprintf(f, "1 3 0 1 -1 0 0 21 0.00000 1 0.000 %ld %ld %ld %ld %ld %ld %ld %ld\n",
		 MAXGRID*(xofs + blocksize*((M->Col)[colnr1]-col0)), 
		 MAXGRID*(yofs + blocksize*(rownr1-row0)),
		 1+(M->Color)[colnr1]*MAXPOINT/COLORS, 1+(M->Color)[colnr1]*MAXPOINT/COLORS,
		 MAXGRID*(xofs + blocksize*((M->Col)[colnr1]-col0)), 
		 MAXGRID*(yofs + blocksize*(rownr1-row0)),
		 MAXGRID*(xofs + blocksize*((M->Col)[colnr1]-col0))+1+(M->Color)[colnr1]*MAXPOINT/COLORS, 
		 MAXGRID*(yofs + blocksize*(rownr1-row0)));
	}
      } else {
        LXfillRectangle(xofs + blocksize*((M->Col)[colnr1]-col0), 
			yofs + blocksize*(rownr1-row0), 
			blocksize, blocksize, (M->Color)[colnr1]);
	if (printing) {
	  fprintf(f, "1 3 0 1 -1 0 0 21 0.00000 1 0.000 %ld %ld %ld %ld %ld %ld %ld %ld\n",
		 MAXGRID*(xofs/blocksize + ((M->Col)[colnr1]-col0)), 
		 MAXGRID*(yofs/blocksize + (rownr1-row0)),
		 1+(M->Color)[colnr1]*MAXPOINT/COLORS, 1+(M->Color)[colnr1]*MAXPOINT/COLORS,
		 MAXGRID*(xofs/blocksize + ((M->Col)[colnr1]-col0)), 
		 MAXGRID*(yofs/blocksize + (rownr1-row0)),
		 MAXGRID*(xofs/blocksize + ((M->Col)[colnr1]-col0))+1+(M->Color)[colnr1]*MAXPOINT/COLORS, 
		 MAXGRID*(yofs/blocksize + (rownr1-row0)));
	}
      }
    }
  }
  LXflushAll();
}

#ifdef _PROTOTYPES
void WriteMatrix ( Display *dpy, Window win, GC gc, SparseMatrix *M, FILE *f )
/* Visualize in the already opened window the matrix M */
#else
void WriteMatrix (dpy, win, gc, M, f )
     Display *dpy;
     Window win;
     GC gc;
     SparseMatrix *M;
     FILE *f;
/* Visualize in the already opened window the matrix M */
#endif
{
  int labelwidth, labelheight;
  double LocMin, LocMax;

  if (printing) { fprintf(f, "#FIG 2.1\n80 2\n"); }
  XClearWindow(dpy, win); posstring[0] = 0; valstring[0] = 0;
  XGetWindowAttributes(dpy, win, &xwa);

  labelwidth = MIN( xwa.width/5, XTextWidth( fontlist[BIGFONT].fontStruct, VERSION, strlen(VERSION) ) );
  labelheight = xwa.height/20;

  /* Normalize showwindow over the matrix */
  if (rdim > M->nRows) rdim = M->nRows;
  if (cdim > M->nRows) cdim = M->nRows;
  if (col0 + cdim > M->nRows) col0 = M->nRows - cdim;
  if (row0 + rdim > M->nRows) row0 = M->nRows - rdim;

  /* Determine blocksize */
  if (rdim*(xwa.width-1-labelwidth-2) < cdim*(xwa.height-1-labelheight-2)) {
    /* Determine blocksize from cdim and xwa.width-labelwidth-2 */
    if (cdim <= xwa.width-1-labelwidth-2) blocksize = (xwa.width-1-labelwidth-2)/cdim;
    else blocksize = -((cdim-1)/(xwa.width-1-labelwidth-2)+1);
  } else {
    /* Determine blocksize from rdim and xwa.height-labelheight-2 */
    if (rdim <= xwa.height-1-labelheight-2) blocksize = (xwa.height-1-labelheight-2)/rdim;
    else blocksize = -((rdim-1)/(xwa.height-1-labelheight-2)+1);
  }

  /* Adjust showwindow sizes */
  if (blocksize > 0) {
    rdim = (xwa.height-1-labelheight-2)/blocksize; 
    cdim = (xwa.width-1-labelwidth-2)/blocksize;    
  } else {
    rdim = (xwa.height-1-labelheight-2)*(-blocksize); 
    cdim = (xwa.width-1-labelwidth-2)*(-blocksize);
  }
  /* printf("BLOCKSIZE=%d\n", blocksize); */
  /* Normalize showwindow over the matrix */
  if (rdim > M->nRows) rdim = M->nRows;
  if (cdim > M->nRows) cdim = M->nRows;
  if (col0 + cdim > M->nRows) col0 = M->nRows - cdim;
  if (row0 + rdim > M->nRows) row0 = M->nRows - rdim;

  if (blocksize > 0) {
    xofs = (xwa.width-1-labelwidth - cdim*blocksize)/2;
    yofs = (xwa.height-1-labelheight - rdim*blocksize)/2;

    LXflushAll();
    LXfillRectangle(xofs - 1, yofs - 1,
		  cdim*blocksize+2, rdim*blocksize+2, 0);
    LXflushAll();

    XSetForeground(dpy, gc, fg);
    XDrawRectangle(dpy, win, gc, xofs - 1, yofs - 1, 
      cdim*blocksize+2, rdim*blocksize+2);
  } else {
    xofs = (xwa.width-1-labelwidth - cdim/(-blocksize))/2;
    yofs = (xwa.height-1-labelheight - rdim/(-blocksize))/2;

    LXflushAll();
    LXfillRectangle(xofs - 1, yofs - 1,
		  cdim/(-blocksize)+2, rdim/(-blocksize)+2, 0);
    LXflushAll();

    XSetForeground(dpy, gc, fg);
    XDrawRectangle(dpy, win, gc, xofs - 1, yofs - 1, 
      cdim/(-blocksize)+2, rdim/(-blocksize)+2);
  }

  /* Make picture */
  if (graphtype == LOCALLOG) GetLocalMinMax(M, &LocMin, &LocMax );
  else { LocMin = M->AbsMin; LocMax = M->AbsMax; }
  InitColors(M, LocMin, LocMax, graphtype );
  MakeLegenda(dpy, win, gc, 
	      xwa.width - labelwidth, 0, labelwidth, xwa.height, 
	      0, xwa.height - labelheight, xwa.width - labelwidth, labelheight, M, LocMin, LocMax, f);
  
  if (blocksize < 0) {
    DoManyToOne(dpy, win, gc, M, row0, col0, rdim, cdim, -blocksize, xofs, yofs);
  } else {
    DoOneToMany(dpy, win, gc, M, row0, col0, rdim, cdim, blocksize, xofs, yofs, f);
  } 
}


  XWMHints xwmh = {
    (InputHint|StateHint),
    True,
    NormalState,
    0,
    0,
    0, 0,
    0,
    0
    };

/* opens a window and shows the matrix M in it. */
#ifdef _PROTOTYPES
void VisualizeSparseMatrix( 
     SparseMatrix *M,
     int reverse )
#else
void VisualizeSparseMatrix( M, reverse )
     SparseMatrix *M;
     int reverse;
#endif
{ 
  typedef struct {
    int32  row0, col0, rdim, cdim, blocksize;
    int xofs, yofs;
  } stackelem;

  stackelem stack[STACKSIZE];
  int sp = 0;
  
  Display *dpy;
  Window win;
  GC gc;
  XEvent event;
  XSizeHints xsh;
  XSetWindowAttributes xswa;
  unsigned int bw = 1;
  int i, CreatingZoomWindow, DoExpose;
  int Col, zx1 = 0, zy1 = 0, zx2 = 0, zy2 = 0, ztx, zty, curc, curr;
  double curv;
  
  char name[80], icon_name[20];

  FILE *file = NULL;
  char filen[256];

  sprintf(name, "Visualize Sparse Matrix - %s", M->file.name);
  sprintf(icon_name, "VSM - %s", M->file.name);
  /* Initialize a window */
  if ((dpy = XOpenDisplay(NULL)) == NULL) {
    fprintf(stderr, "Can't open %s.\n", XDisplayName(NULL));
    exit(1);
  }
  fontlist[BIGFONT].name = BIGFONTNAME;
  if ((fontlist[BIGFONT].fontStruct = XLoadQueryFont(dpy, BIGFONTNAME)) == NULL) {
    fprintf(stderr, "Display %s doesn't know font %s.\n", DisplayString(dpy), BIGFONTNAME);
    exit(1);
  }
  fontlist[SMALLFONT].name = SMALLFONTNAME;
  if ((fontlist[SMALLFONT].fontStruct = XLoadQueryFont(dpy, SMALLFONTNAME)) == NULL) {
    fontlist[SMALLFONT] = fontlist[BIGFONT];
    fprintf( stderr, "Unknown font %s.\nSubstituted by font %s\n", 
	     SMALLFONTNAME, fontlist[SMALLFONT].name );
  }
  fontlist[TINYFONT].name = TINYFONTNAME;
  if ((fontlist[TINYFONT].fontStruct = XLoadQueryFont(dpy, TINYFONTNAME)) == NULL) {
    fontlist[TINYFONT] = fontlist[SMALLFONT];
    fprintf( stderr, "Unknown font %s.\nSubstituted by font %s\n", 
	     TINYFONTNAME, fontlist[TINYFONT].name );
  }
  bd = WhitePixel(dpy, DefaultScreen(dpy));
  fg = BlackPixel(dpy, DefaultScreen(dpy));
  bg = WhitePixel(dpy, DefaultScreen(dpy));
  xsh.flags = (PPosition|PSize|PMinSize);
  xsh.width = 640+84;
  xsh.height = 480+116;
  {
    int baselineskip = fontlist[TINYFONT].fontStruct->max_bounds.ascent - 1;
    xsh.min_height = ( 9 + COLORS )*baselineskip + 11*baselineskip/2 + 7*baselineskip/6 + baselineskip/2;
    /* xsh.min_width = (xsh.min_height*19)/12; */
    xsh.min_width = (xsh.min_height*19)/16;
  }
  xsh.x = (DisplayWidth(dpy, DefaultScreen(dpy)) - xsh.width)/2;
  xsh.y = (DisplayWidth(dpy, DefaultScreen(dpy)) - xsh.height)/2;

  win = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy),
    xsh.x, xsh.y, xsh.width, xsh.height, bw, bd, bg);
  XSetStandardProperties(dpy, win, name, icon_name, None, NULL, 0, &xsh);
  XSetWMHints(dpy, win, &xwmh);

  xswa.colormap = DefaultColormap(dpy, DefaultScreen(dpy));
  xswa.bit_gravity = CenterGravity;
  XChangeWindowAttributes(dpy, win, (CWColormap|CWBitGravity), &xswa);

  gc = DefaultGC(dpy, DefaultScreen(dpy));
  XSetFont(dpy, gc, fontlist[curFont].fontStruct->fid);
  XSetForeground(dpy, gc, fg);
  XSetBackground(dpy, gc, bg);

  LXinit(dpy, win, reverse);

  XSelectInput(dpy, win, 
    ExposureMask|/*StructureNotifyMask|*/ButtonPressMask|
    ButtonReleaseMask|PointerMotionMask|KeyPressMask);

  XMapWindow(dpy, win);

  row0 = 0; col0 = 0; rdim = M->nRows; cdim = M->nRows; 
  CreatingZoomWindow = False; DoExpose = True;
  M->view.min = M->AbsMin; M->view.max = M->AbsMax;

  /* Eventloop with writing the window */
  while (1) {
    XNextEvent(dpy, &event);
    if ((event.type == ConfigureNotify) || (event.type == Expose)) {
      while (XCheckTypedEvent(dpy, ConfigureNotify, &event) ||
             XCheckTypedEvent(dpy, Expose, &event));
      WriteMatrix( dpy, win, gc, M, file ); 
    } else if (event.type == ButtonPress && !CreatingZoomWindow) {
      if (event.xbutton.button == 1) {
        if (event.xbutton.x >= b[1].x && event.xbutton.x <= b[1].width + b[1].x - 1 &&
            event.xbutton.y >= b[1].y && event.xbutton.y <= b[1].height + b[1].y - 1) {
          graphtype = GLOBALLOG; WriteMatrix(dpy, win, gc, M, file);
        } else
        if (event.xbutton.x >= b[2].x && event.xbutton.x <= b[2].width + b[2].x - 1 &&
            event.xbutton.y >= b[2].y && event.xbutton.y <= b[2].height + b[2].y - 1) {
          graphtype = LOCALLOG; WriteMatrix(dpy, win, gc, M, file);
        } else 
        if (event.xbutton.x >= b[3].x && event.xbutton.x <= b[3].width + b[3].x - 1 &&
            event.xbutton.y >= b[3].y && event.xbutton.y <= b[3].height + b[3].y - 1) {
          graphtype = GLOBALSIGN; WriteMatrix(dpy, win, gc, M, file);
        } else
        if (event.xbutton.x >= b[7].x && event.xbutton.x <= b[7].width + b[7].x - 1 &&
            event.xbutton.y >= b[7].y && event.xbutton.y <= b[7].height + b[7].y - 1) {
          matrixToggle = !matrixToggle; WriteMatrix(dpy, win, gc, M, file);
        } else
        if (event.xbutton.x >= b[4].x && event.xbutton.x <= b[4].width + b[4].x - 1 &&
            event.xbutton.y >= b[4].y && event.xbutton.y <= b[4].height + b[4].y - 1) {
          InputReal(dpy, win, gc, b[4].x, b[4].y, b[4].width + b[4].x - 1, b[4].height + b[4].y - 1, &M->view.max, M->AbsMax); 
          if (M->view.min > M->view.max) M->view.max = M->view.min;
	  WriteMatrix(dpy, win, gc, M, file);
        } else
        if (event.xbutton.x >= b[5].x && event.xbutton.x <= b[5].width + b[5].x - 1 &&
            event.xbutton.y >= b[5].y && event.xbutton.y <= b[5].height + b[5].y - 1) {
          InputReal(dpy, win, gc, b[5].x, b[5].y, b[5].width + b[5].x - 1, b[5].height + b[5].y - 1, &M->view.min, M->AbsMin); 
          if (M->view.min > M->view.max) M->view.min = M->view.max;
	  WriteMatrix(dpy, win, gc, M, file);
        } else 
        if (event.xbutton.x >= b[6].x && event.xbutton.x <= b[6].width + b[6].x - 1 &&
            event.xbutton.y >= b[6].y && event.xbutton.y <= b[6].height + b[6].y - 1) {
          printing = 1;
	  { int i = 0;
	    while (i <= 999 && 
		   (sprintf(filen, "vsm-%s-%.3d.fig", M->file.name, i), 
		    access(filen, R_OK) == 0)) i++;
	    if (i > 999) { printf("Too many files?\n"); exit(1); }
	  }
	  if ((file = fopen(filen, "w")) == NULL) {
	    fprintf(stderr, "Can't open %s for writing.\n", filen); exit(1); }
	  WriteMatrix(dpy, win, gc, M, file);
	  fclose(file);
	  printing = 0;
        } else {
          CreatingZoomWindow = True; 
	  zx1 = zx2 = event.xbutton.x; zy1 = zy2 = event.xbutton.y;
        }
      } else if (event.xbutton.button == 3) {
        if (sp > 0) {
          sp--; row0 = stack[sp].row0; col0 = stack[sp].col0; rdim = stack[sp].rdim;
          cdim = stack[sp].cdim; blocksize = stack[sp].blocksize;
          xofs = stack[sp].xofs; yofs = stack[sp].yofs;
          WriteMatrix(dpy, win, gc, M, file);
        }
      }
    } else if (event.type == ButtonRelease) {
      if (event.xbutton.button == 1) {
        if (event.xbutton.x >= b[0].x && event.xbutton.x <= b[0].width + b[0].x - 1 &&
            event.xbutton.y >= b[0].y && event.xbutton.y <= b[0].height + b[0].y - 1) {
          XCloseDisplay(dpy); exit(0);
        }
      }
      if (event.xbutton.button == 1 && CreatingZoomWindow) {
        XSetFunction(dpy, gc, GXxor); XSetForeground(dpy, gc, fg^bg);
        ztx = zx1; zty = zy1;
        if (zx2 < ztx) { i = zx2; zx2 = ztx; ztx = i; }
        if (zy2 < zty) { i = zy2; zy2 = zty; zty = i; }
        XDrawRectangle(dpy, win, gc, ztx, zty, zx2-ztx, zy2-zty);
        XSetFunction(dpy, gc, GXcopy); CreatingZoomWindow = False; 
        if (sp < STACKSIZE) {
          stack[sp].row0 = row0; stack[sp].col0 = col0; stack[sp].rdim = rdim;
          stack[sp].cdim = cdim; stack[sp].blocksize=blocksize; stack[sp].xofs = xofs;
          stack[sp].yofs = yofs; sp++;
          zx2 = event.xbutton.x; zy2 = event.xbutton.y;
          if (zx2 < zx1) { i = zx2; zx2 = zx1; zx1 = i; }
          if (zy2 < zy1) { i = zy2; zy2 = zy1; zy1 = i; }
          if (zx2 - zx1 > 4 && zy2 - zy1 > 4) {
            if (blocksize > 0) {
              curc = col0; curr = row0;
              col0 = curc + (zx1 - xofs + blocksize/2)/blocksize;
              row0 = curr + (zy1 - yofs + blocksize/2)/blocksize;
              cdim = curc + (zx2 - xofs + blocksize/2)/blocksize - col0;
              rdim = curr + (zy2 - yofs + blocksize/2)/blocksize - row0;
              if (cdim == 0) { col0 -= 1; cdim += 2; }
              if (rdim == 0) { row0 -= 1; rdim += 2; }
            } else {
              col0 = col0 + (zx1 - xofs)*(-blocksize);
              row0 = row0 + (zy1 - yofs)*(-blocksize);
              cdim = (zx2 - zx1)*(-blocksize);
              rdim = (zy2 - zy1)*(-blocksize);
            }
            if (col0 < 0) { cdim -= -col0; col0 = 0; }
            if (row0 < 0) { rdim -= -row0; row0 = 0; }
            if (cdim < 1) cdim = 1;
            if (rdim < 1) rdim = 1;
            if (col0 >= M->nRows || row0 >= M->nRows) {
              row0 = 0; col0 = 0; rdim = M->nRows; cdim = M->nRows; }
            WriteMatrix(dpy, win, gc, M, file);
            if (stack[sp-1].row0 == row0 && stack[sp-1].col0 == col0 && 
              stack[sp-1].rdim == rdim && stack[sp-1].cdim == cdim && 
              stack[sp-1].blocksize == blocksize && stack[sp-1].xofs == xofs &&
              stack[sp-1].yofs == yofs) sp--;
          }
        }
      }
    } else if (event.type == MotionNotify) {
        strcpy(oldposstring, posstring);
        strcpy(oldvalstring, valstring);
        if (blocksize > 0) {
          curc = col0 + (event.xmotion.x - xofs)/blocksize;
          curr = row0 + (event.xmotion.y - yofs)/blocksize;
        } else {
          curc = col0 + (event.xmotion.x - xofs)*(-blocksize);
          curr = row0 + (event.xmotion.y - yofs)*(-blocksize);
        }
	valstring[0] = 0;
        if (curc < col0 || curc >= col0+cdim || curr < row0 || curr >= row0+rdim)
          strcpy(posstring, "out of range"); else {
	    if ( M->nRows < 100000 )
	      sprintf(posstring, "%5d %5d", curr+1, curc+1);
	    else if ( M->nRows < 1000000 )
	      sprintf(posstring, "%6d %6d", curr+1, curc+1);
	    else if ( M->nRows < 10000000 )
	      sprintf(posstring, "%7d %7d", curr+1, curc+1);
	    else if ( M->nRows < 100000000 )
	      sprintf(posstring, "%8d %8d", curr+1, curc+1);
	    if (blocksize > 0) {
	      for (Col = (M->Row)[curr]; 
		   Col < (M->Row)[curr+1] && (M->Col)[Col] < curc; Col++);
	      if (Col < (M->Row)[curr+1] && (M->Col)[Col] == curc) 
		curv = (M->Val)[Col]; 
	      else curv = 0.0;
	      if BETWEEN( fabs(curv), MINCUT, MAXCUT ) 
		sprintf(valstring, FORMAT3, curv);
	      else sprintf(valstring, FORMAT4, curv);
	    }
	  }
        if (strcmp(oldposstring, posstring) != 0) {
	  int left = xwa.width - 
	    MIN( xwa.width/5, XTextWidth( fontlist[BIGFONT].fontStruct, VERSION, strlen(VERSION) ) ) +
	    leftmargin;
	  XSetFunction( dpy, gc, GXxor ); XSetForeground( dpy, gc, fg^bg );
	  XDrawString( dpy, win, gc, left, xwa.height-3*baselineskip/2, 
		       oldposstring, strlen(oldposstring) );
	  XDrawString( dpy, win, gc, left, xwa.height-3*baselineskip/2, 
		       posstring, strlen(posstring) );
          XSetFunction( dpy, gc, GXcopy );
        }
        if (strcmp(oldvalstring, valstring) != 0) {
	  int left = xwa.width - 
	    MIN( xwa.width/5, XTextWidth( fontlist[BIGFONT].fontStruct, VERSION, strlen(VERSION) ) ) +
	    leftmargin;
          XSetFunction(dpy, gc, GXxor); XSetForeground(dpy, gc, fg^bg);
	  XDrawString( dpy, win, gc, left, xwa.height - baselineskip/2, 
		       oldvalstring, strlen(oldvalstring));
	  XDrawString( dpy, win, gc, left, xwa.height - baselineskip/2, 
		       valstring, strlen(valstring));
          XSetFunction(dpy, gc, GXcopy);
        }
        if (CreatingZoomWindow) {
          XSetFunction(dpy, gc, GXxor); XSetForeground(dpy, gc, fg^bg);
          ztx = zx1; zty = zy1;
          if (zx2 < ztx) { i = zx2; zx2 = ztx; ztx = i; }
          if (zy2 < zty) { i = zy2; zy2 = zty; zty = i; }
          XDrawRectangle(dpy, win, gc, ztx, zty, zx2-ztx, zy2-zty);
          zx2 = event.xmotion.x; zy2 = event.xmotion.y;
          ztx = zx1; zty = zy1;
          if (zx2 < ztx) { i = zx2; zx2 = ztx; ztx = i; }
          if (zy2 < zty) { i = zy2; zy2 = zty; zty = i; }
          XDrawRectangle(dpy, win, gc, ztx, zty, zx2-ztx, zy2-zty);
          XSetFunction(dpy, gc, GXcopy); 
          zx2 = event.xmotion.x; zy2 = event.xmotion.y;
        }
    } else if (event.type == KeyPress) {
      /* printf("%d\n", XLookupKeysym(&event.xkey, event.xkey.state)); */
    }
  }
}
