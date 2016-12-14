#include "globals.h"

Display *currentDisplay;
Drawable currentDrawable;

/******************************************************************************
***** COLOR ADMINISTRATION
******************************************************************************/

/* ColorTab[0] == background
   ColorTab[1..COLORS] the colors */
typedef  uint32  ColorTable[ COLORS + 1 ];

ColorTable colorTable;

#ifdef _PROTOTYPES
uint32 XAllocColorRGB( Display  *display, 
	               Colormap colormap, 
		       unsigned short red, 
		       unsigned short green, 
		       unsigned short blue )
#else
uint32 XAllocColorRGB( display, colormap, red, green, blue )
  Display *display;
  Colormap colormap; 
  unsigned short red; 
  unsigned short green; 
  unsigned short blue;
#endif
{
  XColor screen_in_out;
  
  screen_in_out.red = red;
  screen_in_out.green = green;
  screen_in_out.blue = blue;
  XAllocColor( display, colormap, &screen_in_out );
  return(screen_in_out.pixel);
}

#define STEP(i,ilow,ihigh,clow,chigh) \
  (clow) + ((i) - (ilow))*((chigh) - (clow))/((ihigh) - (ilow))

#ifdef _PROTOTYPES
void initColorTable( int reverse )
#else
void initColorTable( reverse )
     int reverse;
#endif
{
  int i;
  Colormap colormap = DefaultColormap( currentDisplay, DefaultScreen( currentDisplay ) );
  if (reverse) 
    colorTable[0] = XAllocColorRGB( currentDisplay, colormap, 0xffff, 0xffff, 0xffff );
  else colorTable[0] = XAllocColorRGB( currentDisplay, colormap, 0, 0, 0 );
  for ( i = 1; i <= 1 + COLORS/2; i++)
    colorTable[i] = 
      XAllocColorRGB( currentDisplay, colormap, 
		      STEP( i, 1, 1 + COLORS/2, 0, 0xFFFF ),
		      0, 
		      STEP( i, 1, 1 + COLORS/2, 0xFFFF, 0 ) );
  for ( i = 1 + COLORS/2; i <= COLORS; i++ ) 
    colorTable[i] = 
      XAllocColorRGB( currentDisplay, colormap, 
		      0xFFFF, 
		      STEP( i, 1 + COLORS/2, COLORS, 0, 0xFFFF ),
		      0 );
}

/******************************************************************************
***** RECTANGLE ADMINISTRATION
******************************************************************************/

typedef struct {
  XRectangle *rectangles;
  int nRectangles;
  int maxRectangles;
} RectangleList;

typedef RectangleList RectangleListForColor[ COLORS + 1 ];

RectangleListForColor rectangleListForColor;

#ifdef _PROTOTYPES
void initRectangleList( RectangleList *pRectangleList )
#else
void initRectangleList( pRectangleList )
     RectangleList *pRectangleList;
#endif
{
  int32  maxRectangles = 4*(XMaxRequestSize( currentDisplay ) - 3)
                         / sizeof( XRectangle );
  pRectangleList->rectangles = NEW( maxRectangles, XRectangle );
  pRectangleList->nRectangles = 0;
  pRectangleList->maxRectangles = maxRectangles;
}

#ifdef _PROTOTYPES
void initRectangleListForColor( RectangleListForColor rectangleForColorList ) 
#else
void initRectangleListForColor( rectangleForColorList )
     RectangleListForColor rectangleForColorList;
#endif
{
  Color color;
  for ( color = 0; color <= COLORS; color++ ) {
    initRectangleList( &rectangleListForColor[ color ] );
  }
}

#ifdef _PROTOTYPES
void flushRectangleList( Color color )
#else
void flushRectangleList( color )
     Color color;
#endif
{
  GC gc = DefaultGC( currentDisplay, DefaultScreen( currentDisplay ) );
  RectangleList *pRectangleList = &rectangleListForColor[ color ];
  if ( pRectangleList->nRectangles > 0 ) {
    XSetForeground( currentDisplay, gc, colorTable[ color ] );
    XFillRectangles( currentDisplay, currentDrawable, gc, 
		     pRectangleList->rectangles, 
		     pRectangleList->nRectangles );
    pRectangleList->nRectangles = 0;
  }
}

#ifdef _PROTOTYPES
void LXfillRectangle( short x, short y, 
		      unsigned short width, unsigned short height, 
		      Color color )
#else
void LXfillRectangle( x, y, width, height, color )
  short x, y; 
  unsigned short width, height;
  Color color;
#endif
{
  RectangleList *pRectangleList = &rectangleListForColor[ color ];
  if ( pRectangleList->nRectangles == pRectangleList->maxRectangles ) 
    flushRectangleList( color );
  {
    XRectangle *pXRectangle = &pRectangleList->rectangles[ pRectangleList->nRectangles ];
    pXRectangle->x = x;
    pXRectangle->y = y;
    pXRectangle->width = width;
    pXRectangle->height = height;
  }
  pRectangleList->nRectangles++;
}

/******************************************************************************
***** POINT ADMINISTRATION
******************************************************************************/

typedef struct {
  XPoint *points;
  int nPoints;
  int maxPoints;
} PointList;

typedef PointList PointListForColor[ COLORS + 1 ];

PointListForColor pointListForColor;

#ifdef _PROTOTYPES
void initPointList( PointList *pPointList )
#else
void initPointList( pPointList )
     PointList *pPointList;
#endif
{
  int32  maxPoints = 4*(XMaxRequestSize( currentDisplay ) - 3)
                     / sizeof( XPoint );
  pPointList->points = NEW( maxPoints, XPoint );
  pPointList->nPoints = 0;
  pPointList->maxPoints = maxPoints;
}

#ifdef _PROTOTYPES
void initPointListForColor( PointListForColor pointForColorList ) 
#else
void initPointListForColor( pointForColorList ) 
     PointListForColor pointForColorList;
#endif
{
  Color color;
  for ( color = 0; color <= COLORS; color++ ) {
    initPointList( &pointListForColor[ color ] );
  }
}

#ifdef _PROTOTYPES
void flushPointList( Color color )
#else
void flushPointList( color )
     Color color;
#endif
{
  GC gc = DefaultGC( currentDisplay, DefaultScreen( currentDisplay ) );
  PointList *pPointList = &pointListForColor[ color ];
  if ( pPointList->nPoints > 0 ) {
    XSetForeground( currentDisplay, gc, colorTable[ color ] );
    XDrawPoints( currentDisplay, currentDrawable, gc, 
		 pPointList->points, pPointList->nPoints,
		 CoordModeOrigin );
    pPointList->nPoints = 0;
  }
}

#ifdef _PROTOTYPES
void LXdrawPoint( short x, short y, Color color )
#else
void LXdrawPoint( x, y, color )
     short x, y;
     Color color;
#endif
{
  PointList *pPointList = &pointListForColor[ color ];
  if ( pPointList->nPoints == pPointList->maxPoints ) 
    flushPointList( color );
  {
    XPoint *pXPoint = &pPointList->points[ pPointList->nPoints ];
    pXPoint->x = x;
    pXPoint->y = y;
  }
  pPointList->nPoints++;
}

#ifdef _PROTOTYPES
void LXflushAll()
#else
void LXflushAll()
#endif
{
  Color color;
  for ( color = 0; color <= COLORS; color++ ) { 
    flushPointList( color );
    flushRectangleList( color );
  }
}

#ifdef _PROTOTYPES
void LXinit( Display *display, Drawable drawable, int reverse )
#else
void LXinit( display, drawable, reverse )
     Display *display;
     Drawable drawable;
     int reverse;
#endif
{
  currentDisplay = display; currentDrawable = drawable;
  initRectangleListForColor( rectangleListForColor );
  initPointListForColor( pointListForColor );
  initColorTable( reverse );
}


