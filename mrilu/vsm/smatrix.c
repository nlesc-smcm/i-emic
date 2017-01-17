#include "globals.h"

/******************************************************************************
***** CURRENT FILE ADMINISTRATION
******************************************************************************/

struct {
  FILE *file;
  char name[1024];
} currentFile;

#ifdef _PROTOTYPES
void fatal( int error )
#else
void fatal( error )
     int error;
#endif
{
  char buffer[1024];
  switch ( error ) {
  case 0: sprintf( buffer, "Error opening file \"%s\"", currentFile.name ); break;
  case 1: sprintf( buffer, "Error reading file \"%s\"", currentFile.name ); break;
  case 2: sprintf( buffer, "Error closeing file \"%s\"", currentFile.name ); break;
  case 3: sprintf( buffer, "Error checking file \"%s\"", currentFile.name ); break;
  case 4: sprintf( buffer, "Error removing file \"%s\"", currentFile.name ); break;
  case 5: sprintf( buffer, "Error getting size of file \"%s\"", currentFile.name ); break;
  default: sprintf( buffer, "Unknown error for file \"%s\"", currentFile.name ); break;
  }
  perror( buffer );
  exit( 1 );
}

#ifdef _PROTOTYPES
void CFopen( char *pattern, char *arg, char *type )
#else
void CFopen( pattern, arg, type )
     char *pattern;
     char *arg;
     char *type;
#endif
{
  sprintf( currentFile.name, pattern, arg );
  if ( ( currentFile.file = fopen( currentFile.name, type ) ) == NULL ) fatal(0);
}

#ifdef _PROTOTYPES
void CFremove( char *pattern, char *arg )
#else
void CFremove( pattern, arg )
     char *pattern;
     char *arg;
#endif
{
  sprintf( currentFile.name, pattern, arg );
  if ( unlink( currentFile.name ) < 0 ) fatal(4);
}

#ifdef _PROTOTYPES
void CFreadBuffer( void *ptr, size_t size, size_t nitems )
#else
void CFreadBuffer( ptr, size, nitems )
     void *ptr;
     size_t size;
     size_t nitems;
#endif
{
  if ( fread( ptr, size, nitems, currentFile.file ) == 0 ) fatal(1);
}

#ifdef _PROTOTYPES
double CFreadDouble()
#else
double CFreadDouble()
#endif
{
  double buffer;
  if ( fscanf( currentFile.file, "%lf", &buffer ) == EOF ) fatal(1);
  return buffer;
}

#ifdef _PROTOTYPES
double CFreadInt32(void)
#else
double CFreadInt32()
#endif
{
  return (int32) floor( CFreadDouble() + 0.5 );
}

#ifdef _PROTOTYPES
void CFclose()
#else
void CFclose()
#endif
{
  if ( fclose( currentFile.file ) == EOF ) fatal(2);
}

/******************************************************************************
***** READING AND REMOVING A SPARSE MATRIX
******************************************************************************/

/* check and correct if the column indixes are not ascending */
#ifdef _PROTOTYPES
void sortOneRow( SparseMatrix *M, int left, int right1 )
#else
void sortOneRow( M, left, right1 )
     SparseMatrix *M;
     int left;
     int right1;
#endif
{
  int i, j, ti, smallest;
  double td;
  for ( i = left; i < right1 - 1; i++ ) {
    smallest = i;
    for ( j = i + 1; j < right1; j++ )
      if ( M->Col[smallest] > M->Col[j] )
	smallest = j;
    if ( smallest != i ) { /* swap i, smallest */
      ti = M->Col[i]; M->Col[i] = M->Col[smallest]; M->Col[smallest] = ti;
      td = M->Val[i]; M->Val[i] = M->Val[smallest]; M->Val[smallest] = td;
    }
  }
}

/* remove (double) 0.0 from the sparse matrix, before AbsMax and AbsMin */
#ifdef _PROTOTYPES
void removeZeros( SparseMatrix *M )
#else
void removeZeros( M )
     SparseMatrix *M;
#endif
{
  int32  row, col;
  int32  colwrite = 0;
  int32  start = M->Row[0];
  for ( row = 0; row < M->nRows; row++ ) {
    for ( col = start; col < M->Row[row+1]; col++ ) {
      if ( M->Val[col] != (double) 0.0 ) {
	M->Col[colwrite] = M->Col[col];
	M->Val[colwrite] = M->Val[col];
	colwrite++;
      }
    }
    start = M->Row[row+1];
    M->Row[row+1] = colwrite;
  }
  M->nVals = M->Row[M->nRows];
}

/* make the matrix a square */
#ifdef _PROTOTYPES
void makeSquare( SparseMatrix *M, int32 nCols )
#else
void makeSquare( M, nCols )
     SparseMatrix *M;
     int32 nCols;
#endif
{
  int32  *newRow = NEW( nCols + 1, int32 );
  int i;
  for ( i = 0; i <= M->nRows; i++ )
    newRow[i] = M->Row[i];
  for ( i = M->nRows + 1; i <= nCols; i++ )
    newRow[i] = M->Row[ M->nRows ];
  free( M->Row );
  M->Row = newRow;
  M->nRows = nCols;
}

#ifdef _PROTOTYPES
void CheckSparseMatrix( SparseMatrix *M )
#else
void CheckSparseMatrix( M )
     SparseMatrix *M;
#endif
{
  int row, i, isAscending;
  /* check if all rows or in ascending column order */
  for ( row = 0; row < M->nRows; row++ ) {
    /* I expect the row is already sorted, so first check for this */
    for ( i = M->Row[row], isAscending = True; i < M->Row[row + 1] - 1; i++ )
      if ( M->Col[i] > M->Col[i+1] ) { isAscending = False; break; }
    if ( !isAscending )
      sortOneRow( M, M->Row[row], M->Row[row+1] );
  }
  /* check if there are no 0's in the matrix */
  {
    int hasZeros = False;
    for ( i = 0; i < M->nVals; i++ )
      if ( M->Val[i] == (double) 0.0 )
	hasZeros = True;
    if ( hasZeros )
      removeZeros( M );
  }
  /* check if the matrix is a square, and if not, make it a square */
  {
    int32  maxcol = 0;
    for ( i = 0; i < M->nVals; i++ )
      if ( M->Col[i] > maxcol )
	maxcol = M->Col[i];
    if ( maxcol + 1 > M->nRows ) {
      makeSquare( M, maxcol + 1 );
    }
  }
}

/* Read the SparseMatrix M from the files `C' */
#ifdef _PROTOTYPES
void ReadSparseMatrix( SparseMatrix *M, char *C, int method )
#else
void ReadSparseMatrix( M, C, method )
     SparseMatrix *M;
     char *C;
     int method;
#endif
{
  int32  i;
  double dtmp;
  int32  Offset;
  struct stat buf;
  int32  size1, size2;

  M->file.name = NEW( strlen(C) + 1, char );
  strcpy( M->file.name, C );
  M->file.kind = method;

  /* read M->nRows */
  switch( method ) {
  case METHOD_ASCII:
    CFopen( "%s.beg", C, "r" );
    i = 0; while ( fscanf( currentFile.file, "%lf", &dtmp ) != EOF ) i++;
    M->nRows = i - 1;
    CFclose();
    break;
  case METHOD_ASCII_SINGLE:
    CFopen( "%s", C, "r" );
    M->nRows = CFreadInt32();
    /*CFclose(); */
    break;
  case METHOD_BINARY:
    CFopen( "%s.beg", C, "rb" );
    if ( fstat( fileno( currentFile.file ), &buf ) < 0 ) fatal(5);
    M->nRows = buf.st_size/sizeof(int32) - 1;
    break;
  case METHOD_FORTRAN:
    CFopen( "%s", C, "rb" );
    CFreadBuffer( &size1, sizeof(size1), 1 );
    CFreadBuffer( &M->nRows, sizeof(M->nRows), 1 );
    CFreadBuffer( &size2, sizeof(size2), 1 );
    if ( ( size1 != size2 ) || ( size1 != sizeof(M->nRows) ) ) fatal(3);
    break;
  }

  M->nRowsOrig = M->nRows;

  /* read M->Row */
  M->Row = NEW( M->nRows + 1, int32 );
  switch( method ) {
  case METHOD_ASCII:
    CFopen( "%s.beg", C, "r" );
    for ( i = 0; i < M->nRows + 1; i++ ) M->Row[i] = CFreadInt32();
    CFclose();
    break;
  case METHOD_ASCII_SINGLE:
    for ( i = 0; i < M->nRows + 1; i++ ) M->Row[i] = CFreadInt32();
    break;
  case METHOD_BINARY:
    CFreadBuffer( M->Row, sizeof(int32), M->nRows + 1 );
    CFclose();
    break;
  case METHOD_FORTRAN:
    CFreadBuffer( &size1, sizeof(size1), 1 );
    CFreadBuffer( M->Row, sizeof(int32), M->nRows + 1 );
    CFreadBuffer( &size2, sizeof(size2), 1 );
    if ( ( size1 != size2 ) || ( size1 != ( M->nRows + 1 )*sizeof(int32) ) )
      fatal(3);
    break;
  }

  Offset = M->Row[0];
  for( i = 0; i <= M->nRows; i++ ) M->Row[i] -= Offset;

  M->nVals = M->Row[M->nRows];

  /* read M->Col */
  M->Col = NEW( M->nVals, int32 );
  switch( method ) {
  case METHOD_ASCII:
    CFopen( "%s.jco", C, "r" );
    for ( i = 0; i < M->nVals; i++ ) M->Col[i] = CFreadInt32();
    CFclose();
    break;
  case METHOD_ASCII_SINGLE:
    for ( i = 0; i < M->nVals; i++ ) M->Col[i] = CFreadInt32();
    break;
  case METHOD_BINARY:
    CFopen( "%s.jco", C, "rb" );
    CFreadBuffer( M->Col, sizeof(int32), M->nVals );
    CFclose();
    break;
  case METHOD_FORTRAN:
    CFreadBuffer( &size1, sizeof(size1), 1 );
    CFreadBuffer( M->Col, sizeof(int32), M->nVals );
    CFreadBuffer( &size2, sizeof(size2), 1 );
    if ( ( size1 != size2 ) || ( size1 != ( M->nVals )*sizeof(int32) ) )
      fatal(3);
  }

  for ( i = 0; i < M->nVals; i++ ) M->Col[i] -= Offset;

  /* read M->Val */
  M->Val = NEW( M->nVals, double );
  switch( method ) {
  case METHOD_ASCII:
    CFopen( "%s.co", C, "r" );
    for ( i = 0; i < M->nVals; i++ ) M->Val[i] = CFreadDouble();
    CFclose();
    break;
  case METHOD_ASCII_SINGLE:
    for ( i = 0; i < M->nVals; i++ ) M->Val[i] = CFreadDouble();
    CFclose();
    break;
  case METHOD_BINARY:
    CFopen( "%s.co", C, "rb" );
    CFreadBuffer( M->Val, sizeof(double), M->nVals );
    CFclose();
    break;
  case METHOD_FORTRAN:
    CFreadBuffer( &size1, sizeof(size1), 1 );
    CFreadBuffer( M->Val, sizeof(double), M->nVals );
    CFreadBuffer( &size2, sizeof(size2), 1 );
    if ( ( size1 != size2 ) || ( size1 != ( M->nVals )*sizeof(double) ) ) fatal(3);
    break;
  }

  CheckSparseMatrix( M );

  M->AbsMax = M->AbsMin = fabs( M->Val[0] );
  for( i = 1; i < M->nVals; i++ )  {
    dtmp = fabs( M->Val[i] );
    if ( dtmp > M->AbsMax ) M->AbsMax = dtmp;
    else if ( dtmp < M->AbsMin ) M->AbsMin = dtmp;
  }

  M->Color = NEW( M->nVals, int32 );
} /* ReadSparseMatrix */

/* Remove the file `C' */
#ifdef _PROTOTYPES
void RemoveSparseMatrix ( SparseMatrix *M )
#else
void RemoveSparseMatrix ( M )
     SparseMatrix *M;
#endif
{
  switch( M->file.kind ) {
  case METHOD_ASCII:
  case METHOD_BINARY:
    CFremove( "%s.beg", M->file.name );
    CFremove( "%s.jco", M->file.name );
    CFremove( "%s.co", M->file.name );
    break;
  case METHOD_FORTRAN:
  case METHOD_ASCII_SINGLE:
    CFremove( "%s", M->file.name );
    break;
  }
}


