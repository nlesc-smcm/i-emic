#include "globals.h"

#ifdef _PROTOTYPES
int main ( int argc, char *argv[] )
#else
void main (argc, argv)
     int argc;
     char *argv[];
#endif
{
  SparseMatrix M;
  int errflg = 0;
  int aflg = 0, Aflg = 0, bflg = 0, dflg = 0, rflg = 0;
  extern char *optarg;
  extern int optind;
  int c;
  int method;

  while ( ( c = getopt(argc, argv, "aAbdr") ) != -1 ) {
    switch (c) {
    case 'a': if (bflg||Aflg) errflg++; else aflg++;
      break;
    case 'A': if (bflg||aflg) errflg++; else Aflg++;
      break;
    case 'b': if (aflg||Aflg)	errflg++; else bflg++;
      break;
    case 'd': dflg++;
      break;
    case 'r': rflg++;
      break;
    case '?': errflg++;
      break;
    }
  }
  if (errflg || optind != argc - 1) {
    fprintf(stderr, "%s (%s) -- Visualize Sparse Matrix\n", VERSION, DATE);
    fprintf(stderr, "(University of Groningen, Department of Mathematics)\n\n");
    fprintf(stderr, "USAGE: %s [-a|-A|-b][-d][-r] <matrixid>\n\n", argv[0]);

    fprintf(stderr, "OPTION: -a: input in ascii format\n");
    fprintf(stderr, " FILES: `matrixid'.co `matrixid'.jco `matrixid'.beg\n");
    fprintf(stderr, "OPTION: -A: input in ascii single file format\n");
    fprintf(stderr, " FILES: only `matrixid'\n");
    fprintf(stderr, "OPTION: -b: input in binary format\n");
    fprintf(stderr, " FILES: `matrixid'.co `matrixid'.jco `matrixid'.beg\n");
    fprintf(stderr, "NOT -a, -A OR -b: input in fortran unformatted format\n");
    fprintf(stderr, " FILES: only `matrixid'\n\n");
    fprintf(stderr, "OPTION: -d: remove file(s) afterwards\n\n");
    fprintf(stderr, "OPTION: -r: reverse background color (white)\n\n");
    fprintf(stderr, "MOUSE: Left Button - create zoom window\n");
    fprintf(stderr, "       Right Button - unzoom\n");
    exit(1);
  }

  if (aflg) method = METHOD_ASCII;
  else if (bflg) method = METHOD_BINARY;
  else if (Aflg) method = METHOD_ASCII_SINGLE;
  else method = METHOD_FORTRAN;
  ReadSparseMatrix( &M, argv[optind], method );
  VisualizeSparseMatrix( &M, rflg > 0 );
  if (dflg) RemoveSparseMatrix( &M );
  exit(0);
}

