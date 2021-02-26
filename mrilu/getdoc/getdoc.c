/*                  getdoc.c                                            */
/*                  ============                                            */
/*                                                                          */
/*  Extract Document from Fortran files                                     */
/*        $Id: getdoc.c,v 1.1.1.1 2006/10/02 13:33:06 wubs Exp $    */
/*                                                                          */
/* This program extracts the documentation lines from Fortran programs.     */
/* These documentation lines can consist of comment lines and statement     */
/* lines.  The lines to be extracted start after the comment line:          */
/*   C#begindoc                                                             */
/* and end just before the comment line:                                    */
/*   C#enddoc                                                               */
/* Only the extracted lines are copied, verbatim, to the output file.       */
/*                                                                          */
/* The program can be called in the following tree ways:                    */
/*   getdoc                                                                 */
/*      The lines are read from the standard input file and the extracted   */
/*      lines are copied to the standard output file.                       */
/*   getdoc <inpfilenm>                                                     */
/*      The filename <inpfilenm> should have the suffix ".f" or ".F".       */
/*      The output filename is <inpfilenm> in which the suffix has changed  */
/*      to  ".txt".                                                         */
/*      The lines are read from the input file <inpfilenm> and the extracted*/
/*      lines are copied to the output file.                                */
/*   getdoc <inpfilenm> <outpfilenm>                                        */
/*      The filename <inpfilenm> should have the suffix ".f" or ".F".       */
/*      The lines are read from the input file <inpfilenm> and the extracted*/
/*      lines are copied to the output file <outpfilenm>.                   */
/*                                                                          */
/*        1997-01-10  ,  Doeke de Vries                                     */
/*        1999-02-19     Correction to handle unsigned characters.          */


#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gendefs.h"


#define       BEGINDOC   "c#begindoc"
#define       ENDDOC     "c#enddoc"
#define       NEWBEGINDOC   "!#begindoc"
#define       NEWENDDOC     "!#enddoc"


static FILE*  inpf;           /* Stream input file */
static FILE*  outpf;          /* Stream output file */

static int    lach;           /* Look ahead character from input file.   */
                              /* The output value of function  fgetc  is */
                              /* of type  int.                           */

static void printmsg (
  char* msg
)
/* Write  <errprefix><msg>!  to the file(s) */
{
  fputc ('\n', stderr);
  fputs (msg , stderr);
  fputc ('!' , stderr);
  fputc ('\n', stderr);
}


void fatalerror (
  char* fmt, ...
)
/* Write  "<fmt, ...>!"  to stderr.   Terminate the program. */
{
  va_list args;

  fprintf(stderr, "\n");
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fprintf(stderr, "!\n");

  abort();  /* Can be caught by a handler in order to cleanup */
} 


void internalerror (
  char* mnm,
  char* fmt, ...
)
/* Write  "Internal error in <mnm>: <fmt, ...>!"  to the file(s). */
/* Terminate the program.                                         */
{
  char    errstr[MAX_STRLEN];
  char    str[MAX_STRLEN];
  va_list args;

  va_start(args, fmt);
  vsprintf(str, fmt, args);
  va_end(args);

  sprintf(errstr, "Interne fout in %s: %s", mnm, str);
  printmsg(errstr);
  abort();  /* Can be caught by a handler in order to cleanup */
}


void toomany (
  char* mnm,
  char* str,
  int   nr
)
/* Write  "In <mnm>: too many <str>(<nr>)!"  to the file(s). */
/* Terminate the program.                                    */
{
  char errstr[MAX_STRLEN];

  sprintf(errstr, "In %s: too many %s(%1i)", mnm, str, nr);
  printmsg(errstr);
  abort();  /* Can be caught by a handler in order to cleanup */
}








void init (
  int argc,
  char* argv[]
)
/* Read and check the command line parameters.  Initialise the input and    */
/* output files (define the global variables  inpf  and  outpf).            */
/* Define the lookahead character in  lach  as the first character of the   */
/* file  inpf.                                                              */
{ char *fnm;
  int  len;

  if ( argc > 3 )
    fatalerror ("Usage:  getdoc [ inpfnm [ outpfnm ] ]");

  if ( argc == 1 ) {
    inpf  = stdin;
    outpf = stdout;
  } else {
    /*  argc == 2  or  argc == 3  */
    len = strlen (argv[1]);
    if ( len + (3 - 1) > MAX_STRLEN )
      fatalerror ("Too many characters in file name %s", argv[1]);
    fnm = argv[1];
    if (   ( len < 3  ||  fnm[len-2] != '.'  ||  tolower(fnm[len-1]) != 'f' )
	   && ( len < 5  ||  fnm[len-4] != '.'  ||  tolower(fnm[len-3]) != 'f' ||  fnm[len-2] != '9'||  fnm[len-1] != '0'  ) )
      fatalerror ("No extension .f, .F, .f90 or .F90 in file name %s", fnm);
    inpf = fopen (fnm, "r");
    if ( inpf == NULL )
      fatalerror ("Input file %s: %s", fnm, strerror (errno));

    if ( argc == 2 ) {
      int maxLen = len + 4;
      char *tmpName = NULL;

      tmpName = malloc(maxLen);
      if (tmpName == NULL) {
        fatalerror ("Failed to allocate memory in init");
      }

      strncpy(tmpName, fnm, len);
      strncat(tmpName, "txt", maxLen - 1);

      outpf = fopen (tmpName, "w");
      free(tmpName);
    } else {
      outpf = fopen (argv[2], "w");
    }
    if ( outpf == NULL )
      fatalerror ("Output file %s: %s", fnm, strerror (errno));
  }

  lach = fgetc (inpf);
}


int read_1st_word (
  char firstword[]
)
/* Read the first word of the actual line from the input file.  Store the   */
/* copy of the first word in parameter  firstword  and return the number of */
/* characters stored.                                                       */
/* lach = EOF  or  "lach  is 1st char. on actual line after the 1st word"   */
{ int i = 0;

  while ( ! isspace(lach)  &&  lach != EOF  &&  i < MAX_STRLEN ) {
    firstword[i] = lach;
    lach = fgetc (inpf);
    i += 1;
  }
  firstword[i] = '\0';
  return  i;
}


void skip_rest_line (void)
/* Skip the rest of the actual line from the file  inpf, i.e. read past the */
/* EndOfLine on the actual line or until EndOfFile.                         */
/* lach == EOF  or  "lach is first character of next line"                  */
{
  while ( lach != '\n'  &&  lach != EOF ) {
    lach = fgetc (inpf);
  }
  if ( lach == '\n' ) {
    lach = fgetc (inpf);
  }
}


void copy_word_rest_line (
  char firstword[]
)
/* Copy the word  firstword  and the rest of the actual line from the file  */
/* inpf  to the file  outpf, i.e. copy past the EndOfLine on the actual line*/
/* or until EndOfFile.                                                      */
/* lach == EOF  or  "lach is first character of next line"                  */
{
  fputs (firstword, outpf);
  while ( lach != '\n'  &&  lach != EOF ) {
    fputc (lach, outpf);
    lach = fgetc (inpf);
  }
  fputc ('\n', outpf);
  if ( lach == '\n' ) {
    lach = fgetc (inpf);
  }
}


int main (
  int argc,
  char* argv[]
)
{
  int  len;                   /* Length of a word <= MAX_STRLEN          */
  int  state = 0;             /* state = 0: skip line from input file    */
                              /* state > 0: copy line from input file to */
                              /*            output file                  */
  int  i;    

  char firstword[MAX_STRLEN+1];
  char lowerword[MAX_STRLEN+1];

  Assert (     MAX_STRLEN >= strlen (BEGINDOC)
           &&  MAX_STRLEN >= strlen (ENDDOC)   );

  init (argc, argv);

  while ( lach != EOF ) {
    len = read_1st_word (firstword);
    for ( i = 0;  i < len;  i += 1 )
      lowerword[i] = tolower(firstword[i]);
    lowerword[len] = '\0';
    if ( strcmp(lowerword, BEGINDOC) == 0 || strcmp(lowerword, NEWBEGINDOC) == 0) {
      state += 1;  skip_rest_line ();
    } else if ( strcmp(lowerword, ENDDOC) == 0 || strcmp(lowerword, NEWENDDOC) == 0) {
      if ( state == 0 )
        fatalerror ("No matching  %s  for  %s", BEGINDOC, firstword);
      state -= 1;  skip_rest_line ();
    } else if ( state <= 0 ) {
      skip_rest_line ();
    } else {
      copy_word_rest_line (firstword);
    }
  }  
  if ( state != 0 )
    fatalerror ("No matching  %s  at end of file", ENDDOC);
  fflush (outpf);

  return  EXIT_SUCCESS;
}
