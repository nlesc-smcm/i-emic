/*                              gendefs.c
 *                              =========
 *
 *  General definitions for C-programs.
 *        $Id: gendefs.c,v 1.1.1.1 2006/10/02 13:33:06 wubs Exp $
 *
 *        14-Dec-1993 ,  Doeke de Vries
 *        08-Aug-1996 ,  addition of  __Assert
 *        24-Dec-1996 ,  write error message to the file stderr.
 *        18-Feb-1999 ,  include  stdlib.h  with prototype of  abort.
 */


#include <stdio.h>
#include <stdlib.h>

#include "gendefs.h"


void __Assert (
  char* filenm,
  int   linenr
)
/* Service routine for  Assert macro */
{
  fprintf (stderr, "%s:%u: failed assertion\n", filenm, linenr);
  abort ();
}
