#if !defined(GENDEFS_H)
#define GENDEFS_H

/*                              gendefs.h
 *                              =========
 *
 *  General definitions for C-programs.
 *        $Id: gendefs.h,v 1.1.1.1 2006/10/02 13:33:06 wubs Exp $
 *
 *        14-Dec-1993 ,  Doeke de Vries
 *        08-Aug-1996 ,  addition of  __Assert
 */

#define  EOL           '\n'      /* End Of Line character      */


#define  FALSE         0         /* Constanten van het type    */
#define  TRUE          1         /* Boolean                    */

typedef  int           Boolean;


#define  MAX_STRLEN    128       /* Max. lengte van de strings */


#  undef Assert

#  ifdef OPTIM
#    define Assert(_EX)  ((void)0)
#  else /* not OPTIM */

     extern void __Assert(char *, int);

#    define Assert(_EX) \
          ((_EX) ? (void)0 : __Assert( __FILE__, __LINE__))

#  endif /* not OPTIM */


#endif /* GENDEFS_H */
