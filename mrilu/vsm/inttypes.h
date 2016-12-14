#if !defined(INTTYPES_H)
#define INTTYPES_H

/*                              inttypes.h
 *                              ==========
 *
 * Integer type definitions for C-programs.
 *
 * The following integer types with their maximum value are defined:
 *    int16   integer          with a 16 bit representation  INT16_MAX
 *    uint16  unsigned integer with a 16 bit representation  UINT16_MAX
 *    int32   integer          with a 32 bit representation  INT32_MAX
 *    uint32  unsigned integer with a 32 bit representation  UINT32_MAX
 *
 * The following integer types are defined, if possible:
 *    int64   integer          with a 64 bit representation  INT64_MAX
 *    uint64  unsigned integer with a 64 bit representation  UINT64_MAX
 *
 * 2000-07-17  ,  Doeke de Vries
 */

#include <limits.h>

#if SHRT_MAX == 32767
   typedef  short           int16;
   typedef  unsigned short  uint16;
   #define  INT16_MAX       SHRT_MAX     /* 32767  */
   #define  UINT16_MAX      USHRT_MAX    /* 65535U */
#elif INT_MAX == 32767
   typedef  int             int16;
   typedef  unsigned int    uint16;
   #define  INT16_MAX       INT_MAX      /*  32767  */
   #define  UINT16_MAX      UINT_MAX     /*  65535U */
#endif

#if INT_MAX == 2147483647
   typedef  int             int32;
   typedef  unsigned int    uint32;
   #define  INT32_MAX       INT_MAX      /* 2147483647  */
   #define  UINT32_MAX      UINT_MAX     /* 4294967295U */
#elif LONG_MAX == 2147483647L
   typedef  long            int32;
   typedef  unsigned long   uint32;
   #define  INT32_MAX       LONG_MAX     /* 2147483647L  */
   #define  UINT32_MAX      ULONG_MAX    /* 4294967295UL */
#endif

#if LONG_MAX > 2147483647L
   #if INT_MAX == 9223372036854775807L
      typedef  int            int64;
      typedef  unsigned int   uint64;
      #define  INT64_MAX      INT_MAX     /*  9223372036854775807  */
      #define  UINT64_MAX     UINT_MAX    /* 18446744073709551615U */
   #elif LONG_MAX == 9223372036854775807L
      typedef  long           int64;
      typedef  unsigned long  uint64;
      #define  INT64_MAX      LONG_MAX    /*  9223372036854775807L  */
      #define  UINT64_MAX     ULONG_MAX   /* 18446744073709551615UL */
   #endif
#endif


#endif /* INTTYPES_H */
