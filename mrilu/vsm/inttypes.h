#ifndef INTTYPES_H
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

#include <stdint.h>

typedef  int16_t int16;
typedef  uint16_t uint16;
typedef  int32_t int32;
typedef  uint32_t uint32;
typedef  int64_t  int64;
typedef  uint64_t uint64;

#endif /* INTTYPES_H */
