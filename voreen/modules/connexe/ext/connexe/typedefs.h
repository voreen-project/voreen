/*************************************************************************
 * typedefs.h -
 *
 * $Id: typedefs.h,v 1.2 1999/12/08 16:35:14 greg Exp $
 *
 * Copyright©INRIA 1998
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 *
 * CREATION DATE:
 * June, 9 1998
 *
 *
 *
 */

#ifndef _typedefs_h_
#define _typedefs_h_

#ifdef __cplusplus
extern "C" {
#endif





/* Differents type coding for images and buffers.
 */
typedef enum {
  TYPE_UNKNOWN /* unknown type */,
  UCHAR_TYPE  /* unsigned char */,
  SCHAR_TYPE  /* signed char */,
  USHORT_TYPE /* unsigned short int */,
  SSHORT /* signed short int */,
  INT_TYPE    /* signed int */,
  ULINT_TYPE  /* unsigned long int */,
  FLOAT_TYPE  /* float */,
  DOUBLE_TYPE  /* double */
} ImageType, bufferType;

typedef char               s8;
typedef unsigned char      u8;
typedef short int          s16;
typedef unsigned short int u16;
typedef int                i32;
typedef int                s32;
typedef unsigned long int  u64;
typedef float              r32;
typedef double             r64;





/* Typedef Booleen
 */
typedef enum {
  False = 0,
  True = 1
} typeBoolean;


#ifdef __cplusplus
}
#endif

#endif
