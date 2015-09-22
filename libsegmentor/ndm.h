/* NDM.H */

/* 16-SEP-1991 */

/* Modifications:
		24-SEP-1991  AL
		12-MAR-1992  HT-KGY  -- see MAXDOUBLE
                11-JAN-1994  FARFRBOX,FREEVAL  see relaxtobox in ndm_cm.c
                6- MAY-1997  LG
*/

#ifndef LIBSEGMENTOR_NDM
#define LIBSEGMENTOR_NDM

//#define HUGE 1e50

#ifdef	__cplusplus
extern "C" {
#endif

struct VARDESC
   {double	min,	/* Domain */
   		max,
   		eps;	/* Absolute tolerance of variable */
   };

typedef struct VARDESC vardesc;

#define FUNCONV   0
#define NOFUNCONV 1
#define ERRCONV   0
#define NOERRCONV 2
#define INCONV    0
#define OUTCONV   4
#define TRYOFLOW  8
#define INFCONV   64
#define FLAT      256
#define SUPERFLAT 512
#define FARFRBOX  1024

#define SINGNEWT   -16
#define ININOCONTR -32

#ifdef  __unix
#include    <values.h>
#else
#ifdef  unix
#include    <values.h>
#else
#ifdef  __unix__
#include    <values.h>
#else
#ifndef     MAXDOUBLE
#define     MAXDOUBLE   HUGE
#endif
#endif
#endif
#endif

#define     FREEVAL  MAXDOUBLE

#ifdef  SMALLREAL
#undef  SMALLREAL
#endif
#define SMALLREAL 1.0E-16

#define ON  1
#define OFF 0

#ifdef	__cplusplus
}
#endif

#endif
