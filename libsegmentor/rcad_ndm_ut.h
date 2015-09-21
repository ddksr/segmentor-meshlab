/* RCAD_NDM_UT.H */

#ifndef LIBSEGMENTOR_RCAD_NDM_UT_H
#define LIBSEGMENTOR_RCAD_NDM_UT_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifndef RCAD_MMHEAD
#include <malloc.h>
#endif
#include "nrutil.h"

/* 3-MAY-1991 LG */

/* Modifications:
		18-SEP-1991  AL
		10-FEB-1997  LG
                6- MAY-1997  LG
*/


/* From Press & al.: Num. Recipes in C, Cambridge Univ. Press, 1988 */


#undef NULL 
#define NULL 0

/*#ifdef ANSI*/

extern double **dtensor( int, int, int, int, int);
/************************************************/


void free_dtensor( double**, int, int, int, int, int);
/****************************************************/


extern void ptens( double**, int, int, int);
/******************************************/


extern void etens( double**, int);
/********************************/

typedef void (*ndm_errfu)(int);

/*#else

extern double **dtensor();
void free_dtensor();
extern void ptens();
extern void etens();

typedef void (*ndm_errfu)();
#endif */

extern ndm_errfu ndm_exit;

#ifdef RCAD_MMHEAD
extern mmhead *doub_mh;
#else
#define GIVMEM( a, b) malloc(b)
#define RELMEM( a, b) free(b)
#endif

#ifdef	__cplusplus
}
#endif

#endif
