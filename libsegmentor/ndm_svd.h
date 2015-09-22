/* NDM_SVD.H */

/* 27-Aug-1991 LG */

/* Modifications:
                6- MAY-1997  LG
		18-SEP-1991  AL
*/

#ifndef LIBSEGMENTOR_NDM_SVD_H


#define LIBSEGMENTOR_NDM_SVD_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef ANSI

extern void init_svdcmp( int, int, char**);
/*****************************************/


extern void svbksb( double**, double[], double**, int, int,
                    double[], double[], char**);
/*********************************************************/


extern int svdcmp( double**, int, int, double*, double**, char**);
/*****************************************************************/

#else

extern void init_svdcmp();
extern void svbksb();
extern int svdcmp();

#endif

#ifdef	__cplusplus
}
#endif

#endif







