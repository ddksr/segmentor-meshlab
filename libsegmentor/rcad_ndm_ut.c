/* RCAD_NDM_UT.C */

/* 8-Aug-1991  LG */

/* Modifications:
		18-SEP-1991  AL
		07-NOV-1991  LG
		10-FEB-1997  LG
		21-APR-1997  LG ANSI/ not ANSI
*/


/* From Press & al.: Num. Recipes in C, Cambridge Univ. Press, 1988 */
/* No mmhead just malloc-free. */

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

#include "rcad_ndm_ut.h"


#ifdef RCAD_MMHEAD
static mmhead default_mmh = {0};

mmhead *doub_mh = &default_mmh;
#endif

ndm_errfu ndm_exit = exit;


#ifndef ANSI
double **dtensor( r, ml, mh, nl, nh)
/**********************************/
int r,ml,mh,nl,nh;
#else
double **dtensor(int r,int ml,int mh,int nl,int nh)
/*************************************************/
#endif
/* Any address is of equal length!!!!! */
/* Allocates a r-index tensor of size (mh-ml)*(nh-nl)^(r-1). */
/* Range [ml...mh] [nl...nh] ... [nl...nh] */

{
    int i;
   double *rd;
   double **m;

   if (r <= 1) {
      rd = (double *) GIVMEM( doub_mh,
			        (unsigned) (mh-ml+1)*sizeof(double));
      if (rd == 0) nrerror("allocation failure 1 in dtensor()");
      rd -= ml;
      return((double **) rd);
   } else {
      /* Allocate pointers to rows */
      m = (double **) GIVMEM( doub_mh,
			     (unsigned) (mh-ml+1)*sizeof(double *));
      if (m == 0) nrerror("allocation failure 2 in dtensor()");
      m -= ml;
      /* Allocate rows and set pointers to them */
      for (i = ml; i <= mh; i++)
	 m[i] = (double *) dtensor ( r-1, nl, nh, nl, nh);	/* !!! */
      /* Return pointer to array of pointers of rows */
      return( m);
   }
}

#ifndef ANSI
void free_dtensor( t, r, ml, mh, nl, nh)
/**************************************/

double **t;
int    r, ml, mh, nl, nh;
#else
void free_dtensor(double **t,int r,int ml,int mh,int nl,int nh)
/**************************************/
#endif
/* Any address is of equal length!!!!! */
/* Allocates a r-index tensor of size (mh-ml)*(nh-nl)^(r-1). */
/* Range [ml...mh] [nl...nh] ... [nl...nh] */

{
    int i;
   double *rd;

   if (r <= 1) {
      rd = (double *) t;
      RELMEM ( doub_mh, (char *) (rd+ml) );
   } else {
      for (i = mh; i >= ml; i--)
	 free_dtensor( (double **) (t[i]), r-1, nl, nh, nl, nh); /* !!! */
      RELMEM( doub_mh, (char *) (t+ml) );
   }
}



#ifndef ANSI
void ptens (t, r, m, n)
/*********************/

double **t;
int    r, n, m;
#else
void ptens (double **t,int r,int m,int n)
/****************************************/
#endif
/* Prints a r-index tensor of size m*n^(r-1). */
/* Range [1...m] [1...n] ... [1...n] */

{
    int j;
   double *d;


   if (r <= 1) {
      d = (double *) t;
      for (j=1; j<= n; j++)  printf("%le ", d[j]);
      printf("\n");
   } else {
      for (j=1; j<= m; j++) {
	 d = t[j];
	 ptens( (double **) d, r-1, n, n);
	 if ((r>=3) && (j < m)) printf("\n");
      }
   }

}


#ifndef ANSI
void etens (t, r)
/***************/

double **t;
int r;
#else
void etens (double **t,int r)
/*****************************/
#endif
/* Prints a r-index tensor element of a tensor size m*n^(r-1). */
/* Depends on calling frame */
/* Range [1...m] [1...n] ... [1...n] */

{
   int j;
   int *ri,i;
   double *d;

   ri = &r;		   /* Sic! */
   printf("Tens.elem.");
   for (j=1; j < r; j++) {
      i = ri[j];	   /* !!!! */
      printf("[%d]",i);
      t = (double **) t[i];
   }
   d = (double *) t;
   i = ri[r];		   /* !!!! */
   printf("[%d]= %le ", i, d[i]);
   printf("\n");

}
