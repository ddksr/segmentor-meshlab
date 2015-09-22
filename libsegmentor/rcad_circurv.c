/* rcad_circurv.c */

#define ANSI

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "ndm.h"
#include "ndm_svd.h"
#include "rcad_ndm_ut.h"

#include "rcad_vector2.h"
#include "rcad_circurv.h"


static double
ComputeCirCurv( double *a,     /* P1 - P0 vector */
		double *b,     /* P2 - P0 vector */
		double *n,     /* Out: unit normal vector of the circle
				  plane or zero */
		double *t)     /* Out: tangent of circle in P0 */
  /* It returns the curvature (non-negative) of the circle */
  /* All indices go from 0 up to 2. */
{

  double c[3];
  double laxb;
  double la,lb;
  double k;

  v_mult(n, a, b);
  laxb = s_mult(n,n);
  la   = s_mult(a,a);
  lb   = s_mult(b,b);

  if (laxb < DBL_EPSILON)
  {
    n[0] = n[1] = n[2] = 0.0;
    k = 0.0;
  } else {
    laxb = sqrt(laxb);
    mult(n, 1.0/laxb, n);                        /* a x b /|a x b| */
    sub(c,b,a);                                 /* P2 - P1: third side of the
						   triangle */
    k = 2.0*laxb/(sqrt(la)*sqrt(lb)*vabs(c));   /*  2*|a x b|/(|a||b||c|) */
  }

  mult(c, lb, a);
  mult(t, la, b);
  sub(t, c,  t);                                /* |b|^2*a - |a|^2*b */
    
  laxb = vabs(t);
  if (laxb < DBL_EPSILON)
  {
    /* It may occur only if  a x b = 0 */
    t[0] = t[1] = t[2] = 0.0;
  } else {
    mult( t, 1.0/laxb, t);
  }

  return k;

}


int
rcad_circurv(int        npnts,      /* Number of points */
	     int *np,         /* Point indices in the sd array
				        (npnts pieces) */
	     double *centre,  /* Point where the curvature is to be
					 computed */
	     double *normal,        /* Normal vector guess (input/output)
					 (e.g. for orientation) */
	     rcad_SurfaceData *sd,
	                            /* Array containing point data */
	     double *k1,            /* Output principal curvatures */
	     double *k2,
	     double *dir1,          /* Output principal directions */
	     double *dir2)
{
  /* Neglects 'centre'. Curvature is computed in sd[0]. Orders points
     around sd[0] looking from 'normal' and then the circles
     Pi P0 P(i+cirn) are considered where
                      cirn = (npnts-1)/2,
     the number of circles. 'npnts' must be at least 7. */
  /* Vector package indexes from 0. ndm (Numerical recipes) indexes
     from 1. */

  int i;
  int j;
  int k,ok_sing;
  int cirn;
  int error=0;
  double NormRefDir[3];
  const double conditionRatio= 1.0e-6;
  double  p0[3];
  double  a[3];
  double  b[3];
/*  double  c[3]; */
  double  xu[3];
  double  yu[3];
  double  sing_vals[3];
  double  act[3];
  double  *w;
  char    *svdbuf=0;
  double *right=0;
  double *cir_curvs=0;
  double **cir_norms=0;
  double **cir_tans=0;
  double **cir_aplusb=0;
  double **PlNormMatr=0;
  double **CurvMatr=0;
  double **svdV=0;
  double *pangle=0;
  int    *pind=0;
  int      ind;
  double singular;
  double MaxSingVals;
  double h;
  double cth;
  double sth;
  double scale;
  const double mincurv2=1.0e-8; /* Cut for 0 curvature scaling */

  if (centre) centre = centre; /* to avoid stupid compiler warnings */
  cirn=(npnts-1)/2;   /* Number of circles */

  *k1 = *k2 = 0.0;


  dir1[0] = dir1[1] = dir1[2] = 0.0;
  dir2[0] = dir2[1] = dir2[2] = 0.0;


  if (cirn < 3)
  {
    error = -22;    /* Few (< 3) circles found */
    /* Normal vector remains as computed previously */
  } else {
    
    NormRefDir[0] = normal[0];
    NormRefDir[1] = normal[1];
    NormRefDir[2] = normal[2];

    init_svdcmp( cirn, 1, &svdbuf);

    cir_curvs= dvector(1,cirn);
    right    = dvector(1,cirn);
    cir_norms = dmatrix(1,cirn,1,3);
    cir_tans  = dmatrix(1,cirn,1,3);
    cir_aplusb = dmatrix(1,cirn,1,3);
    PlNormMatr = dmatrix(1,3,1,3);
    svdV       = dmatrix(1,3,1,3);
    CurvMatr   = dmatrix(1,cirn,1,3);

    assign(p0, sd[np[0]].p); 

    /* Ordering points by angle around p0 looking from normal */

    makeCoordinateSystem( NormRefDir, xu, yu);

    pangle     = dvector(1,npnts-1);
    pind       = (int *) malloc(sizeof(int)*npnts);
         /* The 0th is not used */

    for (i=1; i < npnts; ++i)
    {
      pind[i] =
      ind = np[i];
      sub(a,sd[ind].p, p0);
      h = atan2(s_mult(a,yu),s_mult(a,xu));
      for (j=1; j < i; ++j)
      {
	if (h < pangle[j])
	{
	  for (k=i-1; k >= j; --k)
	  {
	    pind[k+1]   = pind[k];
	    pangle[k+1] = pangle[k];
	  }
	  break;
	}
      }
      pind[j]   = ind;
      pangle[j] = h; 
    }

    /* Computation of circle data */

    for (i=1; i <= cirn; ++i)
      {
	sub(a, sd[pind[i]].p, p0);
	sub(b, sd[pind[i+cirn]].p, p0);
	add(cir_aplusb[i]+1, a, b);
	cir_curvs[i] = ComputeCirCurv(a,b, cir_norms[i]+1, cir_tans[i]+1);
      }

    PlNormMatr[1][1] = 0.0;
    PlNormMatr[2][2] = 0.0;
    PlNormMatr[3][3] = 0.0;
    PlNormMatr[1][2] = 0.0;
    PlNormMatr[1][3] = 0.0;
    PlNormMatr[2][3] = 0.0;

    for( i = 1; i <= cirn; ++i)
      {
	w = cir_tans[i];
	PlNormMatr[1][1] += w[1]*w[1];
	PlNormMatr[2][2] += w[2]*w[2];
	PlNormMatr[3][3] += w[3]*w[3];
	PlNormMatr[1][2] += w[1]*w[2];
	PlNormMatr[1][3] += w[1]*w[3];
	PlNormMatr[2][3] += w[2]*w[3];
      }

    PlNormMatr[2][1] = PlNormMatr[1][2];
    PlNormMatr[3][1] = PlNormMatr[1][3];
    PlNormMatr[3][2] = PlNormMatr[2][3];

    ok_sing = 1;

    for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++)
        if (PlNormMatr[i][j] + 0.0 != PlNormMatr[i][j]) /* detect NaN */
	{ fprintf(stderr,"\nNaN detected in rcad_circurv\n");
          ok_sing = 0;
	  }

    if (ok_sing) 
    
    { svdcmp( PlNormMatr, 3, 3, sing_vals-1, svdV, &svdbuf);

      singular=MaxSingVals=fabs(sing_vals[2]);

      j = 2;

      for( i = 0; i < 2; ++i)
      {
	if (fabs(sing_vals[i]) < singular)
	  {
	    singular = fabs(sing_vals[i]);
	    j = i;
	  }
	if (fabs(sing_vals[i]) > MaxSingVals)
	  {
	    MaxSingVals = fabs(sing_vals[i]);
	  }
      }


    if (MaxSingVals < DBL_EPSILON)
    {
      error = -24;
    } else {


      normal[0] = svdV[1][j+1];
      normal[1] = svdV[2][j+1];
      normal[2] = svdV[3][j+1];

      if (s_mult(NormRefDir, normal) < 0.0)
	{
	  mult(normal, -1.0, normal);
	}

      for (i=1; i <= cirn; ++i)
	{
	  mult(act, s_mult(normal, cir_tans[i]+1), cir_tans[i]+1);
	  sub(xu, cir_tans[i]+1, act);
	  h = vabs(xu);
	  if (h > DBL_EPSILON)
	    {
	      mult(xu, 1.0/h, xu);
	      v_mult(yu, normal, xu);
	      break;
	    }
	}

      if (i >= cirn) {
	error = -23;          /* All tangents vanish */
      } else {

	for( j = 1; j <= cirn;  j++)
	  {
	    mult(act, s_mult(normal, cir_tans[j]+1), cir_tans[j]+1);
	    sub(act, cir_tans[j]+1, act);
	    h = vabs(act);
	    if (h > DBL_EPSILON)
	      {
		mult(act, 1.0/h, act);
		h = s_mult(cir_norms[j]+1,normal); /* sine */
		right[j] = scale = sqrt(1.0-h*h)*cir_curvs[j];
		scale = 1.0/sqrt(scale*scale + mincurv2);
		right[j] *= scale;
		if (s_mult(cir_aplusb[j]+1,normal) < 0.0)
		  {
		    right[j] *= -1.0;
		  }
		h = atan2( s_mult(act, yu), s_mult(act, xu));
	   
		CurvMatr[j][2] = -scale*cos(2.0*h);
		CurvMatr[j][3] = -scale*sin(2.0*h);
	      } else {
		right[j] = scale = 0.0;
		CurvMatr[j][2] = CurvMatr[j][3] = 0.0;
	      }
	    CurvMatr[j][1] = -scale;
	  }
    
       for (i = 1; i <= 3; i++)
         for (j = 1; j <= 3; j++)
           if (CurvMatr[i][j] + 0.0 != CurvMatr[i][j]) /* detect NaN */
	   { fprintf(stderr,"\nNaN detected in rcad_circurv\n");
              ok_sing = 0;
              error = -28;
	      } 

          if (ok_sing) svdcmp( CurvMatr, cirn, 3, sing_vals-1, svdV, &svdbuf);

	  MaxSingVals = fabs(sing_vals[0]);

	  for( j = 1; j < 3; ++j)
	    if (fabs(sing_vals[j]) > MaxSingVals)
	      {
	        MaxSingVals = fabs(sing_vals[j]);
	      }

  	  for( j = 1; j <= 3; ++j)
	    if (fabs(sing_vals[j]) < conditionRatio*MaxSingVals)
	      sing_vals[j]=0.0;

	  svbksb( CurvMatr, sing_vals-1, svdV, cirn, 3, right, act-1, &svdbuf);

	  h = sqrt(act[1]*act[1] + act[2]*act[2]);

	  *k1 = act[0] + h;
	  *k2 = act[0] - h;
		   
	  h = 0.5*atan2( act[2], act[1]);
	  cth = cos (h);
	  sth = sin (h);
      
	  mult(dir1, cth, xu);
	  mult(a, sth, yu);
	  add(dir1,dir1,a);
      
	  mult(dir2, -sth, xu);
	  mult(a, cth, yu);
	  add(dir2,dir2,a);
	
      }
    }
    init_svdcmp( cirn, 0, &svdbuf);
        
    } else error = -27;
  }        

    free_dvector(pangle, 1,npnts-1);
    free(pind);

    free_dvector(cir_curvs, 1, cirn);
    free_dmatrix(cir_norms,1,cirn,1,3);
    free_dmatrix(cir_tans, 1, cirn, 1, 3);
    free_dmatrix(cir_aplusb, 1, cirn, 1, 3);
    free_dmatrix(PlNormMatr, 1,3,1,3);
    free_dmatrix(svdV, 1,3,1,3);
    free_dmatrix(CurvMatr,1,cirn,1,3);
    free_dvector(right, 1, cirn);

    

   return error;
}



