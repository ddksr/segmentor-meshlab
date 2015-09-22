/* rcad_circurv.h */

#ifndef LIBSEGMENTOR_RCAD_CIRCURV_H
#define LIBSEGMENTOR_RCAD_CIRCURV_H

#include "rcad_SurfaceData.h"


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
	     double *dir2);

  /* Neglects 'centre'. Curvature is computed in sd[0]. Orders points
     around sd[0] looking from 'normal' and then the circles
     Pi P0 P(i+cirn) are considered where
                      cirn = (npnts-1)/2,
     the number of circles. 'npnts' must be at least 7. */
#endif

