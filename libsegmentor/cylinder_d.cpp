// cylinder_d.C

#include "cylinder_d.h"
#include "cylinder.h"

cylinder_d::cylinder_d(region& r, model *m) : description(r)
{ mregion = new region(r);
  if (m != NULL) mmodel = new cylinder(r,*((cylinder *)m));
   else mmodel = new cylinder(r);
  }
  
cylinder_d::cylinder_d(FILE *f,image *im, image *norm, int camera_set) : description(f,im,norm, camera_set)
{ mmodel = new cylinder(f);
  }  

double cylinder_d::m_dist = 0.0;
double cylinder_d::m_err = 0.0;
