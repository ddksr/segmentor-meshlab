// cone_d.C

#include "cone_d.h"
#include "cone.h"

cone_d::cone_d(region& r, model *m) : description(r)
{ mregion = new region(r);
  if (m != NULL) mmodel = new cone(r,*((cone *)m));
    else mmodel = new cone(r);
  }
  
cone_d::cone_d(FILE *f,image *im, image *norm, int camera_set) : description(f,im,norm, camera_set)
{ mmodel = new cone(f);
  }  

double cone_d::m_dist = 0.0;
double cone_d::m_err = 0.0;
