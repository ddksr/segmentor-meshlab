// torus_d.C

#include "torus_d.h"
#include "torus.h"

torus_d::torus_d(region& r, model *m) : description(r)
{ mregion = new region(r);
  if (m != NULL) mmodel = new torus(r,*((torus *)m));
    else mmodel = new torus(r);
  }
  
torus_d::torus_d(FILE *f,image *im, image *norm, int camera_set) : description(f,im,norm, camera_set)
{ mmodel = new torus(f);
  }  

double torus_d::m_dist = 0.0;
double torus_d::m_err = 0.0;
