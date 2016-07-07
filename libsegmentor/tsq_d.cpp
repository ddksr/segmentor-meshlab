// sq_d.C

#include "tsq_d.h"
#include "tsq.h"

tsq_d::tsq_d(region& r, model *m) : description(r)
{ mregion = new region(r);
  if (m != NULL) mmodel = new tsq((tsq *)m,r);
    else mmodel = new tsq(r);
  }
  
tsq_d::tsq_d(FILE *f,image *im, image *norm, int camera_set) : description(f,im,norm,camera_set)
{
  // legacy
}  

double tsq_d::m_dist = 0.0;
double tsq_d::m_err = 0.0;
