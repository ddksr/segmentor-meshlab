// sq_d.C

#include "asq_d.h"
#include "sq.h"

asq_d::asq_d(region& r, model *m) : description(r)
{ mregion = new region(r);
  if (m != NULL) mmodel = new sq((sq *)m,r);
    else mmodel = new sq(r);
  }
  
asq_d::asq_d(FILE *f,image *im, image *norm, int camera_set) : description(f,im,norm,camera_set)
{ mmodel = new sq(f);
  }  

double asq_d::m_dist = 0.0;
double asq_d::m_err = 0.0;
