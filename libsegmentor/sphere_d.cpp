// sphere_d.C

#include "sphere_d.h"
#include "sphere.h"

sphere_d::sphere_d(region &r) : description(r)
{ mregion = new region(r);
  mmodel = new sphere(r);
  }

sphere_d::sphere_d(FILE *f,image *im, image *norm, int camera_set) : description(f,im,norm,camera_set)
{ mmodel = new sphere(f);
  }

int sphere_d::compatible(int i, int j)
{ int k;
  k = (mmodel -> distance(mregion->get_point(i,j)) <= max_point_distance());
  return(k);
  }

double sphere_d::m_dist = 0.0;
double sphere_d::m_err = 0.0;


