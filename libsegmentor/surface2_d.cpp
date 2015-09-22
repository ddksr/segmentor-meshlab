// surface2_d.C

#include "surface2_d.h"
#include "surface2.h"

surface2_d::surface2_d(region& r) : description(r)
{ mregion = new region(r);
  mmodel = new surface2(r);
  }
  
surface2_d::surface2_d(FILE *f,image *im,image *norm, int camera_set) : description(f,im,norm, camera_set)
{ mmodel = new surface2(f);
  }  
  
int surface2_d::compatible(int i, int j)
{ int k;
//  double x,y,z;
//  struct point p;
  
  k = (mmodel->distance(mregion->get_point(i,j)) <= max_point_distance());
  /*  if (k && normals != NULL)
  { mmodel -> get_normal_vector(x,y,z);
    p = normals->pixel(i,j);
    } */
  return(k);  
  }

double surface2_d::m_dist = 0.0;
double surface2_d::m_err = 0.0;



