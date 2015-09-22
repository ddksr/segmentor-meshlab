// planar_d.C

#include "planar_d.h"
#include "plane.h"

planar_d::planar_d(region& r) : description(r)
{ mregion = new region(r);
  mmodel = new plane(r);
  }
  
planar_d::planar_d(FILE *f,image *im, image *norm, int camera_set) : description(f,im,norm, camera_set)
{ mmodel = new plane(f);
  }  
  
int planar_d::compatible(int i, int j)
{ int k;
  double x,y,z;
  struct point p;
  
  k = (mmodel->distance(mregion->get_point(i,j)) <= max_point_distance());
  if (k && normals != NULL)
  { mmodel -> get_normal_vector(x,y,z);
    p = normals->pixel(i,j);
    k &= (fabs(p.x*x + p.y*y + p.z * z) >= normal_compatibility);
    }
  return(k);  
  }

double planar_d::normal_compatibility = 0.0;
double planar_d::m_dist = 0.0;
double planar_d::m_err = 0.0;
