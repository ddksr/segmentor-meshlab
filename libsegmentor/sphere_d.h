// sphere_d.H

#ifndef LIBSEGMENTOR_SPHERE_DESCRIPTION
#define LIBSEGMENTOR_SPHERE_DESCRIPTION 1

#include "description.h"
#include "sphere.h"

#include <stdio.h>

class sphere_d : public description
{
public:
  static double m_dist,m_err;
  
  sphere_d(region &r);
  sphere_d(FILE *f,image *im, image *norm, int camera_set);
  description* new_description(region& r, model *m = NULL) 
  { if (m) m = m;   // to avoid stupid compiler warnings
    return(new sphere_d(r)); 
    }
  double max_point_distance() { return(m_dist); }
  double max_error() { return(m_err); }
  void set_max_point_distance(double d) { m_dist = d; }
  void set_max_error(double d) { m_err = d; }
  int compatible(int i, int j);
  };

#endif
