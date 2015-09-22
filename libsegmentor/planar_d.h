// planar_d.H

#ifndef LIBSEGMENTOR_PLANAR_DESCRIPTION
#define LIBSEGMENTOR_PLANAR_DESCRIPTION 1

#include "description.h"
#include "plane.h"
#include <stdio.h>

class planar_d : public description
{ 
public:
  static double m_dist,m_err;
  static double normal_compatibility;

  planar_d(region& r);
  planar_d(FILE *f,image *im,image *norm, int camera_set);
  description* new_description(region& r, model *m = NULL) 
  { m = m;     // to avoid stupid compiler warnings.
    return(new planar_d(r)); 
    }
  double max_point_distance() { return(m_dist); }
  double max_error() { return(m_err); }
  void set_max_point_distance(double d) { m_dist = d; }
  void set_max_error(double d) { m_err = d; }
  int compatible(int i, int j);
  };

#endif
