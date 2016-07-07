// sq_d.H

#ifndef LIBSEGMENTOR_TSQ_DESCRIPTION
#define LIBSEGMENTOR_TSQ_DESCRIPTION 1

#include "description.h"

#include <stdio.h>

class tsq_d : public description
{ 
public:
  static double m_dist,m_err;

  tsq_d(region& r, model *m);
  tsq_d(FILE *f,image *im, image *norm, int camera_set);
  description* new_description(region& r, model *m = NULL) { return(new tsq_d(r,m)); }
  double max_point_distance() { return(m_dist); }
  double max_error() { return(m_err); }
  void set_max_point_distance(double d) { m_dist = d; }
  void set_max_error(double d) { m_err = d; }
  };

#endif
