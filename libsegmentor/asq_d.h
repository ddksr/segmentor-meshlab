// sq_d.H

#ifndef LIBSEGMENTOR_ASQ_DESCRIPTION
#define LIBSEGMENTOR_ASQ_DESCRIPTION 1

#include "description.h"
#include "sq.h"

#include <stdio.h>

class asq_d : public description
{ 
public:
  static double m_dist,m_err;

  asq_d(region& r, model *m);
  asq_d(FILE *f,image *im, image *norm, int camera_set);
  description* new_description(region& r, model *m = NULL) { return(new asq_d(r,m)); }
  double max_point_distance() { return(m_dist); }
  double max_error() { return(m_err); }
  void set_max_point_distance(double d) { m_dist = d; }
  void set_max_error(double d) { m_err = d; }
  };

#endif
