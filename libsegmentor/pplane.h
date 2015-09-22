// Projection plane class
// written by Bojan Kverh 20.5.98

#ifndef LIBSEGMENTOR_PPLANE 
#define LIBSEGMENTOR_PPLANE 1

#include "matrix.h"
#include "region.h"

class pplane
{ vector3D n,r1,r2;

public:
  pplane(region &reg, int *ind);
  };

double point_distance(double x1, double y1, double x2, double y2);


#endif
