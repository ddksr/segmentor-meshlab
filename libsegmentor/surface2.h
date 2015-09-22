// surface2.H

#ifndef LIBSEGMENTOR_SURFACE2
#define LIBSEGMENTOR_SURFACE2 1

#include "model.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

class surface2: public model
{ double a, b, c, d, e, f;              // parameters, describing 2nd order surface
  double data[25];                      // data for fast recovery;

public:    
  surface2() {}
  surface2(region& r);
  surface2(FILE *ff);
  model* new_model(region& r) { return(new surface2(r)); }
  model* improve(region& orig, region& added);
  model* improve(model *m, region& orig, region& added) 
  { m = m; // to avoid stupid compiler warnings
    return(improve(orig,added)); 
    }
  double abs_signed_distance(struct point& p) { return(a*p.x*p.x+b*p.x*p.y+c*p.y*p.y+d*p.x+e*p.y+f-p.z); }
  void draw() {}
  void drawGL() {}
  void print(); 
  void fprint(FILE *ff);
  void save(FILE *ff);
  int no_param() { return(6); }
  void get_normal_vector(double &x, double &y, double &z) 
  { if (x > y && y > z) x = x; // to avoid stupid compiler warnings.
    }
  MODELTYPE what_model() { return(CSURFACE2); }
  void parameters(char **name, double *value);
  void set_parameters(double *value);

  };
  
#endif
