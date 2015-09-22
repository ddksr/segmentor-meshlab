// sphere.H

#ifndef LIBSEGMENTOR_SPHERE 
#define LIBSEGMENTOR_SPHERE 1

#include "model.h"

#include <iostream>
#include <stdio.h>
#include <math.h>

class sphere: public model
{ double a,b,c,radius;
  double data[20];

  double sqr(double x) { return(x*x); }
  double euclidian(struct point& p) { return(p.x*p.x+p.y*p.y+p.z*p.z); }

public:
  sphere() {}
  sphere(region& r);
  sphere(FILE *f);
  model* new_model(region &r) { return(new sphere(r)); }
  model* improve(region& orig, region& added);
  model* improve(model *m, region& orig, region& added) 
  { if (m) m = m;  // to avoid stupid compiler warnings
    return(improve(orig,added)); 
    }
  double abs_signed_distance(struct point& p) { return(sqrt(sqr(a-p.x)+sqr(b-p.y)+sqr(c-p.z)) - radius); }
  void draw()
  { switch(draw_type)
    { case M_OPENGL: 
        drawGL();
        break;
      }       
    }
  void drawGL();
  void print();
  void fprint(FILE *f);
  void save(FILE *f);
  int no_param() { return(4); }
  void get_normal_vector(double &x, double &y, double &z);
  MODELTYPE what_model() { return(CSPHERE); }
  void rif_write(FILE *f);
  void parameters(char **name, double *value);
  void set_parameters(double *value);
  struct point project(struct point &p);

};

#endif
