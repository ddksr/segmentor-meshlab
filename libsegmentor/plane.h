// plane.H

#ifndef LIBSEGMENTOR_PLANE
#define LIBSEGMENTOR_PLANE 1

#include "model.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

class plane: public model
{ 
public:
  double a, b, c, d;              // parameters, describing a plane
  double data[10];                  // data for fast recovery;
  symatrix *LCS;

public:
  plane() {}
  plane(double a1, double b1, double c1, double d1) 
  {  a = a1;
      b = b1;
      c = c1;
	d = d1;
      }
  plane(region& r);
  plane(FILE *f);
  model* new_model(region& r) { return(new plane(r)); }
  model* improve(region& orig, region& added);
  model* improve(model *m, region& orig, region& added) 
  { if (m) m = m; // to avoid stupid compiler warnings.
    return(improve(orig,added)); 
    }
  double abs_signed_distance(struct point& p) { return(a*p.x+b*p.y+c*p.z+d); }
  void draw() {}
  void drawGL() {}
  void print(); 
  void fprint(FILE *f);
  void save(FILE *f);
  int no_param() { return(3); }
  void get_normal_vector(double &x, double &y, double &z);
  double model_distance();
  MODELTYPE what_model() { return(CPLANE); }
  void rif_write(FILE *f);
  void parameters(char **name, double *value);
  void set_parameters(double *value);
  struct point project(struct point &p);
  void  calc_lcs();
  vector3D point_to_lcs(struct point& p);
  void rotate(symatrix& A);
  void translate(vector3D& v);
  };
  
#endif
