// sq.H

#ifndef LIBSEGMENTOR_ASUPERQUADRIC
#define LIBSEGMENTOR_ASUPERQUADRIC 1

#include "region.h"
#include "matrix.h"
#include "model.h"

#include <stdio.h>

class asq : public model
{ 
  
  double f(double x, double y, double z) const;

public:
  double a1,a2,a3,e1,e2,px,py,pz,phi,theta,psi,kxy,kf;    // SQ parameters
  hmatrix l_from_g,g_from_l;
  
  
  asq() {}
  asq(region& r);
  asq(asq *m, region& r);
  virtual ~asq() { }
  model* new_model(region& r) { return( new asq(r) ); }
  model* improve(region& orig, region& added);
  model* improve(model *m, region& orig, region& added);
  double abs_signed_distance(struct point& p);
  void print(); 

  int no_param() { return(13); }
  void get_normal_vector(double &x, double &y, double &z) 
  { if (x > y && y > z) x = x;  // to avoid stupid compiler messages
    }
  
  vector normal(double, double) const;
  vector r(double, double) const;
  vector transform(vector& v);

  MODELTYPE what_model() { return(CASQ); }  

  void parameters(char **name, double *value);
  void set_parameters(double *value);

  double map_eta(double);
  double map_omega(double);

  // legacy
  void fprint(FILE *f) {}
  void save(FILE *f) {}
  };


#endif
