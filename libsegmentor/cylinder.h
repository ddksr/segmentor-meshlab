// cylinder.H 


#ifndef LIBSEGMENTOR_CYLINDER
#define LIBSEGMENTOR_CYLINDER 1

#include "model.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

class cylinder : public model 
{  
  int rep; /* == 0 if general quadric rep, == 1 if cylinder rep  == 2 Gabor rep*/

 /* In Gabor rep use a[1 - 3] for polar axes, 
      a[1]  = rho, a[2] =phi, a[3] = theta.
     a[4] for alpha and a[5] for k 
     a[11 - 13] for local Origin */


  /* In Gabor rep use a[1 - 3] for polar axes, 
     a[1]  = rho, a[2] =phi, a[3] = theta.
     a[4] for alpha and a[5] for k 
     a[11 - 13] for local Origin */
 
  /* general quadric rep, a[0] not used!! */
  /* a[1  - 10] quad params, a[11 - 13] c fo g */
  
  // Added by Bojan Kverh 19.2.1998
  // For accurate drawing

  

  struct point norm, disp;
  double r;

public:

  double a[14];
  double zmin,zmax;

  static double max_radius;  

  cylinder() {}
  cylinder(region &reg);
  cylinder(region &reg, cylinder &init);
  cylinder(FILE *f);
  
  model* new_model(region& r) { return(new cylinder(r)); }
  model* improve(region& orig, region& added);
  model* improve(model *m, region& orig, region& added);
  
  double abs_signed_distance(struct point &p); 

  void adjust_z(region& reg);

  void print();
  void fprint(FILE *f);
  void save(FILE *f);
  int no_param() {return(5);};

  void get_normal_vector(double &x, double &y, double &z);
  MODELTYPE what_model() { return(CCYLINDER);}
  void parameters(char **name, double *value);
  void set_parameters(double *value);
  struct point project(struct point& p);
  virtual void rotate(symatrix& A);
  virtual void translate(vector3D& v);
  virtual int model_ok();
};


#endif


