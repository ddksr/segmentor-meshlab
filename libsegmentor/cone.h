#ifndef LIBSEGMENTOR_CONE
#define LIBSEGMENTOR_CONE


#include "matrix.h"
#include "model.h"
#include "region.h"
#include "sline.h"

#include <iostream>

class cone: public model {

public:

  static double max_radius;

  double a[14],zmin, zmax;  
  /* In Gabor rep use a[1 - 3] for polar axes, 
     a[1]  = rho, a[2] =phi, a[3] = theta.
     a[4] for sigma and a[5]  for  tau and a[6] for k 
     a[11 - 13] for local Origin */
 
   
  
   
  
   hmatrix g_from_l;  // homogenuos matrix representation
  hmatrix l_from_g;


  sline sl;
  
  cone() {}
  cone(region &r);
  cone(region &r, cone &sp);
  cone(FILE *f);
  virtual ~cone() { }

  model* new_model(region& r) { return(new cone(r)); }
  model* improve(region& orig, region& added);
  model* improve(model *m, region& orig, region& added);
  

 int cone_fit(region &reg, double **list_x, double *list_y, double **norm, int no, int xspan, int yspan, double *param, double *chisq);

// int conedx(double, double, double, double*, double*, double *, int);

  
  double abs_signed_distance(struct point& ); 
  int no_param() {return(6); }

  void draw();
  void drawGL();
  void adjust_z(region &reg);
 
  void print();
  void fprint(FILE *f);
  void save(FILE *f);
  void get_normal_vector(double &x, double &y, double &z);
  MODELTYPE what_model() { return(CCONE); }
  void parameters(char **name, double *value);
  void set_parameters(double *value);
  struct point project(struct point &p);
  virtual void rotate(symatrix& A);
  virtual void translate(vector3D& v);
  virtual int model_ok();
};


#endif
