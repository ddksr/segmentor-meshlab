// torus.H 


#ifndef LIBSEGMENTOR_TORUS
#define LIBSEGMENTOR_TORUS   

#include "matrix.h"
#include "model.h"
#include "region.h"
#include "sline.h"


class torus: public model {


  // double m_error;  // provide maximum point distance and description error
  // double m_dist;   // for all objects of the class

  double a[14];  
  /* In Gabor rep 
     a[1]  = rho, a[2] =phi, a[3] = theta.
     a[4] for sigma and a[5]  for  tau and a[6] for k 
     a[7] for s
     a[8] for epsilon == -1 for lemon torus == +1 for apple torus
     a[11 - 13] for local Origin */
 
   
public:  
   
  
   hmatrix g_from_l;  // homogenuos matrix representation
  hmatrix l_from_g;

  //void torus_draw(int no_hidding = 0) const;

  sline sl;
  
  torus() {}  
  torus(region &r);
  torus(region &r, torus &sp);
  torus(FILE *f);
  virtual ~torus() { }

  model* new_model(region& r) { return(new torus(r)); }
  model* improve(region& orig, region& added);
  model* improve(model *m, region& orig, region& added);

 int torus_fit(double **list_x, double *list_y, double **norm, int no, int xspan, int yspan, double *param, double *chisq, region &reg);

  // int torusdx(double, double, double, double*, double*, double *, int);

  
  int torus_ok();

  double abs_signed_distance(struct point& );
  
  int no_param() {return(7); }

  void draw();
  void drawGL();

  void print();
  void fprint(FILE *f);
  void save(FILE *f);
  void get_normal_vector(double &x, double &y, double &z) {}
  MODELTYPE what_model() { return(CTORUS); }
  void parameters(char **name, double *value);
  void set_parameters(double *value);
  struct point project(struct point &p);
};


#endif

